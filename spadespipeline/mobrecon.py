#!/usr/bin/env python3
from accessoryFunctions.accessoryFunctions import GenObject, MetadataObject, run_subprocess, strainer, write_to_logfile, make_path, SetupLogging
from argparse import ArgumentParser
from click import progressbar
import multiprocessing
import pandas as pd
import logging
import os

__author__ = 'adamkoziol'


class MobRecon(object):

    def mob_recon(self):
        self.recon()
        self.read_tsv()
        self.summary_reporter()
        self.reporter()

    def recon(self):
        """
        Prep
        """
        logging.info('Running MOB-recon on Assemblies')
        with progressbar(self.metadata) as bar:
            for sample in bar:
                # Create and populate the mob_recon genobject
                setattr(sample, self.analysistype, GenObject())
                sample[self.analysistype].outputdir = os.path.join(sample.general.outputdirectory, self.analysistype)
                sample[self.analysistype].contig_report = os.path.join(sample[self.analysistype].outputdir,
                                                                       'contig_report.txt')
                sample[self.analysistype].report_dict = dict()
                sample[self.analysistype].logout = os.path.join(sample[self.analysistype].outputdir, 'out')
                sample[self.analysistype].logerr = os.path.join(sample[self.analysistype].outputdir, 'err')
                make_path(sample[self.analysistype].outputdir)
                if sample.general.bestassemblyfile != 'NA':
                    # Create the system call
                    sample.commands.mobrecon = 'mob_recon -i {fasta} -o {outdir} --run_typer -n {threads}'\
                        .format(fasta=sample.general.bestassemblyfile,
                                outdir=sample[self.analysistype].outputdir,
                                threads=self.threads)
                    # Ensure that the report doesn't already exist
                    if not os.path.isfile(sample[self.analysistype].contig_report):
                        # Run the analyses
                        out, err = run_subprocess(sample.commands.mobrecon)
                        # Write the outputs to the log file
                        write_to_logfile(out=sample.commands.mobrecon,
                                         err=sample.commands.mobrecon,
                                         logfile=self.logfile,
                                         samplelog=sample.general.logout,
                                         sampleerr=sample.general.logerr,
                                         analysislog=sample[self.analysistype].logout,
                                         analysiserr=sample[self.analysistype].logerr)
                        write_to_logfile(out=out,
                                         err=err,
                                         logfile=self.logfile,
                                         samplelog=sample.general.logout,
                                         sampleerr=sample.general.logerr,
                                         analysislog=sample[self.analysistype].logout,
                                         analysiserr=sample[self.analysistype].logerr)

    def read_tsv(self):
        """
        Read in the .tsv contig report file with pandas, and create a dictionary of all the headers: values
        """
        logging.info('Parsing MOB-recon outputs')
        for sample in self.metadata:
            if os.path.isfile(sample[self.analysistype].contig_report):
                # Read in the .tsv file with pandas. Skip the comment lines
                df = pd.read_csv(sample[self.analysistype].contig_report, delimiter='\t')
                for header in df:
                    # Remove any unwanted whitespace
                    clean_header = header.lstrip().rstrip()
                    # primary_key is the primary key, and value is the value of the cell for that key + header combo
                    for primary_key, value in df[header].items():
                        # Update the dictionary with the new data
                        try:
                            sample[self.analysistype].report_dict[primary_key].update({clean_header: value})
                        # Create the nested dictionary if it hasn't been created yet
                        except KeyError:
                            sample[self.analysistype].report_dict[primary_key] = dict()
                            sample[self.analysistype].report_dict[primary_key].update({clean_header: value})

    def summary_reporter(self):
        """
        Parse individual MOB Recon reports into a summary report
        """
        logging.info('Creating MOB-recon summary report')
        with open(os.path.join(self.reportpath, 'mob_recon_summary.csv'), 'w') as summary:
            data = 'Strain,Location,Contig,Incompatibility,IncompatibilityAccession,RelaxaseType,' \
                   'MashNearestNeighbor,MashNeighborDistance\n'
            for sample in self.metadata:
                # Initialise a dictionary to store results for the COWBAT final report
                sample[self.analysistype].pipelineresults = dict()
                for primarykey, results in sample[self.analysistype].report_dict.items():
                    # Only process results if they are not calculated to be chromosomal
                    if results['cluster_id'] != 'chromosome':
                        data += ','.join(str(result).replace(',', ';') if str(result) != 'nan' else 'ND'
                                         for result in [
                                             sample.name,
                                             results['cluster_id'],
                                             results['contig_id'].split('|')[1],
                                             results['rep_type'],
                                             results['rep_type_accession'],
                                             results['relaxase_type'],
                                             results['mash_nearest_neighbor'],
                                             results['mash_neighbor_distance']]
                                         )
                        data += '\n'
                        # Add the calculated incompatibility to the pipeline results for use in the final COWBAT report
                        sample[self.analysistype].pipelineresults[results['cluster_id']] =  \
                            ';'.join(str(result).replace(',', ';') if str(result) != 'nan' else 'ND'
                                     for result in [
                                         results['rep_type']]
                                     )
            summary.write(data)

    def reporter(self):
        """
        Create a report combining results from resfinder_assembled and mob_recon_summary reports
        """
        logging.info('Creating AMR summary table from ResFinder and MOB-recon outputs')
        with open(os.path.join(self.reportpath, 'amr_summary.csv'), 'w') as amr:
            data = 'Strain,Gene,Allele,Resistance,PercentIdentity,Contig,Location,PlasmidIncompatibilitySets\n'
            for sample in self.metadata:
                # Initialise a dictionary to store a set of all the incompatibility types listed for a contig.
                # As the inc type will only be located on one of possibly several contigs associated with a predicted
                # plasmid, it is nice to know details about the plasmid
                inc_dict = dict()
                for primarykey, results in sample[self.analysistype].report_dict.items():
                    try:
                        inc = results['cluster_id']
                        # Convert the rep_type field (predicted incompatibilities) into a more a consistent
                        # format - pandas will call empty fields 'nan', which is a float
                        rep = str(results['rep_type']).replace(',', ';') if str(results['rep_type']) != 'nan' else 'ND'
                        # Add the incompatibility to the set
                        try:
                            inc_dict[inc].add(rep)
                        except KeyError:
                            inc_dict[inc] = set()
                            inc_dict[inc].add(rep)
                    except KeyError:
                        pass
                #
                for primarykey, results in sample[self.analysistype].report_dict.items():
                    try:
                        contig = results['contig_id'].split('|')[1]
                        # Use the list of results from the resfinder analyses
                        for amr_result in sample.resfinder_assembled.sampledata:
                            # Ensure that the current contig is the same as the one in the resfinder results
                            if contig in amr_result:
                                # Set up the output string
                                data += '{},'.format(sample.name)
                                # Add the
                                data += '{amr},{mob}\n'\
                                    .format(amr=','.join(str(res) if str(res) != 'nan' else 'ND' for res in
                                                         amr_result[0:4]),
                                            mob=','.join(str(res) if str(res) != 'nan' else 'ND' for res in
                                                         [contig, results['cluster_id'],
                                                          ';'.join(sorted(inc_dict[str(results['cluster_id'])]))
                                                          ]
                                                         )
                                            )
                    except KeyError:
                        pass
            amr.write(data)

    def __init__(self, metadata, analysistype, databasepath, threads, logfile, reportpath):
        self.metadata = metadata
        self.analysistype = analysistype
        self.databasepath = os.path.join(databasepath, analysistype)
        self.threads = threads
        self.logfile = logfile
        self.reportpath = reportpath


if __name__ == '__main__':

    def resfinder_extract(reportpath, metadata):
        # A dictionary to store the parsed excel file in a more readable format
        nesteddictionary = dict()
        resfinder_report = os.path.join(reportpath, 'resfinder_blastn.xlsx')
        assert os.path.isfile(resfinder_report), 'Missing ResFinder Report!'
        dictionary = pd.read_excel(resfinder_report).to_dict()
        for header in dictionary:
            # Sample is the primary key, and value is the value of the cell for that primary key + header combination
            for sample, value in dictionary[header].items():
                # Update the dictionary with the new data
                try:
                    nesteddictionary[sample].update({header: value})
                # Create the nested dictionary if it hasn't been created yet
                except KeyError:
                    nesteddictionary[sample] = dict()
                    nesteddictionary[sample].update({header: value})
        for sample in metadata:
            sample.resfinder_assembled = GenObject()
            sample.resfinder_assembled.sampledata = list()
            for line in nesteddictionary:
                sample_name = nesteddictionary[line]['Strain']
                if sample_name == sample.name:
                    sample.resfinder_assembled.sampledata.append(([nesteddictionary[line]['Gene'],
                                                                   nesteddictionary[line]['Allele'],
                                                                   nesteddictionary[line]['Resistance'],
                                                                   nesteddictionary[line]['PercentIdentity'],
                                                                   nesteddictionary[line]['PercentCovered'],
                                                                   nesteddictionary[line]['Contig']]))
        return metadata

    # Parser for arguments
    parser = ArgumentParser(description='Performing the typing component of the COWBAT pipeline on assemblies')
    parser.add_argument('-s', '--sequencepath',
                        required=True,
                        help='Path to folder containing sequencing reads')
    parser.add_argument('-r', '--referencefilepath',
                        required=True,
                        help='Provide the location of the folder containing the pipeline accessory files (reference '
                             'genomes, MLST data, etc.')
    SetupLogging()
    arguments = parser.parse_args()
    # Extract the list of strains in the sequence path, and create a metadata object with necessary values
    strains, metadata = strainer(arguments.sequencepath)
    report_path = os.path.join(arguments.sequencepath, 'reports')
    #
    metadata = resfinder_extract(reportpath=report_path,
                                 metadata=metadata)
    mob = MobRecon(metadata=metadata,
                   analysistype='mobrecon',
                   databasepath=arguments.referencefilepath,
                   threads= multiprocessing.cpu_count() - 1,
                   logfile=os.path.join(arguments.sequencepath, 'log'),
                   reportpath=report_path)
    mob.mob_recon()
