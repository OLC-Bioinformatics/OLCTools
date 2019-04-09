#!/usr/bin/env python3
from accessoryFunctions.accessoryFunctions import GenObject, MetadataObject, run_subprocess, strainer, write_to_logfile, make_path, SetupLogging
from argparse import ArgumentParser
from click import progressbar
import multiprocessing
from glob import glob
import pandas as pd
import logging
import csv
import os

__author__ = 'adamkoziol'


class MobRecon(object):

    def mob_recon(self):
        self.recon()
        self.read_tsv()
        self.summary_reporter()
        if self.matchtype == 'amrsummary':
            self.amrsummary()
        else:
            self.geneseekrsummary()

    def recon(self):
        """
        Prep
        """
        logging.info('Running MOB-recon on Assemblies')
        with progressbar(self.metadata) as bar:
            commands = list()
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
                    sample.commands.mobrecon = 'mob_recon -i {fasta} -o {outdir} --run_typer -n {threads} -d {db}'\
                        .format(fasta=sample.general.bestassemblyfile,
                                outdir=sample[self.analysistype].outputdir,
                                threads=self.threads,
                                db=os.path.join(self.databasepath, 'mob_suite'))
                    # Ensure that the report doesn't already exist
                    if not os.path.isfile(sample[self.analysistype].contig_report):
                        # Run the analyses
                        commands.append(sample.commands.mobrecon)
            p = multiprocessing.Pool(processes=self.threads)
            out_err = p.map(MobRecon.run_cmd, commands)
            p.close()
            p.join()
            # At this point, out_err has a list of tuples with out as index 0 in each tuple and
            # err at index 1 in each tuple. These will be in the same order as the samples, so retrieve them by index.
            index = 0
            for sample in bar:
                # Write the outputs to the log file
                out = out_err[index][0]
                err = out_err[index][1]
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
                index += 1

    @staticmethod
    def run_cmd(command):
        out, err = run_subprocess(command)
        return out, err

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

    def amrsummary(self):
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
                        # Unicycler gives contigs names such as: 3_length=187116_depth=1.60x_circular=true - test
                        # to see if the contig name looks unicycler-like, and set the name appropriately (in this
                        # case, it would be 3)
                        if contig.split('_')[1].startswith('length'):
                            contig = contig.split('_')[0]
                        # Use the list of results from the resfinder analyses
                        for amr_result in sample.resfinder_assembled.sampledata:
                            # Ensure that the current contig is the same as the one in the resfinder results. Ensure
                            # that the slice of the amr result is treated as a string. Unicycler contigs seem to be
                            # treated as integers
                            if contig == str(amr_result[-1]):
                                # Set up the output string
                                data += '{sn},'.format(sn=sample.name)
                                # Add the resistance and MOB recon outputs for the strain
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

    def geneseekrsummary(self):
        """
        Create a report combining GeneSeekr and MOB Recon outputs
        """
        logging.info('Creating predicted plasmid-borne gene summary table')
        with open(os.path.join(self.reportpath, 'plasmid_borne_summary.csv'), 'w') as pbs:
            data = 'Strain,Gene,PercentIdentity,Contig,Location,PlasmidIncompatibilitySets\n'
            for sample in self.metadata:
                # Create a flag to determine whether the strain name needs to be added to the data string if there
                # were no results
                result_bool = False
                # Initialise a dictionary to store a set of all the incompatibility types listed for a contig.
                # As the inc type will only be located on one of possibly several contigs associated with a predicted
                # plasmid, it is nice to know details about the plasmid
                inc_dict = dict()
                # Iterate through all the MOB recon outputs to populate the incompatibility set
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
                for primarykey, results in sample[self.analysistype].report_dict.items():
                    try:
                        contig = results['contig_id'].split('|')[1]
                        # Unicycler gives contigs names such as: 3_length=187116_depth=1.60x_circular=true - test
                        # to see if the contig name looks unicycler-like, and set the name appropriately (in this
                        # case, it would be 3)
                        if contig.split('_')[1].startswith('length'):
                            contig = contig.split('_')[0]
                        for gene, result_dict in sample.geneseekr_results.sampledata.items():
                            if contig == result_dict['query_id']:
                                percent_identity = result_dict['PercentIdentity']
                                # Set up the output string if the percent identity of the match is greater than the
                                # cutoff
                                if float(result_dict['PercentIdentity']) >= self.cutoff:
                                    # As there was at least a single gene passing the threshold, set the boolean to True
                                    result_bool = True
                                    data += '{sn},'.format(sn=sample.name)
                                    data += '{gene},{pi},{contig},{cid},{inc}\n'\
                                        .format(gene=gene,
                                                pi=percent_identity,
                                                contig=contig,
                                                cid=results['cluster_id'],
                                                inc=';'.join(sorted(inc_dict[str(results['cluster_id'])])))
                    except KeyError:
                        pass
                # If there were no results associated with the strain, make the row the strain name only
                if not result_bool:
                    data += '{sn}\n'.format(sn=sample.name)
            # Write the string to the report
            pbs.write(data)

    def __init__(self, metadata, analysistype, databasepath, threads, logfile, reportpath, matchtype='amrsummary',
                 cutoff=70):
        self.metadata = metadata
        self.analysistype = analysistype
        self.databasepath = os.path.join(databasepath, analysistype)
        self.threads = threads
        self.logfile = logfile
        self.reportpath = reportpath
        self.matchtype = matchtype
        self.cutoff = cutoff


if __name__ == '__main__':

    def resfinder_extract(reportpath, metadata):
        """
        Extract the results of the ResFinder analyses, and update the metadata object with these results
        :param reportpath: type STR: Absolute path to folder in which the report is to be found
        :param metadata: type LIST: List of metadata objects
        :return: metadata: Updated metadata object
        """
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
            # Create the necessary GenObjects
            sample.resfinder_assembled = GenObject()
            sample.resfinder_assembled.sampledata = list()
            for line in nesteddictionary:
                # Extract the strain name from the dictionary
                sample_name = nesteddictionary[line]['Strain']
                # Ensure that the loop is on the correct strain
                if sample_name == sample.name:
                    # Append the dictionary to the list of data
                    sample.resfinder_assembled.sampledata.append(([nesteddictionary[line]['Gene'],
                                                                   nesteddictionary[line]['Allele'],
                                                                   nesteddictionary[line]['Resistance'],
                                                                   nesteddictionary[line]['PercentIdentity'],
                                                                   nesteddictionary[line]['PercentCovered'],
                                                                   nesteddictionary[line]['Contig']]))
        return metadata

    def geneseekr_extract(reportpath, metadata):
        """
        Extract the results of the GeneSeekr analyses, and update the metadata object with these results
        :param reportpath: type STR: Absolute path to folder in which the report is to be found
        :param metadata: type LIST: List of metadata objects
        :return: metadata: Updated metadata object
        """
        logging.info('Extracting GeneSeekr results from reports')
        # Load the GeneSeekr outputs from the combined output file - this file contains all genes in the analysis as
        # well as the percent identity (if the gene is not present, it has a percent identity of 0)
        for sample in metadata:
            # Initialise GenObjects
            sample.geneseekr_results = GenObject()
            sample.geneseekr_results.sampledata = dict()
            report = os.path.join(reportpath, 'geneseekr_blastn.csv')
            # Open the BLAST report
            with open(report, 'r') as blast_report:
                # Create a reader object with csv.reader
                reader = csv.reader(blast_report)
                # The headers will be the first line of the report
                headers = next(reader)
                for result in reader:
                    for i, header in enumerate(headers):
                        # Ensure that the current strain matches the strain of interest, and that the iterator is not 0
                        # (keeps from adding the strain name to the dictionary)
                        if result[0] == sample.name and i != 0:
                            # Add gene name: percent identity to the dictionary
                            sample.geneseekr_results.sampledata[headers[i]] = {'PercentIdentity': result[i]}
        # Load the strain-specific BLAST outputs
        for sample in metadata:
            report = os.path.join(reportpath, '{sn}_blastn_geneseekr.tsv'.format(sn=sample.name))
            # Open the report using csv.reader, and set the headers as the first line of the report
            with open(report, 'r') as blast_report:
                reader = csv.reader(blast_report, delimiter='\t')
                headers = next(reader)
                for result in reader:
                    for i, header in enumerate(headers):
                        # Add the raw BLAST outputs (e.g. sample_id, positives, alignment_length, etc.) to the
                        # dictionary (gene name: header: result)
                        sample.geneseekr_results.sampledata[result[1]].update({headers[i]: result[i]})
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
    parser.add_argument('-a', '--analysistype',
                        choices=['amrsummary', 'geneseekr'],
                        default='amrsummary',
                        help='The analysis to perform. Options are "amrsummary" (combine ResFinder results with '
                             'MOBRecon outputs), and "geneseekr" combine geneseekr results with MOBRecon. Default is '
                             'amrsummary')
    parser.add_argument('-c', '--cutoff',
                        default=70,
                        type=int,
                        help='Integer of the cutoff value to use. Default is 70')
    SetupLogging()
    arguments = parser.parse_args()
    # Extract the list of strains in the sequence path, and create a metadata object with necessary values
    strains, metadata_object = strainer(arguments.sequencepath)
    report_path = os.path.join(arguments.sequencepath, 'reports')
    # Update the metadata object with the correct report outputs
    if arguments.analysistype == 'amrsummary':
        metadata_object = resfinder_extract(reportpath=report_path,
                                            metadata=metadata_object)
    else:
        metadata_object = geneseekr_extract(reportpath=report_path,
                                            metadata=metadata_object)
    mob = MobRecon(metadata=metadata_object,
                   analysistype='mobrecon',
                   databasepath=arguments.referencefilepath,
                   threads=multiprocessing.cpu_count() - 1,
                   logfile=os.path.join(arguments.sequencepath, 'log'),
                   reportpath=report_path,
                   matchtype=arguments.analysistype,
                   cutoff=arguments.cutoff)
    mob.mob_recon()
