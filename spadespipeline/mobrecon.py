#!/usr/bin/env python3
from accessoryFunctions.accessoryFunctions import GenObject, run_subprocess, write_to_logfile, make_path
from click import progressbar
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
            data = 'Strain,Gene,Allele,Resistance,PercentIdentity,Contig,Location,Incompatibility\n'
            for sample in self.metadata:
                for primarykey, results in sample[self.analysistype].report_dict.items():
                    try:
                        contig = results['contig_id'].split('|')[1]
                        for amr_l in sample.resfinder_assembled.sampledata:
                            if contig in amr_l:
                                data += '{},'.format(sample.name)
                                data += '{amr},{mob}\n'\
                                    .format(amr=','.join(str(res) if str(res) != 'nan' else 'ND' for res in amr_l[:4]),
                                            mob=','.join(str(res) if str(res) != 'nan' else 'ND' for res in
                                                         [contig, results['cluster_id'], results['rep_type']]))
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
