#!/usr/bin/env python
from accessoryFunctions.accessoryFunctions import GenObject, run_subprocess, write_to_logfile, make_path, \
    MetadataObject
from argparse import ArgumentParser
from click import progressbar
from queue import Queue
import multiprocessing
from glob import glob
import logging
import json
import time
import os
__author__ = 'adamkoziol'


class Sistr(object):

    def main(self):
        """
        Run the functions
        """
        self.sistr()
        self.report()

    def sistr(self):
        """
        Perform sistr analyses on Salmonella
        """
        logging.info('Performing sistr analyses')
        with progressbar(self.metadata) as bar:
            for sample in bar:
                # Create the analysis-type specific attribute
                setattr(sample, self.analysistype, GenObject())
                if sample.general.bestassemblyfile != 'NA':
                    try:
                        # Only process strains that have been determined to be Salmonella
                        if sample.general.referencegenus == 'Salmonella' or sample.general.closestrefseqgenus \
                                == 'Salmonella':
                            # Set and create the path of the directory to store the strain-specific reports
                            sample[self.analysistype].reportdir = os.path.join(sample.general.outputdirectory,
                                                                               self.analysistype)
                            # Name of the .json output file
                            sample[self.analysistype].jsonoutput = os.path.join(sample[self.analysistype].reportdir,
                                                                                '{sn}.json'.format(sn=sample.name))
                            # Set the sistr system call
                            sample.commands.sistr = \
                                'sistr -f json -o {out} -t {threads} -T {tmp_dir} {in_file}'\
                                .format(out=sample[self.analysistype].jsonoutput,
                                        threads=self.cpus,
                                        tmp_dir=os.path.join(sample[self.analysistype].reportdir, 'tmp'),
                                        in_file=sample.general.bestassemblyfile)
                            # Create the log out and err attributes
                            sample[self.analysistype].logout = os.path.join(sample[self.analysistype].reportdir,
                                                                            'logout')
                            sample[self.analysistype].logerr = os.path.join(sample[self.analysistype].reportdir,
                                                                            'logerr')
                            # Only run the analyses if the output json file does not exist
                            if not os.path.isfile(sample[self.analysistype].jsonoutput):
                                out, err = run_subprocess(sample.commands.sistr)
                                write_to_logfile(out=sample.commands.sistr,
                                                 err=sample.commands.sistr,
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
                            self.queue.task_done()
                    except (ValueError, KeyError):
                        pass
        self.queue.join()

    def report(self):
        """
        Create sistr reports
        """
        # Initialise strings to store report data
        header = '\t'.join(self.headers) + '\n'
        data = str()
        for sample in self.metadata:
            if sample.general.bestassemblyfile != 'NA':
                # Each strain is a fresh row
                row = str()
                try:
                    # Read in the output .json file into the metadata
                    sample[self.analysistype].jsondata = json.load(open(sample[self.analysistype].jsonoutput, 'r'))
                    # Set the name of the report.
                    # Note that this is a tab-separated file, as there can be commas in the results
                    sample[self.analysistype].report = os.path.join(sample[self.analysistype].reportdir,
                                                                    '{sn}.tsv'.format(sn=sample.name))
                    # Iterate through all the headers to use as keys in the json-formatted output
                    for category in self.headers:
                        # Tab separate all the results
                        row += '{result}\t'.format(result=sample[self.analysistype].jsondata[0][category])
                        # Create attributes for each category
                        setattr(sample[self.analysistype], category,
                                str(sample[self.analysistype].jsondata[0][category]))
                    # Remove any trailing tabs, and end the results with a newline character
                    row = row.rstrip('\t') + '\n'
                    data += row
                    # Create and write headers and results to the strain-specific report
                    with open(sample[self.analysistype].report, 'w') as strainreport:
                        strainreport.write(header)
                        strainreport.write(row)
                except (KeyError, AttributeError) as e:
                    data += '{sn}\n'.format(sn=sample.name)
            else:
                data += '{sn}\n'.format(sn=sample.name)
        # Create and write headers and cumulative results to the combined report
        with open(os.path.join(self.reportdir, 'sistr.tsv'), 'w') as report:
            report.write(header)
            report.write(data)

    def __init__(self, inputobject, analysistype):
        try:
            self.metadata = inputobject.runmetadata.samples
        except AttributeError:
            self.metadata = inputobject.metadata
        try:
            self.start = inputobject.starttime
        except AttributeError:
            self.start = inputobject.start
        self.cpus = inputobject.cpus
        self.logfile = inputobject.logfile
        self.reportdir = os.path.join(inputobject.reportpath, '')
        make_path(self.reportdir)
        self.analysistype = analysistype
        # self.devnull = open(os.devnull, 'wb')
        self.queue = Queue()
        self.headers = ['genome', 'cgmlst_distance', 'cgmlst_genome_match', 'cgmlst_matching_alleles', 'h1', 'h2',
                        'serogroup', 'serovar', 'serovar_antigen', 'serovar_cgmlst']


if __name__ == '__main__':
    # Parser for arguments
    parser = ArgumentParser(description='Automate sistr analyses on a folder of .fasta files')
    parser.add_argument('path',
                        help='Output directory. Created if it does not already exist.')
    parser.add_argument('-s', '--sequencepath',
                        required=True,
                        help='Path of .fasta files to process.')
    # Get the arguments into an object
    arguments = parser.parse_args()

    # Define the start time
    arguments.starttime = time.time()

    # Find the files
    fastas = sorted(glob(os.path.join(arguments.sequencepath, '*.fa*')))

    # Create a metadata object
    arguments.runmetadata = MetadataObject()
    arguments.runmetadata.samples = list()

    # Have to add logfile to make the SISTR automator happy.
    arguments.logfile = os.path.join(arguments.path, 'sistr')

    for fasta in fastas:
        metadata = MetadataObject()
        metadata.name = os.path.split(fasta)[1].split('.')[0]
        # Initialise the general and run categories
        metadata.general = GenObject()
        metadata.run = GenObject()
        # Set the destination folder
        outputdir = os.path.join(arguments.sequencepath, metadata.name)
        make_path(outputdir)
        # Add the output directory to the metadata
        metadata.general.outputdirectory = outputdir
        metadata.run.outputdirectory = outputdir
        metadata.general.bestassemblyfile = True
        # Initialise an attribute to store commands
        metadata.commands = GenObject()
        # Assume that all samples are Salmonella
        metadata.general.referencegenus = 'Salmonella'
        # Set the .fasta file as the best assembly
        metadata.general.bestassemblyfile = fasta
        # Also set the logging attributes, otherwise everything crashes.
        metadata.general.logerr = os.path.join(arguments.path, 'log_err')
        metadata.general.logout = os.path.join(arguments.path, 'log_out')
        if not os.path.isdir(os.path.join(arguments.path, metadata.name, 'sistr')):
            make_path(os.path.join(arguments.path, metadata.name, 'sistr'))

        # Append to our list of samples
        arguments.runmetadata.samples.append(metadata)

    arguments.cpus = multiprocessing.cpu_count()
    arguments.reportpath = os.path.join(arguments.path, 'reports')

    # Run the script
    sistr = Sistr(inputobject=arguments,
                  analysistype='sistr')

    # Print a bold, blue exit statement
    logging.info('Analyses complete!')
