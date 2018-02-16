#!/usr/bin/env python
from accessoryFunctions.accessoryFunctions import printtime, GenObject, run_subprocess, write_to_logfile, make_path, \
    MetadataObject
from argparse import ArgumentParser
from queue import Queue
import multiprocessing
from glob import glob
import json
import time
import os
__author__ = 'adamkoziol'


class Sistr(object):

    def sistr(self):
        """Perform sistr analyses on Salmonella"""
        printtime('Performing sistr analyses', self.start)
        for sample in self.metadata:
            # Create the analysis-type specific attribute
            setattr(sample, self.analysistype, GenObject())
            if sample.general.bestassemblyfile != 'NA':
                try:
                    # Only process strains that have been determined to be Salmonella
                    if sample.general.referencegenus == 'Salmonella':
                        # Set and create the path of the directory to store the strain-specific reports
                        sample[self.analysistype].reportdir = os.path.join(sample.general.outputdirectory,
                                                                           self.analysistype)
                        # Name of the .json output file
                        sample[self.analysistype].jsonoutput = os.path.join(sample[self.analysistype].reportdir,
                                                                            '{}.json'.format(sample.name))
                        # Set the sistr system call
                        sample.commands.sistr = \
                            'sistr -f json -o {} -t {} -T {} {}'\
                            .format(sample[self.analysistype].jsonoutput,
                                    self.cpus,
                                    os.path.join(sample[self.analysistype].reportdir, 'tmp'),
                                    sample.general.bestassemblyfile)
                        #
                        sample[self.analysistype].logout = os.path.join(sample[self.analysistype].reportdir, 'logout')
                        sample[self.analysistype].logerr = os.path.join(sample[self.analysistype].reportdir, 'logerr')
                        # Only run the analyses if the output json file does not exist
                        if not os.path.isfile(sample[self.analysistype].jsonoutput):
                            out, err = run_subprocess(sample.commands.sistr)
                            write_to_logfile(sample.commands.sistr, sample.commands.sistr, self.logfile,
                                             sample.general.logout, sample.general.logerr,
                                             sample[self.analysistype].logout, sample[self.analysistype].logerr)
                            write_to_logfile(out, err, self.logfile, sample.general.logout, sample.general.logerr,
                                             sample[self.analysistype].logout, sample[self.analysistype].logerr)
                        self.queue.task_done()
                except (ValueError, KeyError):
                    pass
        self.queue.join()
        self.report()

    def report(self):
        """Creates sistr reports"""
        # Initialise strings to store report data
        header = '\t'.join(self.headers) + '\n'
        data = ''
        for sample in self.metadata:
            if sample.general.bestassemblyfile != 'NA':
                # Each strain is a fresh row
                row = ''
                try:
                    # Read in the output .json file into the metadata
                    sample[self.analysistype].jsondata = json.load(open(sample[self.analysistype].jsonoutput, 'r'))
                    # Set the name of the report.
                    # Note that this is a tab-separated file, as there can be commas in the results
                    sample[self.analysistype].report = os.path.join(sample[self.analysistype].reportdir,
                                                                    '{}.tsv'.format(sample.name))
                    # Iterate through all the headers to use as keys in the json-formatted output
                    for category in self.headers:
                        # Tab separate all the results
                        row += '{}\t'.format(sample[self.analysistype].jsondata[0][category])
                        # Create attributes for each category
                        setattr(sample[self.analysistype], category,
                                str(sample[self.analysistype].jsondata[0][category]))
                    # End the results with a newline
                    row += '\n'

                    data += row
                    # Create and write headers and results to the strain-specific report
                    with open(sample[self.analysistype].report, 'w') as strainreport:
                        strainreport.write(header)
                        strainreport.write(row)
                except (KeyError, AttributeError):
                    pass
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
        # Run the analyses
        self.sistr()


if __name__ == '__main__':
    # Parser for arguments
    parser = ArgumentParser(description='Automate sistr analyses on a folder of .fasta files')
    parser.add_argument('path',
                        help='Specify input directory')
    parser.add_argument('-s', '--sequencepath',
                        required=True,
                        help='Path of .fastq(.gz) files to process.')
    # Get the arguments into an object
    arguments = parser.parse_args()

    # Define the start time
    arguments.starttime = time.time()

    # Find the files
    fastas = sorted(glob(os.path.join(arguments.sequencepath, '*.fa*')))

    # Create a metadata object
    arguments.runmetadata = MetadataObject()
    arguments.runmetadata.samples = list()
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
        arguments.runmetadata.samples.append(metadata)

    arguments.cpus = multiprocessing.cpu_count()
    arguments.reportpath = os.path.join(arguments.path, 'reports')

    # Run the script
    Sistr(arguments, 'sistr')

    # Print a bold, green exit statement
    print('\033[92m' + '\033[1m' + "\nElapsed Time: %0.2f seconds" % (time.time() - arguments.starttime) + '\033[0m')
