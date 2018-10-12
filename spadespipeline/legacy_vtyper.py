#!/usr/bin/env python3
from accessoryFunctions.accessoryFunctions import GenObject, make_path, MetadataObject
from Bio import SeqIO, Seq
from argparse import ArgumentParser
from click import progressbar
from itertools import product
from threading import Thread
from subprocess import call
from queue import Queue
from glob import glob
from time import time
import logging
import os

__author__ = 'adamkoziol'


class Vtyper(object):

    def vtyper(self):
        if not os.path.isfile(self.formattedprimers):
            self.epcr_primers()
            self.epcr_primer_file()
        self.epcr_threads()
        self.epcr_parse()
        self.epcr_report()

    def epcr_primers(self):
        """
        Read in the primer file, and create a properly formatted output file that takes any degenerate bases
        into account
        """
        logging.info('Populating primer dictionaries')
        for record in SeqIO.parse(self.primerfile, 'fasta'):
            # from https://stackoverflow.com/a/27552377 - find any degenerate bases in the primer sequence, and
            # create all possibilities as a list
            degenerates = Seq.IUPAC.IUPACData.ambiguous_dna_values
            primerlist = list(map("".join, product(*map(degenerates.get, str(record.seq).upper()))))
            # As the record.id is being updated in the loop below, set the name of the primer here so that will
            # be able to be recalled when setting the new record.ids
            primername = record.id
            # Iterate through all the possible primers created from any degenerate bases
            for primer in primerlist:
                # Split the base name of the target from the direction
                # e.g. vtx1a-F1 is split in vtx1a and F1
                basename, direction = primername.split('-')
                # Populate the dictionaries of forward and reverse primers based on the direction determined above
                if direction.startswith('F'):
                    # Attempt to add the current primer sequence to the dictionary
                    try:
                        self.forward_dict[basename].append(primer)
                    # On a key error, initialise the list of primers
                    except KeyError:
                        self.forward_dict[basename] = list()
                        self.forward_dict[basename].append(primer)
                else:
                    try:
                        self.reverse_dict[basename].append(primer)
                    except KeyError:
                        self.reverse_dict[basename] = list()
                        self.reverse_dict[basename].append(primer)

    def epcr_primer_file(self):
        """
        Create the ePCR-compatible primer file from the dictionaries of primer combinations
        """
        logging.info('Creating re-PCR-compatible primer file')
        with open(self.formattedprimers, 'w') as formatted:
            # Iterate through all the targets
            for basename in sorted(self.forward_dict):
                # Use enumerate to number the iterations for each forward and reverse primer in the lists
                for forward_index, forward_primer in enumerate(self.forward_dict[basename]):
                    for reverse_index, reverse_primer in enumerate(self.reverse_dict[basename]):
                        # Set the name of the primer using the target name, and the indices of the primers
                        # e.g. vtx1a_0_0
                        primer_name = '{bn}_{fi}_{ri}'.format(bn=basename,
                                                              fi=forward_index,
                                                              ri=reverse_index)
                        # Create the string to write to the ePCR-compatible primer file
                        # e.g. vtx1a_0_0	CCTTTCCAGGTACAACAGCGGTT	GGAAACTCATCAGATGCCATTCTGG
                        output_string = '{pn}\t{fp}\t{rp}\n'.format(pn=primer_name,
                                                                    fp=forward_primer,
                                                                    rp=reverse_primer)
                        # Write the string to file
                        formatted.write(output_string)

    def epcr_threads(self):
        """
        Run ePCR in a multi-threaded fashion
        """
        # Create the threads for the ePCR analysis
        for sample in self.metadata:
            if sample.general.bestassemblyfile != 'NA':
                threads = Thread(target=self.epcr, args=())
                threads.setDaemon(True)
                threads.start()
        logging.info('Running ePCR analyses')
        for sample in self.metadata:
            if sample.general.bestassemblyfile != 'NA':
                setattr(sample, self.analysistype, GenObject())
                # Get the primers ready
                sample[self.analysistype].primers = self.formattedprimers
                # Make the output path
                sample[self.analysistype].reportdir = os.path.join(sample.general.outputdirectory,
                                                                   self.analysistype)
                make_path(sample[self.analysistype].reportdir)
                outfile = os.path.join(sample[self.analysistype].reportdir, sample.name)
                # Set the hashing and mapping commands
                sample.commands.famap = '{famap} -b {outfile}.famap {fasta}'\
                    .format(famap=os.path.join(self.homepath, 'ePCR', 'famap'),
                            outfile=outfile,
                            fasta=sample.general.bestassemblyfile)
                sample.commands.fahash = '{fahash} -b {outfile}.hash {outfile}.famap'\
                    .format(fahash=os.path.join(self.homepath, 'ePCR', 'fahash'),
                            outfile=outfile)
                # re-PCR uses the subtyping primers list to search the contigs file using the following parameters
                # -S {hash file} (Perform STS lookup using hash-file), -r + (Enable/disable reverse STS lookup)
                # -m 10000 (Set variability for STS size for lookup), this very large, as I don't necessarily know
                # the size of the amplicon
                # -n 1 (Set max allowed mismatches per primer pair for lookup)
                # -g 0 (Set max allowed indels per primer pair for lookup),
                # -G (Print alignments in comments)
                # -o {output file}
                sample.commands.epcr = \
                    '{rePCR} -S {outfile}.hash -r + -m 10000 -n 1 -g 0 -G -q -o {outfile}.txt {primers}'\
                    .format(rePCR=os.path.join(self.homepath, 'ePCR', 're-PCR'),
                            outfile=outfile,
                            primers=sample[self.analysistype].primers)
                sample[self.analysistype].resultsfile = '{}.txt'.format(outfile)
                # Add the sample object and the output file to the queue
                self.epcrqueue.put((sample, outfile))
        # Join the threads
        self.epcrqueue.join()

    def epcr(self):
        while True:
            sample, linkfile = self.epcrqueue.get()
            # Run the commands if the ePCR output file doesn't exist
            if not os.path.isfile('{}.txt'.format(linkfile)):
                call(sample.commands.famap, shell=True, stdout=self.devnull, stderr=self.devnull)
                call(sample.commands.fahash, shell=True, stdout=self.devnull, stderr=self.devnull)
                call(sample.commands.epcr, shell=True, stdout=self.devnull, stderr=self.devnull)
            # Clean up the temporary files
            try:
                os.remove('{}.famap'.format(linkfile))
            except FileNotFoundError:
                pass
            try:
                os.remove('{}.hash'.format(linkfile))
            except FileNotFoundError:
                pass
            # Signal that the thread is complete
            self.epcrqueue.task_done()

    def epcr_parse(self):
        """
        Parse the ePCR outputs
        """
        logging.info('Parsing ePCR outputs')
        for sample in self.metadata:
            if sample.general.bestassemblyfile != 'NA':
                # Create a set to store all the unique results
                toxin_set = set()
                if os.path.isfile(sample[self.analysistype].resultsfile):
                    with open(sample[self.analysistype].resultsfile) as epcrresults:
                        for result in epcrresults:
                            # Only the lines without a # contain results
                            if "#" not in result:
                                # Split on \t
                                data = result.split('\t')
                                # The subtyping primer pair is the first entry on lines with results
                                vttype = data[0].split('_')[0]
                                # Add the verotoxin subtype to the set of detected subtypes
                                toxin_set.add(vttype)
                # Create a string of the entries in the sorted list of toxins joined with ";"
                sample[self.analysistype].toxinprofile = ";".join(sorted(list(toxin_set))) if toxin_set else 'ND'
            else:
                setattr(sample, self.analysistype, GenObject())
                sample[self.analysistype].toxinprofile = 'NA'

    def epcr_report(self):
        """
        Create a report of the ePCR-calculated toxin profiles
        """
        logging.info('Creating {at} report'.format(at=self.analysistype))
        with open(os.path.join(self.reportpath, '{at}.csv'.format(at=self.analysistype)), 'w') as report:
            data = 'Strain,ToxinProfile\n'
            for sample in self.metadata:
                data += '{sn},{tp}\n'.format(sn=sample.name,
                                             tp=sample[self.analysistype].toxinprofile)
            # Write the data to the report
            report.write(data)

    def __init__(self, inputobject, analysistype):
        self.metadata = inputobject.runmetadata.samples
        self.analysistype = analysistype
        self.start = inputobject.starttime
        self.reportpath = inputobject.reportpath
        make_path(self.reportpath)
        self.devnull = open(os.devnull, 'wb')
        self.epcrqueue = Queue()
        # Extract the path of the current script from the full path + file name
        self.homepath = os.path.split(os.path.abspath(__file__))[0]
        self.primerfile = os.path.join(self.homepath, 'primers.txt')
        self.formattedprimers = os.path.join(self.homepath, 'ssi_subtyping_primers.txt')
        self.forward_dict = dict()
        self.reverse_dict = dict()


if __name__ == '__main__':

    def filer(args):
        """
        Create metadata objects with necessary attributes for each FASTA file found in the sequence path
        :param args: Argument parser object with necessary variables
        :return: samples: List of metadata objects
        """
        # List to store all the metadata objects
        samples = list()
        # Find all the sequence files in the path
        fastas = sorted(glob(os.path.join(args.sequencepath, '*.fa*')))
        for fasta in fastas:
            # Create a metadata object for each sample
            metadata = MetadataObject()
            # Populate the metadata object with the required attributes
            metadata.name = os.path.splitext(os.path.basename(fasta))[0]
            metadata.general = GenObject()
            metadata.commands = GenObject()
            metadata.general.bestassemblyfile = fasta
            metadata.general.outputdirectory = os.path.join(args.sequencepath, metadata.name)
            samples.append(metadata)
        return samples

    def argument_parser():

        # Parser for arguments
        parser = ArgumentParser(description='Perform verotoxin sub-typing on FASTA files')
        parser.add_argument('-s', '--sequencepath',
                            required=True,
                            help='Path to folder containing sequencing reads')
        # Get the arguments into an object
        arguments = parser.parse_args()
        arguments.starttime = time()
        arguments.reportpath = os.path.join(arguments.sequencepath, 'reports')
        arguments.runmetadata = MetadataObject()
        # Create metadata objects for the samples
        arguments.runmetadata.samples = filer(arguments)
        # Perform vtx typing
        vtyper = Vtyper(arguments, 'vtyper_legacy')
        vtyper.vtyper()

    # Run the script
    argument_parser()
