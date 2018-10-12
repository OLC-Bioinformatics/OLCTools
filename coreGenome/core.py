#!/usr/bin/env python
from accessoryFunctions.accessoryFunctions import MetadataObject
from geneseekr.blast import BLAST
from Bio import SeqIO
from csv import DictReader
from glob import glob
import operator
import logging
import os

__author__ = 'adamkoziol'


class CoreGenome(BLAST):

    def parse_results(self):
        self.parser(self.metadata,
                    self.analysistype,
                    self.fieldnames,
                    self.cutoff,
                    self.program)

    @staticmethod
    def parser(metadata, analysistype, fieldnames, cutoff, program):
        for sample in metadata:
            # Initialise a dictionary attribute to store results
            sample[analysistype].blastresults = dict()
            try:
                # Open the sequence profile file as a dictionary
                blastdict = DictReader(open(sample[analysistype].report), fieldnames=fieldnames, dialect='excel-tab')
                resultdict = dict()
                resultset = dict()
                # Initialise a dictionary to store all the target sequences
                sample[analysistype].targetsequence = dict()
                coregenomes = list()
                # A set to store the number of core genes
                coregenes = set()
                # Create a list of all the names of the database files, replace - with _, remove path and extension
                for fasta in sample[analysistype].targets:
                    fastaname = os.path.basename(os.path.splitext(fasta)[0]).replace('-', '_')
                    fastaname = fastaname.split('.')[0]
                    coregenomes.append(fastaname)
                # Go through each BLAST result
                for row in blastdict:
                    # Create the subject length variable - if the sequences are DNA (e.g. blastn), use the subject
                    # length as usual; if the sequences are protein (e.g. tblastx), use the subject length / 3
                    if program == 'blastn' or program == 'blastp' or program == 'blastx':
                        subject_length = float(row['subject_length'])

                    else:
                        subject_length = float(row['subject_length']) / 3
                    # Calculate the percent identity and extract the bitscore from the row
                    # Percent identity is the (length of the alignment - number of mismatches) / total subject length
                    percentidentity = float('{:0.2f}'.format((float(row['positives']) - float(row['gaps'])) /
                                                             subject_length * 100))
                    # If the percent identity is greater than the cutoff
                    if percentidentity >= cutoff:
                        # Split off any | from the sample name
                        target = row['subject_id'].split('|')[0]
                        # As there are variable _ in the name, try to split off the last one only if there are
                        #  multiple and only keep the first part of the split if there is one _ in the name
                        underscored = '_'.join(target.split('_')[:-1]) if len(target.split('_')) > 2 else \
                            target.split('_')[0]
                        try:
                            # Update the dictionary with the reference genome and the target
                            resultset[underscored].add(target)
                        except KeyError:
                            # Initialise the dictionary with the first hit
                            resultset[underscored] = set()
                            resultset[underscored].add(target)
                # Get the number of unique genes per reference genome
                for underscored, target_set in resultset.items():
                    resultdict[underscored] = len(target_set)
                # Sort the dictionary on the number of hits - best at the top
                topcore = sorted(resultdict.items(), key=operator.itemgetter(1), reverse=True)
                # If there are no results, populate negative results
                if not resultdict:
                    sample[analysistype].blastresults = 'NA'
                # If results, add a string of the best number of hits, and a string of the total number of genes
                # This is currently 1013. If this changes, I may re-implement a dynamic method of determining
                # this value
                else:
                    sample[analysistype].blastresults[topcore[0][0]] = (str(topcore[0][1]), str(1013))
            except FileNotFoundError:
                sample[analysistype].blastresults = 'NA'
        return metadata

    def create_reports(self):
        # Create dictionaries
        self.metadata = self.geneseekr.dict_initialise(self.metadata,
                                                       self.analysistype)
        self.reporter(self.metadata, self.analysistype, self.reportpath)

    @staticmethod
    def reporter(metadata, analysistype, reportpath):
        header = 'Strain,ClosestRef,GenesPresent/Total,\n'
        data = str()
        for sample in metadata:
            try:
                if sample[analysistype].blastresults != 'NA':
                    if sample.general.closestrefseqgenus == 'Listeria':
                        # Write the sample name, closest ref genome, and the # of genes found / total # of genes
                        closestref = list(sample[analysistype].blastresults.items())[0][0]
                        coregenes = list(sample[analysistype].blastresults.items())[0][1][0]
                        # Find the closest reference file
                        try:
                            ref = glob(os.path.join(sample[analysistype].targetpath, '{fasta}*'
                                                    .format(fasta=closestref)))[0]
                        except IndexError:
                            # Replace underscores with dashes to find files
                            closestref = closestref.replace('_', '-')
                            ref = glob(os.path.join(sample[analysistype].targetpath, '{fasta}*'
                                                    .format(fasta=closestref)))[0]
                        # Determine the number of core genes present in the closest reference file
                        totalcore = 0
                        for _ in SeqIO.parse(ref, 'fasta'):
                            totalcore += 1
                        # Add the data to the object
                        sample[analysistype].targetspresent = coregenes
                        sample[analysistype].totaltargets = totalcore
                        sample[analysistype].coreresults = '{}/{}'.format(coregenes, totalcore)
                        row = '{},{},{}/{}\n'.format(sample.name, closestref, coregenes, totalcore)
                        # Open the report
                        with open(os.path.join(sample[analysistype].reportdir,
                                               '{}_{}.csv'.format(sample.name, analysistype)), 'w') as report:
                            # Write the row to the report
                            report.write(header)
                            report.write(row)
                        data += row
                    else:
                        sample[analysistype].targetspresent = 'NA'
                        sample[analysistype].totaltargets = 'NA'
                        sample[analysistype].coreresults = 'NA'
            except KeyError:
                sample[analysistype].targetspresent = 'NA'
                sample[analysistype].totaltargets = 'NA'
                sample[analysistype].coreresults = 'NA'
        with open(os.path.join(reportpath, 'coregenome.csv'), 'w') as report:
            # Write the data to the report
            report.write(header)
            report.write(data)


class AnnotatedCore(object):

    def annotatedcore(self):
        """
        Calculates the core genome of organisms using custom databases
        """
        logging.info('Calculating annotated core')
        # Determine the total number of core genes
        self.total_core()
        # Iterate through all the samples, and process all Escherichia
        for sample in self.metadata:
            if sample.general.bestassemblyfile != 'NA':
                # Create a set to store the names of all the core genes in this strain
                sample[self.analysistype].coreset = set()
                if sample.general.referencegenus == 'Escherichia':
                    # Add the Escherichia sample to the runmetadata
                    self.runmetadata.samples.append(sample)
                    # Parse the BLAST report
                    try:
                        report = sample[self.analysistype].report
                        self.blastparser(report, sample)
                    except KeyError:
                        sample[self.analysistype].coreset = list()
        # Create the report
        self.reporter()

    def total_core(self):
        """
        Determine the total number of core genes present
        """
        corefile = os.path.join(self.reffilepath, self.analysistype, 'Escherichia', 'core_combined.fasta')
        for record in SeqIO.parse(corefile, 'fasta'):
            gene_name = record.id.split('-')[0]
            if gene_name not in self.coregenomes:
                self.coregenomes.append(gene_name)

    def blastparser(self, report, sample):
        """
        Parse the number of core genes present in the strain from the BLAST outputs
        :param report: the name and path of the BLAST outputs
        :param sample: the sample object
        """
        try:
            # Open the sequence profile file as a dictionary
            blastdict = DictReader(open(report), fieldnames=self.fieldnames, dialect='excel-tab')
            # Go through each BLAST result
            for row in blastdict:
                # Calculate the percent identity and extract the bitscore from the row
                # Percent identity is the (length of the alignment - number of mismatches) / total subject length
                percentidentity = float('{:0.2f}'.format((float(row['positives']) - float(row['gaps'])) /
                                                         float(row['subject_length']) * 100))
                # Split off any | and - from the sample name
                target = row['subject_id'].split('|')[0].split('-')[0]
                # If the hit passes the cutoff threshold, add it to the set of core genes present
                if percentidentity >= self.cutoff:
                    sample[self.analysistype].coreset.add(target)
        except FileNotFoundError:
            pass

    def reporter(self):
        """
        Create a .csv file with the strain name, and the number of core genes present/the total number of core genes
        """
        with open(os.path.join(self.reportpath, 'Escherichia_core.csv'), 'w') as report:
            data = 'Strain,Genes Present/Total\n'
            for sample in self.runmetadata.samples:
                # Convert the set to a list for JSON serialization
                sample[self.analysistype].coreset = list(sample[self.analysistype].coreset)
                sample[self.analysistype].coreresults = '{}/{}'.format(len(sample[self.analysistype].coreset),
                                                                       len(self.coregenomes))
                # Add strain name, the number of core genes present, and the number of total core genes to the string
                data += '{},{}\n'.format(sample.name, sample[self.analysistype].coreresults)
            report.write(data)

        for sample in self.metadata:
            # Remove the messy blast results and set/list of core genes from the object
            try:
                delattr(sample[self.analysistype], "blastresults")
            except AttributeError:
                pass
            try:
                delattr(sample[self.analysistype], 'coreset')
            except AttributeError:
                pass

    def __init__(self, inputobject, genus_specific=False):
        self.start = inputobject.starttime
        self.commit = inputobject.commit
        self.starttime = inputobject.starttime
        self.homepath = inputobject.homepath
        self.path = inputobject.path
        self.cpus = inputobject.cpus
        self.metadata = inputobject.runmetadata.samples
        self.runmetadata = MetadataObject()
        self.runmetadata.samples = list()
        self.reffilepath = inputobject.reffilepath
        self.reportpath = inputobject.reportpath
        self.logfile = inputobject.logfile
        self.analysistype = 'coregenome'
        self.cutoff = 90
        self.coregenomes = list()
        # Fields used for custom outfmt 6 BLAST output:
        self.fieldnames = ['query_id', 'subject_id', 'positives', 'mismatches', 'gaps',
                           'evalue', 'bit_score', 'subject_length', 'alignment_length',
                           'query_start', 'query_end', 'query_sequence',
                           'subject_start', 'subject_end', 'subject_sequence']
        self.genus_specific = genus_specific
        # Run the analyses
        self.annotatedcore()
