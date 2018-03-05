#!/usr/bin/env python
from accessoryFunctions.accessoryFunctions import MetadataObject, printtime
from spadespipeline.GeneSeekr import GeneSeekr
from csv import DictReader
from glob import glob
import operator
import os

__author__ = 'adamkoziol'


class CoreGenome(GeneSeekr):

    def blastparser(self, report, sample):
        # Open the sequence profile file as a dictionary
        blastdict = DictReader(open(report), fieldnames=self.fieldnames, dialect='excel-tab')
        resultdict = dict()
        coregenomes = list()
        # Create a list of all the names of the database files - glob, replace - with _, remove path and extension
        corefiles = glob(
            os.path.join(self.referencefilepath, self.analysistype, sample.general.referencegenus, '*.tfa'))
        for fasta in corefiles:
            fastaname = os.path.basename(os.path.splitext(fasta)[0]).replace('-', '_')
            fastaname = fastaname.split('.')[0]
            coregenomes.append(fastaname)
        # Initialise resultdict with an integer for every database file
        for genome in coregenomes:
            resultdict[genome] = int()
        # A set to store the number of core genes
        coregenes = set()
        # Go through each BLAST result
        for row in blastdict:
            # Calculate the percent identity and extract the bitscore from the row
            # Percent identity is the (length of the alignment - number of mismatches) / total subject length
            # noinspection PyTypeChecker
            percentidentity = float('{:0.2f}'.format((float(row['positives']) - float(row['gaps'])) /
                                                     float(row['subject_length']) * 100))
            # Split off any | from the sample name
            target = row['subject_id'].split('|')[0]
            # As there are variable numbers of _ in the name, try to split off the last one only if there are multiple
            # and only keep the first part of the split if there is one _ in the name
            underscored = '_'.join(target.split('_')[:-1]) if len(target.split('_')) > 2 else target.split('_')[0]
            try:
                # Since the number of core genes is the same for each reference strain, only need to determine it once
                if underscored == sorted(coregenomes)[0]:
                    coregenes.add(target)
                # If the percent identity is greater than the cutoff - adjust the cutoff to 90% for these analyses
                self.cutoff = 90
                if percentidentity >= self.cutoff:
                    # Update the dictionary with the target and the number of hits
                    resultdict[underscored] += 1
            except (KeyError, IndexError):
                pass
        # Sort the dictionary on the number of hits - best at the top
        topcore = sorted(resultdict.items(), key=operator.itemgetter(1), reverse=True)
        # Initialise a dictionary attribute to store results
        sample[self.analysistype].blastresults = dict()
        # If there are no results, populate negative results
        if not resultdict:
            sample[self.analysistype].blastresults = 'NA'
        # If results, add a string of the best number of hits, and a string of the total number of genes
        else:
            sample[self.analysistype].blastresults[topcore[0][0]] = (str(topcore[0][1]), str(len(coregenes)))

    def reporter(self):
        header = 'Strain,ClosestRef,GenesPresent/Total,\n'
        data = ''
        for sample in self.metadata:
            try:
                if sample[self.analysistype].blastresults != 'NA':
                    if sample.general.closestrefseqgenus == 'Listeria':
                        # Write the sample name, closest ref genome, and the # of genes found / total # of genes
                        closestref = list(sample[self.analysistype].blastresults.items())[0][0]
                        coregenes = list(sample[self.analysistype].blastresults.items())[0][1][0]
                        totalcore = list(sample[self.analysistype].blastresults.items())[0][1][1]
                        # Add the data to the object
                        sample[self.analysistype].targetspresent = coregenes
                        sample[self.analysistype].totaltargets = totalcore
                        sample[self.analysistype].coreresults = '{}/{}'.format(coregenes, totalcore)
                        row = '{},{},{}/{}\n'.format(sample.name, closestref, coregenes, totalcore)
                        # If the script is being run as part of the assembly pipeline, make a report for each sample
                        if self.pipeline:
                            # Open the report
                            with open(os.path.join(sample[self.analysistype].reportdir,
                                                   '{}_{}.csv'.format(sample.name, self.analysistype)), 'w') as report:
                                # Write the row to the report
                                report.write(header)
                                report.write(row)
                        data += row
                    else:
                        sample[self.analysistype].targetspresent = 'NA'
                        sample[self.analysistype].totaltargets = 'NA'
                        sample[self.analysistype].coreresults = 'NA'
            except KeyError:
                sample[self.analysistype].targetspresent = 'NA'
                sample[self.analysistype].totaltargets = 'NA'
                sample[self.analysistype].coreresults = 'NA'
        with open(os.path.join(self.reportpath, 'coregenome.csv'), 'w') as report:
            # Write the data to the report
            report.write(header)
            report.write(data)


class AnnotatedCore(object):

    def annotatedcore(self):
        """
        Calculates the core genome of organisms using custom databases
        """
        printtime('Calculating annotated core', self.start)
        # Iterate through all the samples, and process all Escherichia
        for sample in self.metadata:
            if sample.general.bestassemblyfile != 'NA':
                if sample.general.referencegenus == 'Escherichia':
                    # Add the Escherichia sample to the runmetadata
                    self.runmetadata.samples.append(sample)
                    # Create a set to store the names of all the core genes in this strain
                    sample[self.analysistype].coreset = set()
                    # Parse the BLAST report
                    self.blastparser(sample[self.analysistype].report, sample)
            # Create the report
        self.reporter()

    def blastparser(self, report, sample):
        """
        Parse the number of core genes present in the strain from the BLAST outputs
        :param report: the name and path of the BLAST outputs
        :param sample: the sample object
        """
        # Open the sequence profile file as a dictionary
        blastdict = DictReader(open(report), fieldnames=self.fieldnames, dialect='excel-tab')
        # Create a list of all the names of the database files - glob, remove path and extension
        self.coregenomes = list(map(lambda x: os.path.basename(x).split('.')[0],
                                    glob(os.path.join(self.reffilepath,
                                                      self.analysistype,
                                                      sample.general.referencegenus,
                                                      '*.tfa'))))
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
            # Remove the messy blast results from the object
            try:
                delattr(sample[self.analysistype], "blastresults")
                delattr(sample[self.analysistype], 'coreset')
            except KeyError:
                pass

    def __init__(self, inputobject):
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
        # Run the analyses
        self.annotatedcore()
