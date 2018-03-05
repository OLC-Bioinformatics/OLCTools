#!/usr/bin/env python3
from accessoryFunctions.accessoryFunctions import make_path, run_subprocess, write_to_logfile, logstr, GenObject, \
    printtime, get_version, MetadataObject
from spadespipeline.bowtie import Bowtie2BuildCommandLine, Bowtie2CommandLine
from Bio.Sequencing.Applications import SamtoolsViewCommandline, SamtoolsSortCommandline, SamtoolsIndexCommandline
from Bio.Application import ApplicationError
from Bio import SeqIO
from threading import Lock, Thread
from io import StringIO
from glob import glob
import threading
import shutil
import os
import re

__author__ = 'mike knowles, adamkoziol'


class QualiMap(object):

    def main(self):
        """
        Run the methods in the correct order
        """
        printtime('Aligning reads with bowtie2 for Qualimap', self.start)
        self.bowtie()
        self.indexing()
        self.pilon()
        self.filter()
        self.clear()

    def bowtie(self):
        """
        Create threads and commands for performing reference mapping for qualimap analyses
        """
        for i in range(self.cpus):
            # Send the threads to the merge method. :args is empty as I'm using
            threads = Thread(target=self.align, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        for sample in self.metadata:
            # Initialise the mapping GenObject
            sample.mapping = GenObject()
            # Set an easier to write shortcut for sample.general
            sagen = sample.general
            if sagen.bestassemblyfile != "NA":
                sagen.QualimapResults = os.path.join(sagen.outputdirectory, 'qualimap_results')
                # Set the results folder
                # Create this results folder if necessary
                make_path(sagen.QualimapResults)
                # Set file names
                sagen.sortedbam = os.path.join(sagen.QualimapResults, '{}_sorted.bam'.format(sample.name))
                filenoext = os.path.splitext(sagen.filteredfile)[0]
                sagen.filenoext = filenoext
                sagen.bowtie2results = os.path.join(sagen.QualimapResults, sample.name)
                # Use fancy new bowtie2 wrapper
                bowtie2build = Bowtie2BuildCommandLine(reference=sagen.bestassemblyfile,
                                                       bt2=sagen.bowtie2results)
                sample.mapping.BamFile = sagen.bowtie2results + "_sorted.bam"
                # SAMtools sort v1.3 has different run parameters
                if self.samversion < "1.3":
                    samsort = SamtoolsSortCommandline(input="-", out_prefix=sample.mapping.BamFile)
                else:
                    samsort = SamtoolsSortCommandline(input=sample.mapping.BamFile,
                                                      o=True,
                                                      out_prefix="-")
                samtools = [SamtoolsViewCommandline(b=True, S=True, input_file="-"), samsort]
                indict = {'D': 5, 'R': 1, 'num_mismatches': 0, 'seed_length': 22, 'i_func': "S,0,2.50"}
                #  Update the dictionary with the appropriate parameters for paired- vs. single-ended assemblies
                try:
                    _ = sample.general.mergedreads
                    if len(sample.general.trimmedcorrectedfastqfiles) == 2:
                        indict.update({'m1': sample.general.trimmedcorrectedfastqfiles[0],
                                       'm2': sample.general.trimmedcorrectedfastqfiles[1]})
                    else:
                        indict.update({'U': sample.general.trimmedcorrectedfastqfiles[0]})
                except KeyError:
                    if len(sample.general.assemblyfastq) == 2:
                        indict.update({'m1': sample.general.assemblyfastq[0], 'm2': sample.general.assemblyfastq[1]})
                    else:
                        indict.update({'U': sample.general.assemblyfastq[0]})
                bowtie2align = Bowtie2CommandLine(bt2=sagen.bowtie2results,
                                                  threads=self.threads,
                                                  samtools=samtools,
                                                  **indict)

                # Convert the commands to strings to allow them to be JSON serialized
                sample.commands.bowtie2align = str(bowtie2align)
                sample.commands.bowtie2build = str(bowtie2build)
                self.bowqueue.put((sample, sample.commands.bowtie2build, sample.commands.bowtie2align))
            else:
                sample.commands.samtools = "NA"
                sample.mapping.MeanInsertSize = 'NA'
                sample.mapping.MeanCoveragedata = 'NA'
        self.bowqueue.join()

    def align(self):
        # from subprocess import call
        while True:
            sample, bowtie2build, bowtie2align = self.bowqueue.get()
            if sample.general.bestassemblyfile != 'NA':
                if not os.path.isfile(sample.mapping.BamFile) and not os.path.isfile(sample.mapping.BamFile + ".bz2"):
                    stdout = StringIO()
                    for func in bowtie2build, bowtie2align:
                        stdout.close()
                        out, err = run_subprocess(str(func))
                        self.threadlock.acquire()
                        write_to_logfile(str(func), str(func), self.logfile, sample.general.logout,
                                         sample.general.logerr, None, None)
                        write_to_logfile(out, err, self.logfile, sample.general.logout, sample.general.logerr,
                                         None, None)
                        self.threadlock.release()
            # For different alignment
            sam = sample.general.bowtie2results + ".sam"
            if os.path.isfile(sam):
                # PIPE stdout to stdin of samtools view then sort (only outputting sorted bam)
                # SAMtools sort v1.3 has different run parameters
                if self.samversion < "1.3":
                    samsort = SamtoolsSortCommandline(input="-", out_prefix=sample.mapping.BamFile[:-4])
                else:
                    samsort = SamtoolsSortCommandline(input=sample.mapping.BamFile, o=True, out_prefix="-")
                # Use cStringIO streams to handle bowtie output
                stdout = StringIO()
                for func in [SamtoolsViewCommandline(b=True, S=True, input_file=sample.mapping.BamFile), samsort]:
                    # Use closing contextmanager for handle __exit__() as close()
                    stdout, stderr = map(StringIO, func(stdin=stdout.getvalue()))
                    # Write the standard error to log
                    with open(os.path.join(sample.general.QualimapResults, "samtools.log"), "ab+") as log:
                        log.writelines(logstr(func, stderr.getvalue()))
                    stderr.close()
                stdout.close()
            self.mapper(sample)
            # Signal to the queue that the job is done
            self.bowqueue.task_done()

    def mapper(self, sample):
        """
        Run qualimap and parse the outputs
        :param sample: metadata object
        """
        if sample.general.bestassemblyfile != "NA":
            # Define the Qualimap log and report files
            reportfile = os.path.join(sample.general.QualimapResults, 'genome_results.txt')
            # Define the Qualimap call
            qualimapcall = 'qualimap bamqc -bam {} -outdir {}'.format(sample.general.sortedbam,
                                                                      sample.general.QualimapResults)
            sample.commands.qualimap = qualimapcall
            # Initialise a dictionary to hold the Qualimap results
            qdict = dict()
            # If the report file doesn't exist, run Qualimap, and print logs to the log file
            if not os.path.isfile(reportfile):
                tlock = threading.Lock()
                out, err = run_subprocess(sample.commands.qualimap)
                tlock.acquire()
                write_to_logfile(sample.commands.qualimap, sample.commands.qualimap, self.logfile,
                                 sample.general.logout, sample.general.logerr, None, None)
                write_to_logfile(out, err, self.logfile, sample.general.logout, sample.general.logerr, None, None)
                tlock.release()
            # Initialise a genobject to store the coverage dictionaries
            sample.depth = GenObject()
            sample.depth.length = dict()
            sample.depth.bases = dict()
            sample.depth.coverage = dict()
            sample.depth.stddev = dict()
            try:
                with open(reportfile) as report:
                    # Read the report
                    for line in report:
                        # Sanitise the keys and values using self.analyze
                        key, value = self.analyze(line)
                        # If the keys and values exist, enter them into the dictionary
                        if (key, value) != (None, None):
                            qdict[key] = value
                        if 'Coverage per contig' in line:
                            for contigline in report:
                                try:
                                    _, name, length, bases, coverage, stddev = contigline.rstrip().split('\t')
                                    sample.depth.length.update({name: length})
                                    sample.depth.bases.update({name: bases})
                                    sample.depth.coverage.update({name: coverage})
                                    sample.depth.stddev.update({name: stddev})
                                except ValueError:
                                    pass

            except (IOError, FileNotFoundError):
                pass
            # If there are values in the dictionary
            if qdict:
                # Make new category for Qualimap results and populate this category with the report data
                for attribute in qdict:
                    # Remove the 'X' from the depth values e.g. 40.238X
                    setattr(sample.mapping, attribute, qdict[attribute].rstrip('X'))

    def indexing(self):
        printtime('Indexing sorted bam files', self.start)
        for i in range(self.cpus):
            # Send the threads to
            threads = Thread(target=self.index, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        for sample in self.metadata:
            if sample.general.bestassemblyfile != 'NA':
                bamindex = SamtoolsIndexCommandline(input=sample.mapping.BamFile)
                sample.mapping.sortedbai = sample.mapping.BamFile + '.bai'
                sample.mapping.bamindex = str(bamindex)
                self.indexqueue.put((sample, bamindex))
        self.indexqueue.join()

    def index(self):
        while True:
            try:
                sample, bamindex = self.indexqueue.get()
                # Only make the call if the .bai file doesn't already exist
                if not os.path.isfile(sample.mapping.sortedbai):
                    # Use cStringIO streams to handle bowtie output
                    stdout, stderr = map(StringIO, bamindex(cwd=sample.general.QualimapResults))
                    if stderr:
                        # Write the standard error to log
                        with open(os.path.join(sample.general.QualimapResults,
                                               'indexing_samtools_bam_index.log'), 'a+') as log:
                            log.writelines(logstr(bamindex, stderr.getvalue(), stdout.getvalue()))
                    stderr.close()
            except ApplicationError:
                pass
            self.indexqueue.task_done()

    def pilon(self):
        """
        Run pilon to fix any misassemblies in the contigs - will look for SNPs and indels
        """
        printtime('Improving quality of assembly with pilon', self.start)
        for i in range(self.cpus):
            # Send the threads to the merge method. :args is empty as I'm using
            threads = Thread(target=self.pilonthreads, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        for sample in self.metadata:
            if sample.general.bestassemblyfile != 'NA':
                # Set the name of the unfiltered spades assembly output file
                sample.general.contigsfile = os.path.join(sample.general.spadesoutput, 'contigs.fasta')
                sample.mapping.pilondir = os.path.join(sample.general.QualimapResults, 'pilon')
                make_path(sample.mapping.pilondir)
                # Create the command line command
                sample.mapping.piloncmd = 'pilon --genome {} --bam {} --fix bases --threads {} ' \
                                          '--outdir {} --changes --mindepth 0.25' \
                    .format(sample.general.contigsfile,
                            sample.mapping.BamFile,
                            self.threads,
                            sample.mapping.pilondir)
                self.pilonqueue.put(sample)
        self.pilonqueue.join()

    def pilonthreads(self):
        while True:
            sample = self.pilonqueue.get()
            sample.general.contigsfile = os.path.join(sample.mapping.pilondir, 'pilon.fasta')
            # Only perform analyses if the output file doesn't already exist
            if not os.path.isfile(sample.general.contigsfile):
                command = sample.mapping.piloncmd
                out, err = run_subprocess(command)
                self.threadlock.acquire()
                write_to_logfile(command, command, self.logfile, sample.general.logout, sample.general.logerr, None,
                                 None)
                write_to_logfile(out, err, self.logfile, sample.general.logout, sample.general.logerr, None, None)
                self.threadlock.release()
            self.pilonqueue.task_done()

    def filter(self):
        """
        Filter contigs based on depth
        """

        printtime('Filtering contigs', self.start)
        for i in range(self.cpus):
            # Send the threads to the filter method
            threads = Thread(target=self.filterthreads, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        for sample in self.metadata:
            # Set the name of the unfiltered spades assembly output file
            if sample.general.bestassemblyfile != 'NA':
                sample.general.contigsfile = os.path.join(sample.general.spadesoutput, 'contigs.fasta')
                self.filterqueue.put(sample)
        self.filterqueue.join()

    def filterthreads(self):
        while True:
            sample = self.filterqueue.get()
            # Only run on samples that have been processed with spades
            if os.path.isfile(sample.general.contigsfile) and not os.path.isfile(sample.general.filteredfile):
                # Create a list to store all the records of contigs that pass the minimum depth filtering
                passdepth = list()
                for record in SeqIO.parse(open(sample.general.contigsfile, "rU"), "fasta"):
                    # Extract the values for the mean and standard deviation of the coverage, and split off
                    # the X from the end e.g. 73.07X
                    coveragemean = float(sample.mapping.MeanCoveragedata.split('X')[0])
                    coveragestd = float(sample.mapping.StdCoveragedata.split('X')[0])
                    # Remove the _pilon added to the contig name in order to allow the contig name to match the original
                    # name used as the key in the sample.depth.coverage dictionary
                    contig = record.id.split('_pilon')[0]
                    # Only include contigs with a depth greater or equal to 10
                    if float(sample.depth.coverage[contig]) > (coveragemean - coveragestd * 1.5) \
                            and len(record.seq) > 500:
                        # Replace 'NODE' in the fasta header with the sample name
                        # >NODE_1_length_705814_cov_37.107_ID_4231
                        newid = re.sub("NODE", sample.name, record.id)
                        record.id = str(record.id).replace('NODE', sample.name)
                        record.id = newid
                        # Clear the name and description attributes of the record
                        record.name = ''
                        record.description = ''
                        # Add this record to our list
                        passdepth.append(record)
                # Only create the file if there are contigs that pass the depth filter
                if passdepth:
                    # Open the filtered assembly file
                    with open(sample.general.filteredfile, 'w') as formatted:
                        # Write the records in the list to the file
                        SeqIO.write(passdepth, formatted, 'fasta')
            # If the filtered file was successfully created, copy it to the BestAssemblies folder
            if os.path.isfile(sample.general.filteredfile):
                # Set the assemblies path
                sample.general.bestassembliespath = os.path.join(self.path, 'BestAssemblies')
                # Make the path (if necessary)
                make_path(sample.general.bestassembliespath)
                # Set the name of the file in the best assemblies folder
                bestassemblyfile = os.path.join(sample.general.bestassembliespath, '{}.fasta'.format(sample.name))
                # Add the name and path of the best assembly file to the metadata
                sample.general.bestassemblyfile = bestassemblyfile
                # Copy the filtered file to the BestAssemblies folder
                if not os.path.isfile(bestassemblyfile):
                    shutil.copyfile(sample.general.filteredfile, bestassemblyfile)
            else:
                sample.general.bestassemblyfile = 'NA'
            self.filterqueue.task_done()

    def clear(self):
        """
        Clear out large attributes from the metadata objects
        """
        for sample in self.metadata:
            try:
                delattr(sample.depth, 'bases')
                delattr(sample.depth, 'coverage')
                delattr(sample.depth, 'length')
                delattr(sample.depth, 'stddev')
            except KeyError:
                pass

    @staticmethod
    def analyze(line):
        # Split on ' = '
        if ' = ' in line:
            key, value = line.split(' = ')
            # Replace occurrences of
            key = key.replace('number of ', "").replace("'", "").title().replace(" ", "")
            # Should we keep comma separation?
            value = value.replace(",", "").replace(" ", "").rstrip()
        # Otherwise set the keys and values to None
        else:
            key, value = None, None
        return key, value

    def __init__(self, inputobject):
        from queue import Queue
        self.metadata = inputobject.runmetadata.samples
        self.start = inputobject.starttime
        self.cpus = inputobject.cpus
        self.threadlock = Lock()
        try:
            self.threads = int(self.cpus / len(self.metadata)) if self.cpus / len(self.metadata) > 1 else 1
        except TypeError:
            self.threads = self.cpus
        self.logfile = inputobject.logfile
        self.path = inputobject.path
        self.samversion = get_version(['samtools']).decode('utf-8').split('\n')[2].split()[1]
        # Initialise queues
        self.mapqueue = Queue(maxsize=self.cpus)
        self.qqueue = Queue(maxsize=self.cpus)
        self.bowqueue = Queue(maxsize=self.cpus)
        self.pilonqueue = Queue(maxsize=self.cpus)
        self.indexqueue = Queue(maxsize=self.cpus)
        self.filterqueue = Queue(maxsize=self.cpus)


if __name__ == '__main__':
    class Parser(object):

        def associate(self):
            # Get the sequences in the sequences folder into a list. Note that they must have a file extension that
            # begins with .fa
            self.strains = [fasta for fasta in sorted(glob('{}*.fa*'.format(self.assemblypath)))
                            if '.fastq' not in fasta]
            for strain in self.strains:
                # Extract the name of the strain from the path and file extension
                strainname = os.path.split(strain)[1].split('.')[0]
                # Find the corresponding fastq files for each strain
                fastq = sorted(glob('{}{}*fastq*'.format(self.fastqpath, strainname)))
                # Ensure that fastq files are present for each assembly
                assert fastq, 'Cannot find fastq files for strain {}'.format(strainname)
                # Create the object
                metadata = MetadataObject()
                # Set the .name attribute to be the file name
                metadata.name = strainname
                # Create the .general attribute
                metadata.general = GenObject()
                # Set the path of the assembly file
                metadata.general.bestassembliespath = self.assemblypath
                # Populate the .fastqfiles category of :self.metadata
                metadata.general.trimmedfastqfiles = fastq
                # Create the output directory path
                metadata.general.outputdirectory = os.path.join(self.path, strainname)
                metadata.mapping = GenObject()
                # Append the metadata for each sample to the list of samples
                self.samples.append(metadata)

        def __init__(self):
            from argparse import ArgumentParser
            import multiprocessing
            parser = ArgumentParser(description='Calculates coverage depth by mapping FASTQ reads against assemblies')
            parser.add_argument('-p', '--path',
                                default=os.getcwd(),
                                help='Specify the path of the folder that either contains the files of interest, or'
                                     'will be used to store the outputs')
            parser.add_argument('-a', '--assemblies',
                                help='Path to a folder of assemblies. If not provided, the script will look for .fa'
                                     'or .fasta files in the path')
            parser.add_argument('-f', '--fastq',
                                help='Path to a folder of fastq files. If not provided, the script will look for '
                                     'fastq or .fastq.gz files in the path')
            parser.add_argument('-t', '--threads',
                                help='Number of threads. Default is the number of cores in the system')
            # Get the arguments into an object
            args = parser.parse_args()
            # Define variables from the arguments - there may be a more streamlined way to do this
            # Add trailing slashes to the path variables to ensure consistent formatting (os.path.join)
            self.path = os.path.join(args.path, '')
            self.assemblypath = os.path.join(args.assemblies, '') if args.assemblies else self.path
            self.fastqpath = os.path.join(args.fastq, '') if args.fastq else self.path
            # Use the argument for the number of threads to use, or default to the number of cpus in the system
            self.cpus = args.threads if args.threads else multiprocessing.cpu_count()
            # Initialise variables
            self.strains = []
            self.samples = []
            self.logfile = os.path.join(self.path, 'logfile.txt')

            # Associate the assemblies and fastq files in a metadata object
            self.associate()

    class MetadataInit(object):
        def __init__(self, start):
            # Run the parser
            self.runmetadata = Parser()
            # Get the appropriate variables from the metadata file
            self.path = self.runmetadata.path
            self.assemblypath = self.runmetadata.assemblypath
            self.fastqpath = self.runmetadata.fastqpath
            self.starttime = start
            self.cpus = self.runmetadata.cpus
            self.logfile = self.runmetadata.logfile
            # Run the analyses - the extra set of parentheses is due to using the __call__ method in the class
            QualiMap(self)

    # Run the class
    from time import time
    starttime = time()
    MetadataInit(starttime)
    printtime('Assembly and characterisation complete', starttime)
