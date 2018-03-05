#!/usr/bin/env python
from accessoryFunctions.accessoryFunctions import GenObject, MetadataObject, printtime, make_path, \
    run_subprocess, write_to_logfile
try:
    from confindr import confindr
except ImportError:
    import confindr
from biotools import bbtools
from Bio.SeqUtils import GC
from Bio import SeqIO
from subprocess import CalledProcessError
from queue import Queue
from glob import glob
import threading
import pandas
import shutil
import os
__author__ = 'adamkoziol'


class Quality(object):

    def validate_fastq(self):
        """
        Runs reformat.sh on the FASTQ files. If a CalledProcessError arises, do not proceed with the assembly of
        these files
        """
        printtime('Validating FASTQ files', self.start)
        validated_reads = list()
        for sample in self.metadata:
            # Try to run reformat.sh on the reads - on any errors try to run repair.sh
            try:
                out, err, cmd = bbtools.validate_reads(forward_in=sample.general.fastqfiles[0],
                                                       returncmd=True)
                write_to_logfile(out, err, self.logfile, sample.general.logout, sample.general.logerr, None, None)
                # Add the sample to the list of samples with FASTQ files that pass this validation step
                validated_reads.append(sample)
            except CalledProcessError:
                # Try to use reformat.sh to repair the reads - if this fails, discard the sample from the analyses
                try:
                    printtime('Errors detected in FASTQ files for sample {sample}. Please check the following log files'
                              'for details {log} {logout} {logerr}. Using reformat.sh to attempt to repair issues'
                              .format(sample=sample.name,
                                      log=self.logfile,
                                      logout=sample.general.logout,
                                      logerr=sample.general.logerr), self.start)
                    # Set the file names for the repaired files
                    outputfile1 = os.path.join(sample.general.outputdirectory, '{}_repaired_R1.fastq.gz'
                                               .format(sample.name))
                    outputfile2 = os.path.join(sample.general.outputdirectory, '{}_repaired_R2.fastq.gz'
                                               .format(sample.name))
                    # Run repair.sh
                    out, err, cmd = bbtools.repair_reads(forward_in=sample.general.fastqfiles[0],
                                                         forward_out=outputfile1,
                                                         returncmd=True)
                    write_to_logfile(out, err, self.logfile, sample.general.logout, sample.general.logerr, None, None)
                    # Update the fastqfiles attribute to point to the repaired files
                    sample.general.fastqfiles = [outputfile1, outputfile2]
                    # Add the sample object to the list of samples passing the FASTQ validation step
                    validated_reads.append(sample)
                except CalledProcessError:
                    # Write that there was an error detected in the FASTQ files
                    write_to_logfile('An error was detected in FASTQ files. These files will not be processed further',
                                     'An error was detected in FASTQ files. These files will not be processed further',
                                     self.logfile,
                                     sample.general.logout,
                                     sample.general.logerr, None, None)
                    # Set the .fastqfiles attribute to 'NA' to remove this strain from the analyses
                    sample.general.fastqfiles = 'NA'
        # Overwrite self.metadata with objects that do not fail the validation
        self.metadata = validated_reads

    def fastqcthreader(self, level):
        printtime('Running quality control on {} fastq files'.format(level), self.start)
        for sample in self.metadata:
            if type(sample.general.fastqfiles) is list:
                # Create and start threads for each fasta file in the list
                threads = threading.Thread(target=self.fastqc, args=())
                # Set the daemon to true - something to do with thread management
                threads.setDaemon(True)
                # Start the threading
                threads.start()
        # Iterate through strains with fastq files to set variables to add to the multithreading queue
        for sample in self.metadata:
            fastqccall = str()
            fastqcreads = str()
            # Check to see if the fastq files exist
            if level == 'Trimmed':
                # Set the appropriate method to read in the trimmed fastq files - if uncompressed, use cat, if
                # compressed with bzip, use bunzip2, and if compressed with gzip use gunzip
                if '.bz2' or '.gz' not in sample.general.trimmedfastqfiles[0]:
                    reader = 'cat'
                elif '.gz' in sample.general.trimmedfastqfiles[0]:
                    reader = 'gunzip --to-stdout'
                else:
                    'bunzip2 --stdout'
                # Try except loop to allow for missing samples
                try:
                    fastqfiles = sample.general.trimmedfastqfiles
                except KeyError:
                    fastqfiles = ""
                    pass
            elif level == 'trimmedcorrected':
                reader = 'gunzip --to-stdout'
                # Try except loop to allow for missing samples
                try:
                    fastqfiles = sample.general.trimmedcorrectedfastqfiles
                except KeyError:
                    fastqfiles = ""
                    pass
            elif level == 'normalised':
                reader = 'gunzip --to-stdout'
                # Try except loop to allow for missing samples
                try:
                    fastqfiles = sample.general.normalisedreads
                except KeyError:
                    fastqfiles = ""
                    pass
            elif level == 'merged':
                reader = 'gunzip --to-stdout'
                # Try except loop to allow for missing samples
                try:
                    fastqfiles = [sample.general.mergedreads]
                except KeyError:
                    fastqfiles = ""
                    pass
            else:
                reader = 'cat' if '.gz' not in sample.general.fastqfiles[0] else 'gunzip --to-stdout'
                fastqfiles = sample.general.fastqfiles
            # As the metadata can be populated with 'NA' (string) if there are no fastq files, only process if
            # :fastqfiles is a list
            if type(fastqfiles) is list:
                # Set the output directory location
                outdir = os.path.join(sample.general.outputdirectory, 'fastqc', level)
                make_path(outdir)
                # Separate system calls for paired and unpaired fastq files
                if len(fastqfiles) == 2:
                    # Combine the fastq files, so that the paired files are processed together
                    # Call fastqc with -q (quiet), -o (output directory), -t (number of threads) flags
                    fastqccall = '{} {} {} | fastqc -q -t {} stdin -o {}'\
                        .format(reader, fastqfiles[0], fastqfiles[1], self.threads, outdir)
                    # Also perform QC on the individual forward and reverse reads
                    fastqcreads = "fastqc {} {} -q -o {} -t {}"\
                        .format(fastqfiles[0], fastqfiles[1], outdir, self.threads)
                elif len(fastqfiles) == 1:
                    fastqccall = '{} {} | fastqc -q -t {} stdin -o {}'\
                        .format(reader, fastqfiles[0], self.threads, outdir)
                    fastqcreads = "fastqc {} -q -o {} -t {}".format(fastqfiles[0], outdir, self.threads)
                # Add the arguments to the queue
                sample.commands.fastqc = fastqccall
                self.qcqueue.put((sample, fastqccall, outdir, fastqcreads))
        # Wait on the trimqueue until everything has been processed
        self.qcqueue.join()

    def fastqc(self):
        """Run fastqc system calls"""
        while True:  # while daemon
            threadlock = threading.Lock()
            # Unpack the variables from the queue
            (sample, systemcall, outputdir, fastqcreads) = self.qcqueue.get()
            # Check to see if the output HTML file already exists
            try:
                _ = glob(os.path.join(outputdir, '*.html'))[0]
            except IndexError:
                # Make the output directory
                make_path(outputdir)
                # Run the system calls
                outstr = str()
                errstr = str()
                out, err = run_subprocess(systemcall)
                outstr += out
                errstr += err
                out, err = run_subprocess(fastqcreads)
                outstr += out
                errstr += err
                # Acquire thread lock, and write the logs to file
                threadlock.acquire()
                write_to_logfile(systemcall, systemcall, self.logfile, sample.general.logout, sample.general.logerr,
                                 None, None)
                write_to_logfile(fastqcreads, fastqcreads, self.logfile, sample.general.logout, sample.general.logerr,
                                 None, None)
                write_to_logfile(outstr, errstr, self.logfile, sample.general.logout, sample.general.logerr, None, None)
                threadlock.release()
                # Rename the outputs
                try:
                    shutil.move(os.path.join(outputdir, 'stdin_fastqc.html'),
                                os.path.join(outputdir, '{}_fastqc.html'.format(sample.name)))
                    shutil.move(os.path.join(outputdir, 'stdin_fastqc.zip'),
                                os.path.join(outputdir, '{}_fastqc.zip'.format(sample.name)))
                except IOError:
                    pass
            # Signal to qcqueue that job is done
            self.qcqueue.task_done()

    def trimquality(self):
        """Uses bbduk from the bbmap tool suite to quality and adapter trim"""
        printtime("Trimming fastq files", self.start)
        # Create and start threads for each strain with fastq files
        for sample in self.metadata:
            if type(sample.general.fastqfiles) is list:
                # Create and start threads for each fasta file in the list
                threads = threading.Thread(target=self.bbduker, args=())
                # Set the daemon to true - something to do with thread management
                threads.setDaemon(True)
                # Start the threading
                threads.start()
        # Iterate through strains with fastq files to set variables to add to the multithreading queue
        for sample in self.metadata:
            # As the metadata can be populated with 'NA' (string) if there are no fastq files, only process if
            # :fastqfiles is a list
            if type(sample.general.fastqfiles) is list:
                # Check to see if the fastq files exist
                fastqfiles = sorted(sample.general.fastqfiles)
                # Define the output directory
                outputdir = sample.general.outputdirectory
                # Define the name of the trimmed fastq files
                cleanforward = os.path.join(outputdir, '{}_R1_trimmed.fastq.gz'.format(sample.name))
                cleanreverse = os.path.join(outputdir, '{}_R2_trimmed.fastq.gz'.format(sample.name))
                min_len = 50
                if self.numreads == 2:
                    # Separate system calls for paired and unpaired fastq files
                    # http://seqanswers.com/forums/showthread.php?t=42776
                    # BBduk 37.23 doesn't need the ktrim=l/mink=11 parameters, so they have been removed.
                    if len(fastqfiles) == 2:
                        # Incorporate read length into the minlength parameter - set it to 50 unless one or more of the
                        # reads has a lower calculated length than 50
                        try:
                            lesser_length = min(int(sample.run.forwardlength), int(sample.run.reverselength))
                            min_len = 50 if lesser_length >= 50 else lesser_length
                            bbdukcall = "bbduk.sh -Xmx1g in1={in1} in2={in2} out1={out1} out2={out2} qtrim=w " \
                                        "trimq=10 k=25 minlength={ml} ref=adapters tbo" \
                                .format(in1=fastqfiles[0],
                                        in2=fastqfiles[1],
                                        out1=cleanforward,
                                        out2=cleanreverse,
                                        ml=min_len)
                        except ValueError:
                            bbdukcall = ''

                    elif len(fastqfiles) == 1:
                        try:
                            lesser_length = int(sample.run.forwardlength)
                        except ValueError:
                            lesser_length = int(sample.run.reverselength)
                        min_len = 50 if lesser_length >= 50 else lesser_length
                        bbdukcall = "bbduk.sh -Xmx1g in={in1} out={out1} qtrim=w trimq=10 k=25 minlength={ml} " \
                                    "ref=adapters" \
                            .format(in1=fastqfiles[0],
                                    out1=cleanforward,
                                    ml=min_len)
                    else:
                        bbdukcall = ""
                # Allows for exclusion of the reverse reads if desired
                else:
                    bbdukcall = "bbduk.sh -Xmx1g in={in1} out={out1} qtrim=w trimq=10 k=25 minlength={ml} " \
                                "ref=adapters"\
                        .format(in1=fastqfiles[0],
                                out1=cleanforward,
                                ml=min_len)
                    # There is a check to ensure that the trimmed reverse file is created. This will change the file
                    # being looked for to the forward file
                    cleanreverse = cleanforward
                    if self.forwardlength != 'full':
                        bbdukcall += ' forcetrimright={}'.format(str(self.forwardlength))
                sample.commands.bbduk = bbdukcall
                # Add the arguments to the queue
                self.trimqueue.put((sample, bbdukcall, cleanreverse))
        # Wait on the trimqueue until everything has been processed
        self.trimqueue.join()
        # Add all the trimmed files to the metadata
        printtime('Fastq files trimmed', self.start)

    def bbduker(self):
        """Run bbduk system calls"""
        while True:  # while daemon
            # Unpack the variables from the queue
            (sample, systemcall, reversename) = self.trimqueue.get()
            # Check to see if the forward file already exists
            if systemcall:
                threadlock = threading.Lock()
                if not os.path.isfile(reversename) and not os.path.isfile('{}.bz2'.format(reversename)):
                    # Run the call
                    out, err = run_subprocess(systemcall)
                    threadlock.acquire()
                    write_to_logfile(systemcall, systemcall, self.logfile, sample.general.logout, sample.general.logerr,
                                     None, None)
                    write_to_logfile(out, err, self.logfile, sample.general.logout, sample.general.logerr, None, None)
                    threadlock.release()
                # Define the output directory
                outputdir = sample.general.outputdirectory
                # Add the trimmed fastq files to a list
                trimmedfastqfiles = sorted(glob(os.path.join(outputdir, '*trimmed.fastq.gz')))
                # Populate the metadata if the files exist
                sample.general.trimmedfastqfiles = trimmedfastqfiles if trimmedfastqfiles else 'NA'
            # Signal to trimqueue that job is done
            self.trimqueue.task_done()

    def contamination_finder(self):
        """
        Helper function to get confindr integrated into the assembly pipeline
        """
        printtime('Calculating contamination in reads', self.start)
        reportpath = os.path.join(self.path, 'confindr')
        report = os.path.join(reportpath, 'confindr_report.csv')
        if not os.path.isfile(report):
            # Create an object to store attributes to pass to confinder
            args = MetadataObject
            args.input_directory = self.path
            args.output_name = reportpath
            args.databases = os.path.join(self.reffilepath, 'ConFindr', 'databases')
            args.forward_id = '_R1'
            args.reverse_id = '_R2'
            args.threads = self.cpus
            args.kmer_size = 31
            args.number_subsamples = 3
            args.subsample_depth = 20
            args.kmer_cutoff = 2
            try:
                shutil.rmtree(args.output_name)
            except IOError:
                pass
            make_path(reportpath)
            # Open the output report file.
            with open(os.path.join(report), 'w') as f:
                f.write('Sample,Genus,NumContamSNVs,NumUniqueKmers,ContamStatus\n')
            for sample in self.metadata:
                if len(sample.general.trimmedcorrectedfastqfiles) == 2:
                    confindr.find_contamination(sample.general.trimmedcorrectedfastqfiles, args)
                elif len(sample.general.trimmedcorrectedfastqfiles) == 1:
                    confindr.find_contamination_unpaired(args, sample.general.trimmedcorrectedfastqfiles[0])
            printtime('Contamination detection complete!', self.start)
        # Load the confindr report into a dictionary using pandas
        # https://stackoverflow.com/questions/33620982/reading-csv-file-as-dictionary-using-pandas
        confindr_results = pandas.read_csv(report, index_col=0).T.to_dict()
        # Find the results for each of the samples
        for sample in self.metadata:
            # Create a GenObject to store the results
            sample.confinder = GenObject()
            # Iterate through the dictionary to find the outputs for each sample
            for line in confindr_results:
                # If the current line corresponds to the sample of interest
                if sample.name in line:
                    # Set the values using the appropriate keys as the attributes
                    sample.confinder.genus = confindr_results[line]['Genus']
                    sample.confinder.num_contaminated_snvs = confindr_results[line]['NumContamSNVs']
                    sample.confinder.unique_kmers = confindr_results[line]['NumUniqueKmers']
                    try:
                        sample.confinder.cross_contamination = confindr_results[line]['CrossContamination']
                    except KeyError:
                        sample.confinder.cross_contamination = str()
                    sample.confinder.contam_status = confindr_results[line]['ContamStatus']
                    if sample.confinder.contam_status is True:
                        sample.confinder.contam_status = 'Contaminated'
                    elif sample.confinder.contam_status is False:
                        sample.confinder.contam_status = 'Clean'
        # Copy the report to the folder containing all reports for the pipeline
        try:
            shutil.copyfile(report, os.path.join(self.path, 'reports', 'confindr_report.csv'))
        except IOError:
            pass

    def estimate_genome_size(self):
        """
        Use kmercountexact from the bbmap suite of tools to estimate the size of the genome
        """
        printtime('Estimating genome size using kmercountexact', self.start)
        for sample in self.metadata:
            # Initialise the name of the output file
            sample[self.analysistype].peaksfile = os.path.join(sample[self.analysistype].outputdir, 'peaks.txt')
            # Run the kmer counting command
            out, err, cmd = bbtools.kmercountexact(forward_in=sorted(sample.general.fastqfiles)[0],
                                                   peaks=sample[self.analysistype].peaksfile,
                                                   returncmd=True,
                                                   threads=self.cpus)
            # Set the command in the object
            sample[self.analysistype].kmercountexactcmd = cmd
            # Extract the genome size from the peaks file
            sample[self.analysistype].genomesize = bbtools.genome_size(sample[self.analysistype].peaksfile)
            write_to_logfile(out, err, self.logfile, sample.general.logout, sample.general.logerr, None, None)

    def error_correction(self):
        """
        Use tadpole from the bbmap suite of tools to perform error correction of the reads
        """
        printtime('Error correcting reads', self.start)
        for sample in self.metadata:
            sample.general.trimmedcorrectedfastqfiles = [fastq.split('.fastq.gz')[0] + '_trimmed_corrected.fastq.gz'
                                                         for fastq in sorted(sample.general.fastqfiles)]
            try:
                out, err, cmd = bbtools.tadpole(forward_in=sorted(sample.general.trimmedfastqfiles)[0],
                                                forward_out=sample.general.trimmedcorrectedfastqfiles[0],
                                                returncmd=True,
                                                mode='correct',
                                                threads=self.cpus)
                # Set the command in the object
                sample[self.analysistype].errorcorrectcmd = cmd
                write_to_logfile(out, err, self.logfile, sample.general.logout, sample.general.logerr, None, None)
            except CalledProcessError:
                sample.general.trimmedcorrectedfastqfiles = sample.general.trimmedfastqfiles
            except KeyError:
                sample.general.trimmedcorrectedfastqfiles = list()

    def normalise_reads(self):
        """
        Use bbnorm from the bbmap suite of tools to perform read normalisation
        """
        printtime('Normalising reads to a kmer depth of 100', self.start)
        for sample in self.metadata:
            # Set the name of the normalised read files
            sample.general.normalisedreads = [fastq.split('.fastq.gz')[0] + '_normalised.fastq.gz'
                                              for fastq in sorted(sample.general.fastqfiles)]
            try:
                # Run the normalisation command
                out, err, cmd = bbtools.bbnorm(forward_in=sorted(sample.general.trimmedcorrectedfastqfiles)[0],
                                               forward_out=sample.general.normalisedreads[0],
                                               returncmd=True,
                                               threads=self.cpus)
                sample[self.analysistype].normalisecmd = cmd
                write_to_logfile(out, err, self.logfile, sample.general.logout, sample.general.logerr, None, None)
            except CalledProcessError:
                sample.general.normalisedreads = sample.general.trimmedfastqfiles
            except IndexError:
                sample.general.normalisedreads = list()

    def merge_pairs(self):
        """
        Use bbmerge from the bbmap suite of tools to merge paired-end reads
        """
        printtime('Merging paired reads', self.start)
        for sample in self.metadata:
            # Can only merge paired-end
            if len(sample.general.fastqfiles) == 2:
                # Set the name of the merged, and unmerged files
                sample.general.mergedreads = \
                    os.path.join(sample.general.outputdirectory, '{}_paired.fastq.gz'.format(sample.name))
                sample.general.unmergedforward = \
                    os.path.join(sample.general.outputdirectory, '{}_unpaired_R1.fastq.gz'.format(sample.name))
                sample.general.unmergedreverse = \
                    os.path.join(sample.general.outputdirectory, '{}_unpaired_R2.fastq.gz'.format(sample.name))
                try:
                    # Run the merging command
                    out, err, cmd = bbtools.bbmerge(forward_in=sample.general.normalisedreads[0],
                                                    merged_reads=sample.general.mergedreads,
                                                    returncmd=True,
                                                    outu1=sample.general.unmergedforward,
                                                    outu2=sample.general.unmergedreverse,
                                                    threads=self.cpus)
                    sample[self.analysistype].bbmergecmd = cmd
                    write_to_logfile(out, err, self.logfile, sample.general.logout, sample.general.logerr, None, None)
                except CalledProcessError:
                    delattr(sample.general, 'mergedreads')
                    delattr(sample.general, 'unmergedforward')
                    delattr(sample.general, 'unmergedreverse')
                except IndexError:
                    delattr(sample.general, 'mergedreads')
                    delattr(sample.general, 'unmergedforward')
                    delattr(sample.general, 'unmergedreverse')
            else:
                sample.general.mergedreads = sorted(sample.general.trimmedcorrectedfastqfiles)[0]

    def __init__(self, inputobject):
        self.metadata = inputobject.runmetadata.samples
        self.cpus = inputobject.cpus
        try:
            self.threads = int(self.cpus / len(self.metadata)) if self.cpus / len(self.metadata) > 1 else 1
        except TypeError:
            self.threads = self.cpus
        # self.devnull = open(os.devnull, 'wb')
        self.qcqueue = Queue(maxsize=self.cpus)
        self.trimqueue = Queue(maxsize=self.cpus)
        self.correctqueue = Queue(maxsize=self.cpus)
        self.start = inputobject.starttime
        try:
            self.forwardlength = inputobject.forwardlength
            self.reverselength = inputobject.reverselength
        except AttributeError:
            self.forwardlength = 'full'
            self.reverselength = 'full'
        self.numreads = inputobject.numreads
        self.logfile = inputobject.logfile
        self.path = inputobject.path
        self.analysistype = 'quality'
        self.reffilepath = inputobject.reffilepath
        # Initialise the quality attribute in the metadata object
        for sample in self.metadata:
            setattr(sample, self.analysistype, GenObject())


class QualityFeatures(object):

    def main(self):
        """
        Run all the methods required for pipeline outputs
        """
        self.fasta_records()
        self.fasta_stats()
        self.find_largest_contig()
        self.find_genome_length()
        self.find_num_contigs()
        self.find_n50()
        self.clear_attributes()

    def fasta_records(self):
        """
        Use SeqIO to create dictionaries of all records for each FASTA file
        """
        for sample in self.metadata:
            # Create the analysis-type specific attribute
            setattr(sample, self.analysistype, GenObject())
            # Create a dictionary of records for each file
            try:
                record_dict = SeqIO.to_dict(SeqIO.parse(sample.general.bestassemblyfile, "fasta"))
            except FileNotFoundError:
                record_dict = dict()
            # Set the records dictionary as the attribute for the object
            sample[self.analysistype].record_dict = record_dict

    def fasta_stats(self):
        """
        Parse the lengths of all contigs for each sample, as well as the total GC%
        """
        for sample in self.metadata:
            # Initialise variables to store appropriate values parsed from contig records
            contig_lengths = list()
            fasta_sequence = str()
            for contig, record in sample[self.analysistype].record_dict.items():
                # Append the length of the contig to the list
                contig_lengths.append(len(record.seq))
                # Add the contig sequence to the string
                fasta_sequence += record.seq
            # Set the reverse sorted (e.g. largest to smallest) list of contig sizes as the value
            sample[self.analysistype].contig_lengths = sorted(contig_lengths, reverse=True)
            try:
                # Calculate the GC% of the total genome sequence using GC - format to have two decimal places
                sample[self.analysistype].gc = float('{:0.2f}'.format(GC(fasta_sequence)))
            except TypeError:
                sample[self.analysistype].gc = 'NA'

    def find_largest_contig(self):
        """
        Determine the largest contig for each strain
        """
        # for file_name, contig_lengths in contig_lengths_dict.items():
        for sample in self.metadata:
            # As the list is sorted in descending order, the largest contig is the first entry in the list
            sample[self.analysistype].longest_contig = sample[self.analysistype].contig_lengths

    def find_genome_length(self):
        """
        Determine the total length of all the contigs for each strain
        """
        for sample in self.metadata:
            # Use the sum() method to add all the contig lengths in the list
            sample[self.analysistype].genome_length = sum(sample[self.analysistype].contig_lengths)

    def find_num_contigs(self):
        """
        Count the total number of contigs for each strain
        """
        for sample in self.metadata:
            # Use the len() method to count the number of entries in the list
            sample[self.analysistype].num_contigs = len(sample[self.analysistype].contig_lengths)

    def find_n50(self):
        """
        Calculate the N50 for each strain. N50 is defined as the largest contig such that at least half of the total
        genome size is contained in contigs equal to or larger than this contig
        """
        for sample in self.metadata:
            # Initialise the N50 attribute in case there is no assembly, and the attribute is not created in the loop
            sample[self.analysistype].n50 = '-'
            # Initialise a variable to store a running total of contig lengths
            currentlength = 0
            for contig_length in sample[self.analysistype].contig_lengths:
                # Increment the current length with the length of the current contig
                currentlength += contig_length
                # If the current length is now greater than the total genome / 2, the current contig length is the N50
                if currentlength >= sample[self.analysistype].genome_length * 0.5:
                    # Populate the dictionary, and break the loop
                    sample[self.analysistype].n50 = contig_length
                    break

    def clear_attributes(self):
        """
        Remove the record_dict attribute from the object, as SeqRecords are not JSON-serializable
        """
        for sample in self.metadata:
            try:
                delattr(sample[self.analysistype], 'record_dict')
            except KeyError:
                pass

    def __init__(self, inputobject):
        self.metadata = inputobject.runmetadata.samples
        self.start = inputobject.starttime
        self.analysistype = 'quality_features'


class GenomeQAML(object):

    def main(self):
        """

        """
        self.run_qaml()
        self.parse_qaml()

    def run_qaml(self):
        """
        Create and run the GenomeQAML system call
        """
        printtime('Running GenomeQAML quality assessment', self.start)
        qaml_call = 'classify.py -t {tf} -r {rf}'\
            .format(tf=self.qaml_path,
                    rf=self.qaml_report)
        make_path(self.reportpath)
        # Only attempt to assess assemblies if the report doesn't already exist
        if not os.path.isfile(self.qaml_report):
            # Run the system calls
            out, err = run_subprocess(qaml_call)
            # Acquire thread lock, and write the logs to file
            self.threadlock.acquire()
            write_to_logfile(qaml_call, qaml_call, self.logfile)
            write_to_logfile(out, err, self.logfile)
            self.threadlock.release()

    def parse_qaml(self):
        """
        Parse the GenomeQAML report, and populate metadata objects
        """
        printtime('Parsing GenomeQAML outputs', self.start)
        # A dictionary to store the parsed excel file in a more readable format
        nesteddictionary = dict()
        # Use pandas to read in the CSV file, and convert the pandas data frame to a dictionary (.to_dict())
        dictionary = pandas.read_csv(self.qaml_report).to_dict()
        # Iterate through the dictionary - each header from the CSV file
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
        # Get the results into the metadata object
        for sample in self.metadata:
            # Initialise the plasmid extractor genobject
            setattr(sample, self.analysistype, GenObject())
            # Initialise the list of all plasmids
            sample[self.analysistype].prediction = str()
            # Iterate through the dictionary of results
            for line in nesteddictionary:
                # Extract the sample name from the dictionary
                name = nesteddictionary[line]['Sample']
                # Ensure that the names match
                if name == sample.name:
                    # Append the plasmid name extracted from the dictionary to the list of plasmids
                    sample[self.analysistype].prediction = nesteddictionary[line]['Predicted_Class']

    def __init__(self, inputobject):
        self.path = inputobject.path
        self.reportpath = inputobject.reportpath
        self.analysistype = 'GenomeQAML'
        self.qaml_path = os.path.join(self.path, 'BestAssemblies')
        self.qaml_report = os.path.join(self.reportpath, 'QAMLReport.csv')
        self.start = inputobject.starttime
        self.cpus = inputobject.cpus
        self.logfile = inputobject.logfile
        self.threadlock = threading.Lock()
        self.metadata = inputobject.runmetadata.samples
