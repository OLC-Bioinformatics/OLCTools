#!/usr/bin/env python
from accessoryFunctions.accessoryFunctions import GenObject, make_path, run_subprocess, write_to_logfile
import accessoryFunctions.metadataprinter as metadataprinter
from biotools import bbtools
from Bio.SeqUtils import GC
from Bio import SeqIO
from subprocess import CalledProcessError
from click import progressbar
from queue import Queue
from glob import glob
import threading
import logging
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
        logging.info('Validating FASTQ files')
        with progressbar(self.metadata) as bar:
            for sample in bar:
                # Tiny files can pass the validation tests - ensure that they don't
                size = os.path.getsize(sample.general.fastqfiles[0])
                if size >= 1000000:
                    # Try to run reformat.sh on the reads - on any errors try to run repair.sh
                    try:
                        out, err, cmd = bbtools.validate_reads(forward_in=sample.general.fastqfiles[0],
                                                               returncmd=True)
                        write_to_logfile(out, err, self.logfile, sample.general.logout, sample.general.logerr,
                                         None, None)
                    except CalledProcessError:
                        # Set the file names for the reformatted and repaired files
                        outputfile1 = os.path.join(sample.general.outputdirectory, '{}_reformatted_R1.fastq.gz'
                                                   .format(sample.name))
                        repair_file1 = os.path.join(sample.general.outputdirectory, '{}_repaired_R1.fastq.gz'
                                                    .format(sample.name))
                        if len(sample.general.fastqfiles) == 2:
                            outputfile2 = os.path.join(sample.general.outputdirectory, '{}_reformatted_R2.fastq.gz'
                                                       .format(sample.name))
                            repair_file2 = os.path.join(sample.general.outputdirectory, '{}_repaired_R2.fastq.gz'
                                                        .format(sample.name))
                        else:
                            outputfile2 = str()
                            repair_file2 = str()
                        # Use reformat.sh to repair the reads. If this fails, discard the sample from the analyses
                        try:
                            logging.warning('Errors detected in FASTQ files for sample {sample}. '
                                            'Please check the following files for details {log} {logout} {logerr}. '
                                            'The pipeline will use reformat.sh to attempt to repair issues'
                                            .format(sample=sample.name,
                                                    log=self.logfile,
                                                    logout=sample.general.logout,
                                                    logerr=sample.general.logerr))
                            if not os.path.isfile(outputfile1):
                                # Run reformat.sh
                                out, err, cmd = bbtools.reformat_reads(forward_in=sample.general.fastqfiles[0],
                                                                       forward_out=outputfile1,
                                                                       returncmd=True)
                                write_to_logfile(out, err, self.logfile, sample.general.logout, sample.general.logerr,
                                                 None, None)
                                # Run repair.sh (if necessary)
                                if outputfile2:
                                    out, err, cmd = bbtools.repair_reads(forward_in=outputfile1,
                                                                         forward_out=repair_file1,
                                                                         returncmd=True)
                                    write_to_logfile(out, err, self.logfile, sample.general.logout,
                                                     sample.general.logerr, None, None)
                            # Ensure that the output file(s) exist before declaring this a success
                            if os.path.isfile(outputfile1):
                                # Update the fastqfiles attribute to point to the repaired files
                                sample.general.fastqfiles = [repair_file1, repair_file2] if repair_file2 \
                                    else [outputfile1]
                        except CalledProcessError:
                            # The file(s) can be created even if there is STDERR from reformat.sh
                            if os.path.isfile(outputfile1) and outputfile2:
                                try:
                                    out, err, cmd = bbtools.repair_reads(forward_in=outputfile1,
                                                                         forward_out=repair_file1,
                                                                         returncmd=True)
                                    write_to_logfile(out, err, self.logfile, sample.general.logout,
                                                     sample.general.logerr, None, None)
                                    # Update the fastqfiles attribute to point to the repaired files
                                    sample.general.fastqfiles = [repair_file1, repair_file2] if repair_file2 else \
                                        [repair_file1]
                                except CalledProcessError:
                                    # Write in the logs that there was an error detected in the FASTQ files
                                    write_to_logfile('An error was detected in the FASTQ files for sample {}. '
                                                     'These files will not be processed further'.format(sample.name),
                                                     'An error was detected in the FASTQ files for sample {}. '
                                                     'These files will not be processed further'.format(sample.name),
                                                     self.logfile,
                                                     sample.general.logout,
                                                     sample.general.logerr,
                                                     None,
                                                     None)
                                    # Update metadata objects with error
                                    self.error(sample, 'fastq_error')
                            else:
                                # Write in the logs that there was an error detected in the FASTQ files
                                write_to_logfile('An error was detected in the FASTQ files for sample {}. '
                                                 'These files will not be processed further'.format(sample.name),
                                                 'An error was detected in the FASTQ files for sample {}. '
                                                 'These files will not be processed further'.format(sample.name),
                                                 self.logfile,
                                                 sample.general.logout,
                                                 sample.general.logerr,
                                                 None,
                                                 None)

                                # Update metadata objects with error
                                self.error(sample, 'fastq_error')
                else:
                    # Update metadata objects with error
                    self.error(sample, 'files_too_small')
        # Print the metadata to file
        metadataprinter.MetadataPrinter(self)

    @staticmethod
    def error(sample, message):
        """
        Check to see if the run GenObject exists. If so, update the run.Description to reflect the error
        :param sample: metadata sample object
        :param message: error message to add to the sample.run.Description attribute
        """
        # Set the .fastqfiles attribute to an empty list
        sample.general.fastqfiles = list()
        # Ensure that the run attribute exists
        if GenObject.isattr(sample, 'run'):
            # If the Description attribute exists, overwrite it, otherwise create and populate it
            if GenObject.isattr(sample.run, 'status'):
                sample.run.status = message
            else:
                setattr(sample.run, 'status', message)
        # Otherwise create and populate the attribute
        else:
            setattr(sample, 'run', GenObject())
            sample.run.Description = message

    def fastqcthreader(self, level):
        logging.info('Running quality control on {} fastq files'.format(level))
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
                    reader = 'bunzip2 --stdout'
                # Try except loop to allow for missing samples
                try:
                    fastqfiles = sample.general.trimmedfastqfiles
                except AttributeError:
                    fastqfiles = None
            elif level == 'trimmedcorrected':
                reader = 'gunzip --to-stdout'
                # Try except loop to allow for missing samples
                try:
                    fastqfiles = sample.general.trimmedcorrectedfastqfiles
                except AttributeError:
                    fastqfiles = None
            elif level == 'normalised':
                reader = 'gunzip --to-stdout'
                # Try except loop to allow for missing samples
                try:
                    fastqfiles = sample.general.normalisedreads
                except AttributeError:
                    fastqfiles = None
            elif level == 'merged':
                reader = 'gunzip --to-stdout'
                # Try except loop to allow for missing samples
                try:
                    fastqfiles = [sample.general.mergedreads]
                except AttributeError:
                    fastqfiles = None
            else:
                try:
                    reader = 'cat' if '.gz' not in sample.general.fastqfiles[0] else 'gunzip --to-stdout'
                    fastqfiles = sample.general.fastqfiles
                except IndexError:
                    reader = str()
                    fastqfiles = None
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
        # Wait on the queue until everything has been processed
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
        logging.info("Trimming fastq files")
        # Iterate through strains with fastq files
        with progressbar(self.metadata) as bar:
            for sample in bar:
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
                    # Incorporate read length into the minlength parameter - set it to 50 unless one or more of the
                    # reads has a lower calculated length than 50
                    try:
                        lesser_length = min(int(sample.run.forwardlength), int(sample.run.reverselength))
                    except ValueError:
                        lesser_length = int(sample.run.forwardlength)
                    min_len = 50 if lesser_length >= 50 else lesser_length
                    # Initialise a variable to store the number of bases to automatically trim from the beginning of
                    # each read, as these bases tend to have lower quality scores. If trimming the reads will cause
                    trim_left = 0
                    # If, for some reason, only the reverse reads are present, use the appropriate output file name
                    try:
                        if 'R2' in fastqfiles[0]:
                            if not os.path.isfile(cleanreverse):
                                out, \
                                    err, \
                                    bbdukcall = bbtools.bbduk_trim(forward_in=fastqfiles[0],
                                                                   reverse_in=None,
                                                                   forward_out=cleanreverse,
                                                                   trimq=10,
                                                                   minlength=min_len,
                                                                   forcetrimleft=trim_left,
                                                                   returncmd=True)
                            else:
                                bbdukcall = str()
                                out = str()
                                err = str()
                        else:
                            if not os.path.isfile(cleanforward):
                                    out, \
                                        err, \
                                        bbdukcall = bbtools.bbduk_trim(forward_in=fastqfiles[0],
                                                                       forward_out=cleanforward,
                                                                       trimq=10,
                                                                       minlength=min_len,
                                                                       forcetrimleft=trim_left,
                                                                       returncmd=True)
                            else:
                                bbdukcall = str()
                                out = str()
                                err = str()
                    except (IndexError, CalledProcessError):
                        bbdukcall = str()
                        out = str()
                        err = str()
                    # Write the command, stdout, and stderr to the logfile
                    write_to_logfile(bbdukcall, bbdukcall, self.logfile, sample.general.logout, sample.general.logerr,
                                     None, None)
                    write_to_logfile(out, err, self.logfile, sample.general.logout, sample.general.logerr, None, None)
                    # Add the trimmed fastq files to a list
                    trimmedfastqfiles = sorted(glob(os.path.join(sample.general.outputdirectory, '*trimmed.fastq.gz')))
                    # Populate the metadata if the files exist
                    sample.general.trimmedfastqfiles = trimmedfastqfiles if trimmedfastqfiles else list()
        # Add all the trimmed files to the metadata
        logging.info('Fastq files trimmed')

    def contamination_finder(self, input_path=None, report_path=None):
        """
        Helper function to get confindr integrated into the assembly pipeline
        """
        logging.info('Calculating contamination in reads')
        if input_path is not None:
            input_dir = input_path
        else:
            input_dir = self.path
        if report_path is not None:
            reportpath = report_path
        else:
            reportpath = os.path.join(input_dir, 'confindr')
        confindr_report = os.path.join(input_dir, 'confindr', 'confindr_report.csv')
        pipeline_report = os.path.join(reportpath, 'confindr_report.csv')
        # Only proceed if the confindr report doesn't exist
        if not os.path.isfile(confindr_report):
            # # Create an object to store attributes to pass to confinder
            # Clear and recreate the output folder
            try:
                shutil.rmtree(reportpath)
            except IOError:
                pass
            make_path(reportpath)
            # Run confindr
            systemcall = 'confindr.py -i {input_dir} -o {output_dir} -d {database_dir} -bf 0.05'\
                .format(input_dir=input_dir,
                        output_dir=os.path.join(input_dir, 'confindr'),
                        database_dir=os.path.join(self.reffilepath, 'ConFindr', 'databases'))
            # Run the call
            out, err = run_subprocess(systemcall)
            write_to_logfile(systemcall, systemcall, self.logfile, None, None, None, None)
            write_to_logfile(out, err, self.logfile, None, None, None, None)
            logging.info('Contamination detection complete!')
        # Load the confindr report into a dictionary using pandas
        # https://stackoverflow.com/questions/33620982/reading-csv-file-as-dictionary-using-pandas
        confindr_results = pandas.read_csv(confindr_report, index_col=0).T.to_dict()
        # Find the results for each of the samples
        for sample in self.metadata:
            # Create a GenObject to store the results
            sample.confindr = GenObject()
            # Iterate through the dictionary to find the outputs for each sample
            for line in confindr_results:
                # If the current line corresponds to the sample of interest
                if sample.name in line:
                    # Set the values using the appropriate keys as the attributes
                    sample.confindr.genus = confindr_results[line]['Genus'] if type(confindr_results[line]['Genus']) \
                                                                               is not float else 'ND'
                    sample.confindr.num_contaminated_snvs = confindr_results[line]['NumContamSNVs']
                    sample.confindr.contam_status = confindr_results[line]['ContamStatus']
                    # Don't break parsing previous ConFindr reports that lack the percent contamination calculations
                    try:
                        sample.confindr.percent_contam = confindr_results[line]['PercentContam'] if \
                            str(confindr_results[line]['PercentContam']) != 'nan' else 0
                    except KeyError:
                        sample.confindr.percent_contam = 'ND'
                    try:
                        sample.confindr.percent_contam_std = \
                            confindr_results[line]['PercentContamStandardDeviation'] if \
                            str(confindr_results[line]['PercentContamStandardDeviation']) != 'nan' else 0
                    except KeyError:
                        sample.confindr.percent_contam_std = 'ND'
                    if sample.confindr.contam_status is True:
                        sample.confindr.contam_status = 'Contaminated'
                    elif sample.confindr.contam_status is False:
                        sample.confindr.contam_status = 'Clean'
        # Re-write the output to be consistent with the rest of the pipeline
        with open(pipeline_report, 'w') as csv:
            data = 'Strain,Genus,NumContamSNVs,ContamStatus,PercentContam,PercentContamSTD\n'
            for sample in self.metadata:
                data += '{str},{genus},{numcontamsnv},{status},{pc},{pcs}\n'.format(
                    str=sample.name,
                    genus=sample.confindr.genus,
                    numcontamsnv=sample.confindr.num_contaminated_snvs,
                    status=sample.confindr.contam_status,
                    pc=sample.confindr.percent_contam,
                    pcs=sample.confindr.percent_contam_std
                )
            csv.write(data)

    def estimate_genome_size(self):
        """
        Use kmercountexact from the bbmap suite of tools to estimate the size of the genome
        """
        logging.info('Estimating genome size using kmercountexact')
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
        logging.info('Error correcting reads')
        for sample in self.metadata:
            sample.general.trimmedcorrectedfastqfiles = [fastq.split('.fastq.gz')[0] + '_trimmed_corrected.fastq.gz'
                                                         for fastq in sorted(sample.general.fastqfiles)]
            try:
                if not os.path.isfile(sample.general.trimmedcorrectedfastqfiles[0]):
                    try:
                        out, err, cmd = bbtools.tadpole(forward_in=sorted(sample.general.trimmedfastqfiles)[0],
                                                        forward_out=sample.general.trimmedcorrectedfastqfiles[0],
                                                        returncmd=True,
                                                        mode='correct',
                                                        threads=self.cpus)
                        # Set the command in the object
                        sample[self.analysistype].errorcorrectcmd = cmd
                        write_to_logfile(out=out,
                                         err=err,
                                         logfile=self.logfile,
                                         samplelog=sample.general.logout,
                                         sampleerr=sample.general.logerr,
                                         analysislog=None,
                                         analysiserr=None)
                    except IndexError:
                        sample.general.trimmedcorrectedfastqfiles = list()
            except CalledProcessError:
                sample.general.trimmedcorrectedfastqfiles = sample.general.trimmedfastqfiles
            except AttributeError:
                sample.general.trimmedcorrectedfastqfiles = list()
            except IndexError:
                sample.general.trimmedcorrectedfastqfiles = list()

    def normalise_reads(self):
        """
        Use bbnorm from the bbmap suite of tools to perform read normalisation
        """
        logging.info('Normalising reads to a kmer depth of 100')
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
        logging.info('Merging paired reads')
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
                    # Run the merging command - forward_in=sample.general.normalisedreads[0],
                    out, err, cmd = bbtools.bbmerge(forward_in=sorted(sample.general.trimmedcorrectedfastqfiles)[0],
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
        self.perform_pilon()
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

    def perform_pilon(self):
        """
        Determine if pilon polishing should be attempted. Do not perform polishing if confindr determines that the
        sample is contaminated or if there are > 500 contigs
        """
        for sample in self.metadata:
            try:
                if sample[self.analysistype].num_contigs > 500 or sample.confindr.contam_status == 'Contaminated':
                    sample.general.polish = False
                else:
                    sample.general.polish = True
            except AttributeError:
                sample.general.polish = True

    def clear_attributes(self):
        """
        Remove the record_dict attribute from the object, as SeqRecords are not JSON-serializable. Also remove
        the contig_lengths and longest_contig attributes, as they are large lists that make the .json file ugly
        """
        for sample in self.metadata:
            try:
                delattr(sample[self.analysistype], 'record_dict')
                delattr(sample[self.analysistype], 'contig_lengths')
                delattr(sample[self.analysistype], 'longest_contig')
            except AttributeError:
                pass

    def __init__(self, inputobject, analysis):
        self.metadata = inputobject.runmetadata.samples
        self.start = inputobject.starttime
        self.analysistype = 'quality_features_{analysis}'.format(analysis=analysis)


class GenomeQAML(object):

    def main(self):
        """
        Run the methods
        """
        self.run_qaml()
        self.parse_qaml()

    def run_qaml(self):
        """
        Create and run the GenomeQAML system call
        """
        logging.info('Running GenomeQAML quality assessment')
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
        logging.info('Parsing GenomeQAML outputs')
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
