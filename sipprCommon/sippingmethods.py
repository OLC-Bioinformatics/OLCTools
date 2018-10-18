#!/usr/bin/env python
from accessoryFunctions.accessoryFunctions import combinetargets, GenObject, make_path, logstr, \
    write_to_logfile, run_subprocess
from accessoryFunctions.metadataprinter import MetadataPrinter
from sipprCommon.bowtie import Bowtie2CommandLine, Bowtie2BuildCommandLine
import sipprCommon.editsamheaders
from Bio.Sequencing.Applications import SamtoolsFaidxCommandline, SamtoolsIndexCommandline, \
    SamtoolsSortCommandline, SamtoolsViewCommandline
from Bio.Application import ApplicationError
from Bio import SeqIO
from collections import Counter
from click import progressbar
from threading import Thread
from io import StringIO
from queue import Queue
from glob import glob
import logging
import psutil
import numpy
import pysam
import os

__author__ = 'adamkoziol'


class Sippr(object):

    def main(self):
        """
        Run the methods in the correct order for pipelines
        """
        # Find the target files
        self.targets()
        # Use bbduk to bait the FASTQ reads matching the target sequences
        self.bait()
        # If desired, use bbduk to bait the target sequences with the previously baited FASTQ files
        if self.revbait:
            self.reversebait()
        # Run the bowtie2 read mapping module
        self.mapping()
        # Use samtools to index the sorted bam file
        self.indexing()
        # Parse the results
        self.parsing()
        # Filter out any sequences with cigar features such as internal soft-clipping from the results
        self.clipper()

    def targets(self):
        """
        Search the targets folder for FASTA files, create the multi-FASTA file of all targets if necessary, and
        populate objects
        """
        logging.info('Performing analysis with {} targets folder'.format(self.analysistype))
        if self.pipeline:
            for sample in self.runmetadata:
                setattr(sample, self.analysistype, GenObject())
                if sample.general.bestassemblyfile != 'NA':
                    sample[self.analysistype].runanalysis = True
                    # Set attributes
                    try:
                        sample[self.analysistype].targetpath = \
                            os.path.join(self.targetpath, self.analysistype, sample.mash.closestrefseqgenus, '')
                    except AttributeError:
                        sample[self.analysistype].targetpath = \
                            os.path.join(self.targetpath, self.analysistype, sample.general.closestrefseqgenus, '')
                    # There is a relatively strict databasing scheme necessary for the custom targets. Eventually,
                    # there will be a helper script to combine individual files into a properly formatted combined file
                    try:
                        sample[self.analysistype].baitfile = glob(os.path.join(sample[self.analysistype].targetpath,
                                                                               '*.fasta'))[0]
                    # If the fasta file is missing, raise a custom error
                    except IndexError:
                        # Combine any .tfa files in the directory into a combined targets .fasta file
                        tfafiles = glob(os.path.join(sample[self.analysistype].targetpath, '*.tfa'))
                        if tfafiles:
                            combinetargets(tfafiles, sample[self.analysistype].targetpath)
                        try:
                            self.baitfile = glob(os.path.join(sample[self.analysistype].targetpath, '*.fasta'))[0]
                        except IndexError as e:
                            # noinspection PyPropertyAccess
                            e.args = [
                                'Cannot find the combined fasta file in {}. Please note that the file must have a '
                                '.fasta extension'.format(sample[self.analysistype].targetpath)]
                            if os.path.isdir(sample[self.analysistype].targetpath):
                                raise
                            else:
                                sample[self.analysistype].runanalysis = False

                else:
                    sample[self.analysistype].runanalysis = False
            for sample in self.runmetadata:
                if sample.general.bestassemblyfile != 'NA' and sample[self.analysistype].runanalysis:
                    # Set the necessary attributes
                    sample[self.analysistype].outputdir = os.path.join(sample.run.outputdirectory, self.analysistype)
                    sample[self.analysistype].logout = os.path.join(sample[self.analysistype].outputdir, 'logout.txt')
                    sample[self.analysistype].logerr = os.path.join(sample[self.analysistype].outputdir, 'logerr.txt')
                    sample[self.analysistype].baitedfastq = \
                        os.path.join(sample[self.analysistype].outputdir,
                                     '{}_targetMatches.fastq.gz'.format(self.analysistype))
        else:
            # There is a relatively strict databasing scheme necessary for the custom targets. Eventually, there will
            # be a helper script to combine individual files into a properly formatted combined file
            try:
                self.baitfile = glob(os.path.join(self.targetpath, '*.fasta'))[0]
            # If the fasta file is missing, raise a custom error
            except IndexError:
                # Combine any .tfa files in the directory into a combined targets .fasta file
                tfafiles = glob(os.path.join(self.targetpath, '*.tfa'))
                if tfafiles:
                    combinetargets(tfafiles, self.targetpath)
                try:
                    self.baitfile = glob(os.path.join(self.targetpath, '*.fasta'))[0]
                except IndexError as e:
                    # noinspection PyPropertyAccess
                    e.args = ['Cannot find the combined fasta file in {}. Please note that the file must have a '
                              '.fasta extension'.format(self.targetpath)]
                    raise
            # Set all the necessary attributes
            for sample in self.runmetadata:
                setattr(sample, self.analysistype, GenObject())
                # Set attributes
                sample[self.analysistype].runanalysis = True
                sample[self.analysistype].baitfile = self.baitfile
                sample[self.analysistype].hashfile = self.hashfile
                sample[self.analysistype].hashcall = self.hashcall
                sample[self.analysistype].targetpath = self.targetpath
                sample[self.analysistype].outputdir = os.path.join(sample.run.outputdirectory, self.analysistype)
                sample[self.analysistype].logout = os.path.join(sample[self.analysistype].outputdir, 'logout.txt')
                sample[self.analysistype].logerr = os.path.join(sample[self.analysistype].outputdir, 'logerr.txt')
                sample[self.analysistype].baitedfastq = \
                    os.path.join(sample[self.analysistype].outputdir,
                                 '{}_targetMatches.fastq.gz'.format(self.analysistype))

    def bait(self, maskmiddle='f', k='19'):
        """
        Use bbduk to perform baiting
        :param maskmiddle: boolean argument treat the middle base of a kmer as a wildcard; increases sensitivity
        in the presence of errors.
        :param k: keyword argument for length of kmers to use in the analyses
        """
        logging.info('Performing kmer baiting of fastq files with {} targets'.format(self.analysistype))
        # There seems to be some sort of issue with java incorrectly calculating the total system memory on certain
        # computers. For now, calculate the memory, and feed it into the bbduk call
        with progressbar(self.runmetadata) as bar:
            for sample in bar:
                if sample.general.bestassemblyfile != 'NA' and sample[self.analysistype].runanalysis:
                    # Create the folder (if necessary)
                    make_path(sample[self.analysistype].outputdir)
                    # Make the system call
                    if len(sample.general.fastqfiles) == 2:
                        # Create the command to run the baiting - paired inputs and a single, zipped output
                        sample[self.analysistype].bbdukcmd = \
                            'bbduk.sh -Xmx{mem} ref={ref} in1={in1} in2={in2} k={kmer} maskmiddle={mm} ' \
                            'threads={c} outm={om}' \
                            .format(mem=self.mem,
                                    ref=sample[self.analysistype].baitfile,
                                    in1=sample.general.trimmedcorrectedfastqfiles[0],
                                    in2=sample.general.trimmedcorrectedfastqfiles[1],
                                    kmer=k,
                                    mm=maskmiddle,
                                    c=str(self.cpus),
                                    om=sample[self.analysistype].baitedfastq)
                    else:
                        sample[self.analysistype].bbdukcmd = \
                            'bbduk.sh -Xmx{mem} ref={ref} in={in1} k={kmer} maskmiddle={mm} ' \
                            'threads={cpus} outm={outm}' \
                            .format(mem=self.mem,
                                    ref=sample[self.analysistype].baitfile,
                                    in1=sample.general.trimmedcorrectedfastqfiles[0],
                                    kmer=k,
                                    mm=maskmiddle,
                                    cpus=str(self.cpus),
                                    outm=sample[self.analysistype].baitedfastq)
                    # Run the system call (if necessary)
                    if not os.path.isfile(sample[self.analysistype].baitedfastq):
                        out, err = run_subprocess(sample[self.analysistype].bbdukcmd)
                        write_to_logfile(sample[self.analysistype].bbdukcmd,
                                         sample[self.analysistype].bbdukcmd,
                                         self.logfile, sample.general.logout, sample.general.logerr,
                                         sample[self.analysistype].logout, sample[self.analysistype].logerr)
                        write_to_logfile(out,
                                         err,
                                         self.logfile, sample.general.logout, sample.general.logerr,
                                         sample[self.analysistype].logout, sample[self.analysistype].logerr)

    def reversebait(self, maskmiddle='f', k=19):
        """
        Use the freshly-baited FASTQ files to bait out sequence from the original target files. This will reduce the
        number of possibly targets against which the baited reads must be aligned
        """
        logging.info('Performing reverse kmer baiting of targets with FASTQ files')
        with progressbar(self.runmetadata) as bar:
            for sample in bar:
                if sample.general.bestassemblyfile != 'NA' and sample[self.analysistype].runanalysis:
                    outfile = os.path.join(sample[self.analysistype].outputdir, 'baitedtargets.fa')
                    sample[self.analysistype].revbbdukcmd = \
                        'bbduk.sh -Xmx{mem} ref={ref} in={in1} k={kmer} threads={cpus} mincovfraction={mcf} ' \
                        'maskmiddle={mm} outm={outm}' \
                        .format(mem=self.mem,
                                ref=sample[self.analysistype].baitedfastq,
                                in1=sample[self.analysistype].baitfile,
                                kmer=k,
                                cpus=str(self.cpus),
                                mcf=self.cutoff,
                                mm=maskmiddle,
                                outm=outfile)
                    # Run the system call (if necessary)
                    if not os.path.isfile(outfile):
                        out, err = run_subprocess(sample[self.analysistype].revbbdukcmd)
                        write_to_logfile(sample[self.analysistype].bbdukcmd,
                                         sample[self.analysistype].bbdukcmd,
                                         self.logfile, sample.general.logout, sample.general.logerr,
                                         sample[self.analysistype].logout, sample[self.analysistype].logerr)
                        write_to_logfile(out,
                                         err,
                                         self.logfile, sample.general.logout, sample.general.logerr,
                                         sample[self.analysistype].logout, sample[self.analysistype].logerr)
                    # Set the baitfile to use in the mapping steps as the newly created outfile
                    sample[self.analysistype].baitfile = outfile

    def subsample_reads(self):
        """
        Subsampling of reads to 20X coverage of rMLST genes (roughly).
        To be called after rMLST extraction and read trimming, in that order.
        """
        logging.info('Subsampling {} reads'.format(self.analysistype))
        with progressbar(self.runmetadata) as bar:
            for sample in bar:
                if sample.general.bestassemblyfile != 'NA':
                    # Create the name of the subsampled read file
                    sample[self.analysistype].subsampledreads = os.path.join(
                        sample[self.analysistype].outputdir,
                        '{}_targetMatches_subsampled.fastq.gz'.format(self.analysistype))
                    # Set the reformat.sh command. It will be run multiple times, overwrite previous iterations
                    # each time. Use samplebasestarget to provide an approximate number of bases to include in the
                    # subsampled reads e.g. for rMLST: 700000 (approx. 35000 bp total length of genes x 20X coverage)
                    sample[self.analysistype].subsamplecmd = \
                        'reformat.sh in={} out={} overwrite samplebasestarget=700000' \
                        .format(sample[self.analysistype].baitedfastq,
                                sample[self.analysistype].subsampledreads)
                    if not os.path.isfile(sample[self.analysistype].subsampledreads):
                        # Run the call
                        out, err = run_subprocess(sample[self.analysistype].subsamplecmd)
                        write_to_logfile(sample[self.analysistype].subsamplecmd,
                                         sample[self.analysistype].subsamplecmd,
                                         self.logfile, sample.general.logout, sample.general.logerr,
                                         sample[self.analysistype].logout, sample[self.analysistype].logerr)
                        write_to_logfile(out,
                                         err,
                                         self.logfile, sample.general.logout, sample.general.logerr,
                                         sample[self.analysistype].logout, sample[self.analysistype].logerr)
                    # Update the variable to store the baited reads
                    sample[self.analysistype].baitedfastq = sample[self.analysistype].subsampledreads

    def mapping(self):
        logging.info('Performing reference mapping')
        for i in range(len(self.runmetadata)):
            # Start threads
            threads = Thread(target=self.map, args=())
            # Set the daemon to True - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        for sample in self.runmetadata:
            if sample.general.bestassemblyfile != 'NA' and sample[self.analysistype].runanalysis:
                # Set the path/name for the sorted bam file to be created
                sample[self.analysistype].sortedbam = os.path.join(sample[self.analysistype].outputdir,
                                                                   '{}_sorted.bam'.format(self.analysistype))
                # Remove the file extension of the bait file for use in the indexing command
                sample[self.analysistype].baitfilenoext = sample[self.analysistype].baitfile.split('.fasta')[0]
                # Use bowtie2 wrapper to create index the target file
                bowtie2build = Bowtie2BuildCommandLine(reference=sample[self.analysistype].baitfile,
                                                       bt2=sample[self.analysistype].baitfilenoext,
                                                       **self.builddict)
                # Use samtools wrapper to set up the bam sorting command
                samsort = SamtoolsSortCommandline(input=sample[self.analysistype].sortedbam,
                                                  o=True,
                                                  out_prefix="-")
                # Determine the location of the SAM header editing script
                scriptlocation = sipprCommon.editsamheaders.__file__
                samtools = [
                    # When bowtie2 maps reads to all possible locations rather than choosing a 'best' placement, the
                    # SAM header for that read is set to 'secondary alignment', or 256. Please see:
                    # http://davetang.org/muse/2014/03/06/understanding-bam-flags/ The script below reads in
                    # the stdin and subtracts 256 from headers which include 256
                    'python3 {}'.format(scriptlocation),
                    # Use samtools wrapper to set up the samtools view
                    SamtoolsViewCommandline(b=True,
                                            S=True,
                                            h=True,
                                            input_file="-"),
                    samsort]
                # Add custom parameters to a dictionary to be used in the bowtie2 alignment wrapper
                indict = {'--very-sensitive-local': True,
                          '-U': sample[self.analysistype].baitedfastq,
                          '-a': True,
                          '--threads': self.threads,
                          '--local': True}
                # Create the bowtie2 reference mapping command
                bowtie2align = Bowtie2CommandLine(bt2=sample[self.analysistype].baitfilenoext,
                                                  threads=self.threads,
                                                  samtools=samtools,
                                                  **indict)
                # Create the command to faidx index the bait file
                sample[self.analysistype].faifile = sample[self.analysistype].baitfile + '.fai'
                samindex = SamtoolsFaidxCommandline(reference=sample[self.analysistype].baitfile)
                # Add the commands (as strings) to the metadata
                sample[self.analysistype].samindex = str(samindex)
                # Add the commands to the queue. Note that the commands would usually be set as attributes of
                # the sample but there was an issue with their serialization when printing out the metadata
                if not os.path.isfile(sample[self.analysistype].baitfilenoext + '.1' + self.bowtiebuildextension):
                    try:
                        stdoutbowtieindex, stderrbowtieindex = \
                            map(StringIO, bowtie2build(cwd=sample[self.analysistype].targetpath))
                        # Write any error to a log file
                        if stderrbowtieindex:
                            # Write the standard error to log, bowtie2 puts alignment summary here
                            with open(os.path.join(sample[self.analysistype].targetpath,
                                                   '{}_bowtie_index.log'.format(self.analysistype)), 'a+') as log:
                                log.writelines(logstr(bowtie2build, stderrbowtieindex.getvalue(),
                                                      stdoutbowtieindex.getvalue()))
                        # Close the stdout and stderr streams
                        stdoutbowtieindex.close()
                        stderrbowtieindex.close()
                    except ApplicationError:
                        pass
                self.mapqueue.put((sample, bowtie2build, bowtie2align, samindex))
        self.mapqueue.join()

    def map(self):
        while True:
            try:
                # Get the necessary values from the queue
                sample, bowtie2build, bowtie2align, samindex = self.mapqueue.get()
                # Use samtools faidx to index the bait file - this will be used in the sample parsing
                if not os.path.isfile(sample[self.analysistype].faifile):
                    stdoutindex, stderrindex = map(StringIO, samindex(cwd=sample[self.analysistype].targetpath))
                    # Write any error to a log file
                    if stderrindex:
                        # Write the standard error to log, bowtie2 puts alignment summary here
                        with open(os.path.join(sample[self.analysistype].targetpath,
                                               '{}_samtools_index.log'.format(self.analysistype)), 'a+') as log:
                            log.writelines(logstr(samindex, stderrindex.getvalue(), stdoutindex.getvalue()))
                    # Close the stdout and stderr streams
                    stdoutindex.close()
                    stderrindex.close()
                # Only run the functions if the sorted bam files and the indexed bait file do not exist
                if not os.path.isfile(sample[self.analysistype].sortedbam):
                    # Set stdout to a stringIO stream
                    stdout, stderr = map(StringIO, bowtie2align(cwd=sample[self.analysistype].outputdir))
                    if stderr:
                        # Write the standard error to log, bowtie2 puts alignment summary here
                        with open(os.path.join(sample[self.analysistype].outputdir,
                                               '{}_bowtie_samtools.log'.format(self.analysistype)), 'a+') as log:
                            log.writelines(logstr([bowtie2align], stderr.getvalue(), stdout.getvalue()))
                    stdout.close()
                    stderr.close()
            except ApplicationError:
                pass
            self.mapqueue.task_done()

    def indexing(self):
        logging.info('Indexing sorted BAM files')
        for i in range(len(self.runmetadata)):
            # Send the threads to
            threads = Thread(target=self.index, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        for sample in self.runmetadata:
            if sample.general.bestassemblyfile != 'NA' and sample[self.analysistype].runanalysis:
                bamindex = SamtoolsIndexCommandline(input=sample[self.analysistype].sortedbam)
                sample[self.analysistype].sortedbai = sample[self.analysistype].sortedbam + '.bai'
                sample[self.analysistype].bamindex = str(bamindex)
                self.indexqueue.put((sample, bamindex))
        self.indexqueue.join()

    def index(self):
        while True:
            try:
                sample, bamindex = self.indexqueue.get()
                # Only make the call if the .bai file doesn't already exist
                if not os.path.isfile(sample[self.analysistype].sortedbai):
                    # Use cStringIO streams to handle bowtie output
                    stdout, stderr = map(StringIO, bamindex(cwd=sample[self.analysistype].outputdir))
                    if stderr:
                        # Write the standard error to log
                        with open(os.path.join(sample[self.analysistype].outputdir,
                                               '{}_samtools_bam_index.log'.format(self.analysistype)), 'a+') as log:
                            log.writelines(logstr(bamindex, stderr.getvalue(), stdout.getvalue()))
                    stderr.close()
            except ApplicationError:
                pass
            self.indexqueue.task_done()

    def parsing(self):
        logging.info('Loading sorted BAM files')
        for i in range(len(self.runmetadata)):
            # Send the threads to
            threads = Thread(target=self.reduce, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        for sample in self.runmetadata:
            if sample.general.bestassemblyfile != 'NA' and sample[self.analysistype].runanalysis:
                # Get the fai file into a dictionary to be used in parsing results
                try:
                    with open(sample[self.analysistype].faifile, 'r') as faifile:
                        for line in faifile:
                            data = line.split('\t')
                            try:
                                sample[self.analysistype].faidict[data[0]] = int(data[1])
                            except (AttributeError, KeyError):
                                sample[self.analysistype].faidict = dict()
                                sample[self.analysistype].faidict[data[0]] = int(data[1])
                    self.parsequeue.put(sample)
                except FileNotFoundError:
                    pass
        self.parsequeue.join()
        self.parsebam()

    def reduce(self):
        """
        Use pysam to parse the sorted bam file to determine the sequence of the query, as well as any features present
        """
        while True:
            sample = self.parsequeue.get()
            # Initialise variables to store query and reference sequences, as well as any feature such as insertions
            # or locations with internal soft clipped reads
            sample[self.analysistype].sequence = dict()
            sample[self.analysistype].referencesequences = dict()
            sample[self.analysistype].features = dict()
            # Load the baitfile using SeqIO to get the reference sequences
            self.record_dict[sample.name] = SeqIO.to_dict(SeqIO.parse(sample[self.analysistype].baitfile, 'fasta'))
            # Create a pysam alignment file object of the sorted bam file
            bamfile = pysam.AlignmentFile(sample[self.analysistype].sortedbam, 'rb')
            # Iterate through each read in the bamfile
            for record in bamfile.fetch():
                # Initialise the starting position of the read to 0
                readpos = 0
                # Set the name and starting position of the match with the reference
                contig = record.reference_name
                refpos = record.reference_start
                # Iterate through the cigar tuples in the cigartuples attribute
                for cigartuple in record.cigartuples:
                    # Split the cigar tuple into the cigar feature (0: match, 1: insertion, 2: deletion, 3: skip,
                    # 4: soft clipping, 5: hard clipping, 6: padding) and the length of the cigar feature
                    cigartype, cigarlength = cigartuple
                    # Treat the different cigar features accordingly. For matches, add the query sequence at each
                    # reference position to the dictionary
                    if cigartype == 0:  # match
                        # Variable to store the final position reached in the read
                        final = 0
                        for i in range(readpos, readpos + cigarlength):
                            try:
                                sample[self.analysistype].sequence[contig][refpos].append(record.query_sequence[i])
                            except KeyError:
                                try:
                                    sample[self.analysistype].sequence[contig][refpos] = list()
                                    sample[self.analysistype].sequence[contig][refpos].append(record.query_sequence[i])
                                except KeyError:
                                    sample[self.analysistype].sequence[contig] = dict()
                                    sample[self.analysistype].sequence[contig][refpos] = list()
                                    sample[self.analysistype].sequence[contig][refpos].append(record.query_sequence[i])
                            # Increment the reference position for each bases in the query
                            refpos += 1
                            # Increment the final position - I used a variable here rather than incrementing readpos,
                            # as I didn't want any issues with incrementing readpos while in a loop reference readpos
                            final = i + 1
                        # Set the final position of the read
                        readpos = final
                    # Add the query sequence from the insertion, but don't increment the refpos; will increase the
                    # length of the query sequence compared to the reference
                    elif cigartype == 1:  # insertions
                        final = 0
                        for i in range(readpos, readpos + cigarlength):
                            try:
                                sample[self.analysistype].sequence[contig][refpos].append(record.query_sequence[i])
                            except KeyError:
                                try:
                                    sample[self.analysistype].sequence[contig][refpos] = list()
                                    sample[self.analysistype].sequence[contig][refpos].append(record.query_sequence[i])
                                except KeyError:
                                    sample[self.analysistype].sequence[contig] = dict()
                                    sample[self.analysistype].sequence[contig][refpos] = list()
                                    sample[self.analysistype].sequence[contig][refpos].append(record.query_sequence[i])
                            # Don't increment refpos, as this insertion is occurring between reference bases
                            final = i + 1
                        # Add the insertion feature to the dictionary
                        for i in range(readpos, readpos + cigarlength):
                            try:
                                sample[self.analysistype].features[contig][refpos].append('insertion')
                            except KeyError:
                                try:
                                    sample[self.analysistype].features[contig][refpos] = list()
                                    sample[self.analysistype].features[contig][refpos].append('insertion')
                                except KeyError:
                                    sample[self.analysistype].features[contig] = dict()
                                    sample[self.analysistype].features[contig][refpos] = list()
                                    sample[self.analysistype].features[contig][refpos].append('insertion')
                        # Set the final read position
                        readpos = final
                    # Add gaps (-) to the query sequence, but don't increment the readpos; will ensure that the
                    # reference and query sequences stay the same length
                    elif cigartype == 2:  # deletion
                        for i in range(readpos, readpos + cigarlength):
                            try:
                                sample[self.analysistype].sequence[contig][refpos].append('-')
                            except KeyError:
                                try:
                                    sample[self.analysistype].sequence[contig][refpos] = list()
                                    sample[self.analysistype].sequence[contig][refpos].append('-')
                                except KeyError:
                                    sample[self.analysistype].sequence[contig] = dict()
                                    sample[self.analysistype].sequence[contig][refpos] = list()
                                    sample[self.analysistype].sequence[contig][refpos].append('-')
                            # Don't increment read pos, as the deletion is occurring between query bases
                            refpos += 1
                    # Don't worry about skips yet. Have not found this cigar feature in any datasets so far
                    elif cigartype == 3:  # skip
                        pass
                    # An issue that was occurring with the parsing was internal soft clipping. Essentially, 5' reads
                    # would be soft right-clipped, and 3' reads would be soft left-clipped. At the resulting junction
                    # between the two clipped reads, the sequence data would look good, but only because of this
                    # undesired clipping. Add the internal soft clip feature to the dictionary
                    elif cigartype == 4:  # soft clipping
                        record_length = float(len(str(self.record_dict[sample.name][record.reference_name].seq)))
                        record_length_ninety = record_length * 0.95
                        # Determine if a soft clip is internal
                        if float(record.reference_start) >= (record_length - record_length_ninety) \
                                and float(record.reference_end) <= record_length_ninety:
                            try:
                                sample[self.analysistype].features[contig][refpos].append('internal soft clip')
                            except KeyError:
                                try:
                                    sample[self.analysistype].features[contig][refpos] = list()
                                    sample[self.analysistype].features[contig][refpos].append('internal soft clip')
                                except KeyError:
                                    sample[self.analysistype].features[contig] = dict()
                                    sample[self.analysistype].features[contig][refpos] = list()
                                    sample[self.analysistype].features[contig][refpos].append('internal soft clip')
                        # Increment the readpos by the length of the soft clipping feature
                        readpos += cigarlength
                    # Don't worry about hard clipping. Have not found this cigar feature in any datasets so far
                    elif cigartype == 5:  # hard clipping
                        pass
                    # Don't worry about padding yet. Have not found this cigar feature in any datasets so far
                    elif cigartype == 6:  # padding
                        pass
            self.parsequeue.task_done()

    def parsebam(self):
        """
        Parse the dictionaries of the sorted bam files extracted using pysam
        """
        logging.info('Parsing BAM')
        with progressbar(self.runmetadata) as bar:
            for sample in bar:
                # Initialise dictionaries to store parsed data
                matchdict = dict()
                depthdict = dict()
                seqdict = dict()
                snplocationsdict = dict()
                gaplocationsdict = dict()
                maxdict = dict()
                mindict = dict()
                deviationdict = dict()
                sample[self.analysistype].results = dict()
                sample[self.analysistype].avgdepth = dict()
                sample[self.analysistype].resultssnp = dict()
                sample[self.analysistype].snplocations = dict()
                sample[self.analysistype].resultsgap = dict()
                sample[self.analysistype].gaplocations = dict()
                sample[self.analysistype].sequences = dict()
                sample[self.analysistype].maxcoverage = dict()
                sample[self.analysistype].mincoverage = dict()
                sample[self.analysistype].standarddev = dict()
                # Iterate through each contig in the dictionary
                try:
                    for contig, poslist in sorted(sample[self.analysistype].sequence.items()):
                        # Use the record_dict dictionary with the contig as the key in order to pull out the
                        # reference sequence
                        refseq = self.record_dict[sample.name][contig].seq
                        # Initialise the reference position to 0
                        refpos = 0
                        # Initialise dictionaries with the contig name
                        matchdict[contig] = int()
                        depthdict[contig] = int()
                        seqdict[contig] = str()
                        snplocationsdict[contig] = list()
                        gaplocationsdict[contig] = list()
                        maxdict[contig] = int()
                        mindict[contig] = int()
                        deviationdict[contig] = list()
                        # Iterate through all the reference positions in the reference sequence
                        for pos, baselist in poslist.items():
                            # If the query position is equal to the reference position, proceed with the loop
                            if pos == refpos:
                                # Use the counter function to count the number of times each base appears in the list of
                                # query bases
                                counted = Counter(baselist)
                                # Sort the bases based on how common they are in the list
                                maxbases = counted.most_common()
                                # Set the query base as the most common - note that this is fairly simplistic - the base
                                # with the highest representation (or the first base in case of a tie) is used
                                querybase = maxbases[0][0]
                                # Extract the corresponding base in the reference sequence
                                refbase = refseq[pos]
                                # The depth of the current position is the length of the list of bases
                                depth = len(baselist)
                                # Update the sequence, depth, and deviation (list of depths) dictionaries
                                seqdict[contig] += querybase
                                depthdict[contig] += depth
                                deviationdict[contig].append(depth)
                                # Set the maximum and minimum depths observed
                                if depth > maxdict[contig]:
                                    maxdict[contig] = depth
                                # Initialise the minimum depth dictionary here
                                if contig not in mindict:
                                    mindict[contig] = depth
                                if depth < mindict[contig]:
                                    mindict[contig] = depth
                                # If the reference and query bases match, increment the number of matches
                                if querybase == refbase:
                                    matchdict[contig] += 1
                                # Using the NCBI 16S database, I observed that degenerate nucleotides were used. This
                                # allows for matches to occur to these bases
                                elif self.analysistype == 'sixteens_full' and refbase not in ['A', 'C', 'G', 'T']:
                                    # If the query base matches the corresponding IUPAC code e.g. A or G will match R,
                                    # increment the number of matches
                                    if querybase in self.iupac[refbase]:
                                        matchdict[contig] += 1
                                    # Otherwise treat the base as a mismatch, and add the base position to the list of
                                    # SNPs
                                    else:
                                        snplocationsdict[contig].append(pos)
                                # If the reference and query bases don't match there could be a couple of reasons
                                else:
                                    # If the bases simply don't match, add the position to the list of SNPs
                                    if querybase != '-':
                                        snplocationsdict[contig].append(pos)
                                    # However, if the query base is a gap (-), add the position to the dictionary of gap
                                    # locations
                                    else:
                                        gaplocationsdict[contig].append(pos)
                            # If the reference position does not equal the query position, assume a gap
                            else:
                                # Add the location of the gap to the dictionary
                                gaplocationsdict[contig].append(pos)
                            # Increment the reference position
                            refpos += 1
                except (AttributeError, KeyError):
                    pass
                # Iterate through all the results, and filter out sequences that do not meet the depth and/or
                # the sequence identity thresholds
                for allele in seqdict:
                    # If the length of the match is greater or equal to the length of the gene/allele (multiplied by the
                    # cutoff value) as determined using faidx indexing, then proceed
                    if matchdict[allele] >= sample[self.analysistype].faidict[allele] * self.cutoff:
                        # Calculate the average depth by dividing the total number of reads observed by the
                        # length of the gene
                        averagedepth = float(depthdict[allele]) / float(matchdict[allele])
                        percentidentity = \
                            float(matchdict[allele]) / float(sample[self.analysistype].faidict[allele]) * 100
                        # Only report a positive result if this average depth is greater than the desired average depth
                        # and if the percent identity is greater or equal to the cutoff
                        if averagedepth > self.averagedepth and percentidentity >= float(self.cutoff * 100):
                            # Populate resultsdict with the gene/allele name, the percent identity, and the avg depth
                            sample[self.analysistype].results.update({allele: '{:.2f}'.format(percentidentity)})
                            sample[self.analysistype].avgdepth.update({allele: '{:.2f}'.format(averagedepth)})
                            # Add the results to dictionaries
                            sample[self.analysistype].resultssnp.update({allele: len(snplocationsdict[allele])})
                            sample[self.analysistype].snplocations.update({allele: snplocationsdict[allele]})
                            sample[self.analysistype].resultsgap.update({allele: len(gaplocationsdict[allele])})
                            sample[self.analysistype].gaplocations.update({allele: gaplocationsdict[allele]})
                            sample[self.analysistype].sequences.update({allele: seqdict[allele]})
                            sample[self.analysistype].maxcoverage.update({allele: maxdict[allele]})
                            sample[self.analysistype].mincoverage.update({allele: mindict[allele]})
                            sample[self.analysistype] \
                                .standarddev.update({allele: '{:.2f}'.format(numpy.std(deviationdict[allele], ddof=1))})

    def clipper(self):
        """
        Filter out results based on the presence of cigar features such as internal soft-clipping
        """
        for sample in self.runmetadata:
            # Create a dictionary to store all the samples that do not have features
            replacementresults = dict()
            try:
                # SixteenS analyses seem to fail if results are filtered out
                if self.analysistype != 'sixteens_full' and self.analysistype != 'resfinder':
                    # Iterate through all the baited genes
                    for gene in sample[self.analysistype].faidict:
                        try:
                            percentidentity = sample[self.analysistype].results[gene]
                            try:
                                # Create a list to store whether a feature is present in enough reads to discard the
                                # sample
                                passingfeature = list()
                                for location, feature in sample[self.analysistype].features[gene].items():
                                    # If the feature is present in under 30% of the reads, set the passing variable
                                    # to true
                                    if len(feature) < int(float(sample[self.analysistype].avgdepth[gene])) * 0.3:
                                        passingfeature.append(True)
                                    # Otherwise set it to false
                                    else:
                                        passingfeature.append(False)
                                # If all the features are 'true' (present in fewer than 30% of the reads), add this
                                # contig to the list of passing results
                                if all(passingfeature):
                                    replacementresults[gene] = percentidentity
                            # If the allele does not have any features, it is added to the passing list
                            except KeyError:
                                replacementresults[gene] = percentidentity
                        except KeyError:
                            pass
                    # Update the .results attribute with the filtered dictionary
                    sample[self.analysistype].results = replacementresults
            except AttributeError:
                pass
        # Remove the features attribute - it takes up a lot of room in the .json file
        for sample in self.runmetadata:
            try:
                delattr(sample[self.analysistype], 'features')
            except AttributeError:
                pass

    # noinspection PyDefaultArgument
    def __init__(self, inputobject, cutoff=0.98, averagedepth=2):
        self.path = inputobject.path
        self.sequencepath = inputobject.sequencepath
        self.targetpath = inputobject.targetpath
        self.reportpath = inputobject.reportpath
        try:
            self.runmetadata = inputobject.runmetadata.samples
        except AttributeError:
            self.runmetadata = inputobject.runmetadata
        self.start = inputobject.starttime
        self.analysistype = inputobject.analysistype
        self.cpus = inputobject.cpus
        self.threads = inputobject.threads
        self.pipeline = inputobject.pipeline
        self.homepath = inputobject.homepath
        self.taxonomy = inputobject.taxonomy
        self.logfile = inputobject.logfile
        try:
            self.portallog = inputobject.portallog
        except AttributeError:
            self.portallog = ''
        self.cutoff = cutoff
        self.mem = int(0.85 * float(psutil.virtual_memory().total))
        self.builddict = dict()
        self.bowtiebuildextension = '.bt2'
        try:
            self.averagedepth = inputobject.averagedepth
        except AttributeError:
            self.averagedepth = averagedepth
        self.baitfile = str()
        self.hashfile = str()
        self.hashcall = str()
        self.devnull = open(os.devnull, 'wb')  # define /dev/null
        self.baitqueue = Queue(maxsize=self.cpus)
        self.mapqueue = Queue(maxsize=self.cpus)
        self.indexqueue = Queue(maxsize=self.cpus)
        self.parsequeue = Queue(maxsize=self.cpus)
        self.iupac = {
            'R': ['A', 'G'],
            'Y': ['C', 'T'],
            'S': ['G', 'C'],
            'W': ['A', 'T'],
            'K': ['G', 'T'],
            'M': ['A', 'C'],
            'B': ['C', 'G', 'T'],
            'D': ['A', 'G', 'T'],
            'H': ['A', 'C', 'T'],
            'V': ['A', 'C', 'G'],
            'N': ['A', 'C', 'G', 'T'],
            '-': ['-']
        }
        # Always perform reverse baiting - may want to change this later, so will keep this variable for now
        self.revbait = True
        self.record_dict = dict()
        # Run the analyses
        self.main()
        # Print the metadata
        printer = MetadataPrinter(self)
        printer.printmetadata()
