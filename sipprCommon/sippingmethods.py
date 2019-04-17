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
import multiprocessing
from collections import Counter
from click import progressbar
from threading import Thread
from io import StringIO
from queue import Queue
from glob import glob
import tempfile
import logging
import psutil
import numpy
import pysam
import json
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
        self.parsebam()
        # Filter out any sequences with cigar features such as internal soft-clipping from the results
        # self.clipper()

    def targets(self):
        """
        Search the targets folder for FASTA files, create the multi-FASTA file of all targets if necessary, and
        populate objects
        """
        logging.info('Performing analysis with {at} targets folder'.format(at=self.analysistype))
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
                                'Cannot find the combined fasta file in {path}. Please note that the file must have a '
                                '.fasta extension'.format(path=sample[self.analysistype].targetpath)]
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
                                     '{at}_targetMatches.fastq.gz'.format(at=self.analysistype))
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
                    e.args = ['Cannot find the combined fasta file in {path}. Please note that the file must have a '
                              '.fasta extension'.format(path=self.targetpath)]
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
                                 '{at}_targetMatches.fastq.gz'.format(at=self.analysistype))

    def bait(self, maskmiddle='f', k='19'):
        """
        Use bbduk to perform baiting
        :param maskmiddle: boolean argument treat the middle base of a kmer as a wildcard; increases sensitivity
        in the presence of errors.
        :param k: keyword argument for length of kmers to use in the analyses
        """
        logging.info('Performing kmer baiting of fastq files with {at} targets'.format(at=self.analysistype))
        # There seems to be some sort of issue with java incorrectly calculating the total system memory on certain
        # computers. For now, calculate the memory, and feed it into the bbduk call
        if self.kmer_size is None:
            kmer = k
        else:
            kmer = self.kmer_size
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
                                    kmer=kmer,
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
                                    kmer=kmer,
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
        if self.kmer_size is None:
            kmer = k
        else:
            kmer = self.kmer_size
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
                                kmer=kmer,
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
        logging.info('Subsampling {at} reads'.format(at=self.analysistype))
        with progressbar(self.runmetadata) as bar:
            for sample in bar:
                if sample.general.bestassemblyfile != 'NA':
                    # Create the name of the subsampled read file
                    sample[self.analysistype].subsampledreads = os.path.join(
                        sample[self.analysistype].outputdir,
                        '{at}_targetMatches_subsampled.fastq.gz'.format(at=self.analysistype))
                    # Set the reformat.sh command. It will be run multiple times, overwrite previous iterations
                    # each time. Use samplebasestarget to provide an approximate number of bases to include in the
                    # subsampled reads e.g. for rMLST: 700000 (approx. 35000 bp total length of genes x 20X coverage)
                    sample[self.analysistype].subsamplecmd = \
                        'reformat.sh in={bf} out={ssr} overwrite samplebasestarget=700000' \
                        .format(bf=sample[self.analysistype].baitedfastq,
                                ssr=sample[self.analysistype].subsampledreads)
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
        for i in range(self.cpus):
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
                                                                   '{at}_sorted.bam'.format(at=self.analysistype))
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
                    'python3 {sl}'.format(sl=scriptlocation),
                    # Use samtools wrapper to set up the samtools view
                    SamtoolsViewCommandline(b=True,
                                            S=True,
                                            h=True,
                                            F=4,
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
                                                   '{at}_bowtie_index.log'.format(at=self.analysistype)), 'a+') as log:
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
                                               '{at}_samtools_index.log'.format(at=self.analysistype)), 'a+') as log:
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
                                               '{at}_bowtie_samtools.log'.format(at=self.analysistype)), 'a+') as log:
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
                sample[self.analysistype].sortedbai = sample[self.analysistype].sortedbam + '.bai'
                self.indexqueue.put(sample)
        self.indexqueue.join()

    def index(self):
        while True:
            try:
                sample = self.indexqueue.get()
                bamindex = SamtoolsIndexCommandline(input=sample[self.analysistype].sortedbam)
                sample[self.analysistype].bamindex = str(bamindex)
                # Only make the call if the .bai file doesn't already exist
                if not os.path.isfile(sample[self.analysistype].sortedbai):
                    # Use cStringIO streams to handle bowtie output
                    stdout, stderr = map(StringIO, bamindex(cwd=sample[self.analysistype].outputdir))
                    if stderr:
                        # Write the standard error to log
                        with open(os.path.join(sample[self.analysistype].outputdir, '{at}_samtools_bam_index.log'
                                  .format(at=self.analysistype)), 'a+') as log:
                            log.writelines(logstr(bamindex, stderr.getvalue(), stdout.getvalue()))
                    stderr.close()
            except ApplicationError:
                pass
            self.indexqueue.task_done()

    @staticmethod
    def parse_one_sample(json_file, sample_name, best_assembly_file, analysistype, iupac, cutoff,
                         desired_average_depth, allow_soft_clips):
        with open(json_file) as f:
            sample = json.load(f)
        if 'faidict' not in sample:
            sample['faidict'] = dict()
        sample['name'] = sample_name
        if best_assembly_file != 'NA' and sample['runanalysis']:
            # Get the fai file into a dictionary to be used in parsing results
            try:
                with open(sample['faifile'], 'r') as faifile:
                    for line in faifile:
                        data = line.split('\t')
                        try:
                            sample['faidict'][data[0]] = int(data[1])
                        except (AttributeError, KeyError):
                            sample['faidict'] = dict()
                            sample['faidict'][data[0]] = int(data[1])
            except FileNotFoundError:
                pass
        sample['results'] = dict()
        sample['avgdepth'] = dict()
        sample['resultssnp'] = dict()
        sample['snplocations'] = dict()
        sample['resultsgap'] = dict()
        sample['gaplocations'] = dict()
        sample['sequences'] = dict()
        sample['maxcoverage'] = dict()
        sample['mincoverage'] = dict()
        sample['standarddev'] = dict()
        # Initialise dictionaries to store parsed data
        matchdict = dict()
        depthdict = dict()
        seqdict = dict()
        snplocationsdict = dict()
        gaplocationsdict = dict()
        maxdict = dict()
        mindict = dict()
        deviationdict = dict()
        has_clips_dict = dict()
        if 'baitfile' in sample and 'sortedbam' in sample:
            # Iterate through each contig in our target file.
            for contig in SeqIO.parse(sample['baitfile'], 'fasta'):
                # analysis since it probably isn't actually there.
                bamfile = pysam.AlignmentFile(sample['sortedbam'], 'rb')
                # Initialise dictionaries with the contig name
                matchdict[contig.id] = int()
                depthdict[contig.id] = int()
                seqdict[contig.id] = str()
                snplocationsdict[contig.id] = list()
                gaplocationsdict[contig.id] = list()
                maxdict[contig.id] = int()
                mindict[contig.id] = int()
                deviationdict[contig.id] = list()
                has_clips_dict[contig.id] = False
                # Settings used here are important for making output match up with bamfile visualised in tablet
                for column in bamfile.pileup(contig.id,
                                             stepper='samtools',
                                             ignore_orphans=False,
                                             min_base_quality=0,
                                             fastafile=pysam.FastaFile(sample['baitfile'])):
                    # Find all the attributes!

                    # Read depth - just get number of aligned reads at that position.
                    depth = column.get_num_aligned()
                    # Need to have at least some bases aligned for this to work at all.
                    if depth == 0:
                        seqdict[contig.id] += '-'
                        deviationdict[contig.id].append(depth)
                    try:  # This almost always works, except for when we have very high depth.
                        # Get list of bases for our column, marking ends and adding indels samtools style.
                        """
                        From http://www.htslib.org/doc/samtools.html
                        In the pileup format (without -u or -g), each line represents a genomic position, 
                        consisting of chromosome name, 1-based coordinate, reference base, the number of reads 
                        covering the site, read bases, base qualities and alignment mapping qualities. 
                        Information on match, mismatch, indel, strand, mapping quality and start and end of a read 
                        are all encoded at the read base column. At this column, a dot stands for a match to the 
                        reference base on the forward strand, a comma for a match on the reverse strand, a '>' or 
                        '<' for a reference skip, `ACGTN' for a mismatch on the forward strand and `acgtn' for a 
                        mismatch on the reverse strand. A pattern `\\+[0-9]+[ACGTNacgtn]+' indicates there is an 
                        insertion between this reference position and the next reference position. The length of
                        the insertion is given by the integer in the pattern, followed by the inserted sequence. 
                        Similarly, a pattern `-[0-9]+[ACGTNacgtn]+' represents a deletion from the reference. 
                        The deleted bases will be presented as `*' in the following lines. Also at the read 
                        base column, a symbol `^' marks the start of a read. The ASCII of the character following
                        `^' minus 33 gives the mapping quality. A symbol `$' marks the end of a read segment.
                        """
                        baselist = column.get_query_sequences(mark_ends=True, add_indels=True)
                        # Make sure everything in baselist is upper case - double check at some point that the
                        # lower case letters don't really mean anything.
                        baselist = [x.upper() for x in baselist]
                        # Use the counter function to count the number of times each base appears in the list of
                        # query bases. Mark ends will marks end/start of reads with $/^, and also takes care of clips
                        counted = Counter(baselist)
                        # Set the query base as the most common - note that this is fairly simplistic - the base
                        # with the highest representation (or the first base in case of a tie) is used
                        maxbases = counted.most_common()
                        querybase = maxbases[0][0]

                        # The $ and ^ characters represent end and start of reads, respectively, or just reads
                        # that have been clipped. If we find these as most common base anywhere not 10 bases from
                        # start of a target or 10 bases from end of a target, we likely have internal clipping.
                        # Set our has_clips_dict for this contig to True so that we know to ignore it.
                        if '$' in querybase or '^' in querybase:
                            if 10 <= column.reference_pos <= len(contig.seq) - 10:
                                has_clips_dict[contig.id] = True
                            # Samtools puts the mapping quality as an ascii char as well as the ^ for read starts
                            # So you end up with ^!A (or similar).
                            # Only take the base itself - last in the string for begin clip, first for end clip.
                            if '$' in querybase:
                                querybase = querybase[0]
                            elif '^' in querybase:
                                querybase = querybase[-1]
                        reference_base = contig.seq[column.reference_pos]
                        # Deletions. See the quoted text above. Briefly, the first reference position with of a deletion
                        # will look something like T-2TA for position 1237 in reference gene C.coliRM4661_23S_1 for
                        # sample 2018-LET-0007. This indicates that the query base at this position is 'T', and there is
                        # a two bp deletion of the reference bases 'T' and 'A'. Split of the extra information, and
                        # save the query base sequence
                        if len(querybase) > 1:
                            if '-' in querybase:
                                querybase = querybase.split('-')[0]
                            # Insertions. Similar to insertions, query base will be T+2GC for position 1191 in reference
                            # C.jejuniNCTC11168 for sample 2018-LET-0009
                            if '+' in querybase:
                                querybase = querybase.split('+')[0]
                        # Continuing the above example, positions 1238 and 1239 are not present in the sample, and are
                        # returned as '*'. Change this '*' to a '-'
                        if '*' in querybase:
                            querybase = '-'
                        # Populate the data dictionaries that we'll need later.
                        # matchdict keeps track of how many identities we have.
                        if reference_base == querybase:
                            matchdict[contig.id] += 1
                        # Using the NCBI 16S database, I observed that degenerate nucleotides were used. This
                        # allows for matches to occur to these bases
                        elif analysistype == 'sixteens_full' and reference_base not in ['A', 'C', 'G', 'T']:
                            # If the query base matches the corresponding IUPAC code e.g. A or G will match R,
                            # increment the number of matches
                            if querybase in iupac[reference_base]:
                                matchdict[contig.id] += 1
                            # Otherwise treat the base as a mismatch, and add the base position to the list of
                            # SNPs
                            else:
                                snplocationsdict[contig.id].append(column.reference_pos + 1)

                        # depthdict keeps track of how many bases total are aligned against the target gene.
                        depthdict[contig.id] += depth

                        # seqdict is our query sequence
                        seqdict[contig.id] += querybase
                        # snplocationsdict keeps track of where snps are in sequence. Append reference position
                        # if ref and query don't match. Add 1, since reference_pos is 0 based.
                        if reference_base != querybase and '-' not in querybase:
                            snplocationsdict[contig.id].append(column.reference_pos + 1)

                        # gaplocations, much the same as snplocations. Need to check that querybase shows as a gap
                        if '-' in querybase:
                            gaplocationsdict[contig.id].append(column.reference_pos + 1)

                        # Also need to keep track of min and max depth for the gene.
                        if depth > maxdict[contig.id]:
                            maxdict[contig.id] = depth

                        if depth < mindict[contig.id]:
                            mindict[contig.id] = depth

                        # Finally, deviationdict just has the depths at every position so we can calculate stdev later
                        deviationdict[contig.id].append(depth)
                    # High depth columns can break the get_query_sequences, since there's a hard-coded 10000 limit
                    # there - appears that read starts count as 3 each if mark_ends is on since then you get the
                    # the base, the mapping quality, and a char telling you it's the start of a read.
                    # To get around this, iterate over the pileupreads for this column manually.
                    # https://github.com/pysam-developers/pysam/issues/727
                    # TODO: Unfortunate amounts of code duplication here - get a function or something written.
                    except AssertionError:  # Very high depth makes us hit an AssertionError - do some more manual
                        # parsing then
                        baselist = list()
                        start_end_count = 0
                        for pileupread in column.pileups:
                            if pileupread.query_position is not None:
                                # Figure out what base is present, and if it's at start or end of read
                                baselist.append(pileupread.alignment.query_sequence[pileupread.query_position])
                                if pileupread.is_head == 1 or pileupread.is_tail == 1:
                                    start_end_count += 1
                        counted = Counter(baselist)
                        # Set the query base as the most common - note that this is fairly simplistic - the base
                        # with the highest representation (or the first base in case of a tie) is used
                        maxbases = counted.most_common()
                        querybase = maxbases[0][0]
                        reference_base = contig.seq[column.reference_pos]
                        if start_end_count >= 0.5 * depth:
                            has_clips_dict[contig.id] = True
                        # Deletions. See above for additional details
                        if len(querybase) > 1:
                            querybase = querybase.split('-')[0]
                        if '*' in querybase:
                            querybase = '-'
                        # Populate the data dictionaries that we'll need later.
                        # matchdict keeps track of how many identities we have.
                        if reference_base == querybase:
                            matchdict[contig.id] += 1
                        # Using the NCBI 16S database, I observed that degenerate nucleotides were used. This
                        # allows for matches to occur to these bases
                        elif analysistype == 'sixteens_full' and reference_base not in ['A', 'C', 'G', 'T']:
                            # If the query base matches the corresponding IUPAC code e.g. A or G will match R,
                            # increment the number of matches
                            if querybase in iupac[reference_base]:
                                matchdict[contig] += 1
                            # Otherwise treat the base as a mismatch, and add the base position to the list of
                            # SNPs
                            else:
                                snplocationsdict[contig].append(column.reference_pos + 1)

                        # depthdict keeps track of how many bases total are aligned against the target gene.
                        depthdict[contig.id] += depth

                        # seqdict is our query sequence
                        seqdict[contig.id] += querybase

                        # snplocationsdict keeps track of where snps are in sequence. Append reference position
                        # if ref and query don't match. Add 1, since reference_pos is 0 based.
                        if reference_base != querybase and querybase != '-':
                            snplocationsdict[contig.id].append(column.reference_pos + 1)

                        # gaplocations, much the same as snplocations. Need to check that querybase shows as a gap
                        if querybase == '-':
                            gaplocationsdict[contig.id].append(column.reference_pos + 1)

                        # Also need to keep track of min and max depth for the gene.
                        if depth > maxdict[contig.id]:
                            maxdict[contig.id] = depth

                        if depth < mindict[contig.id]:
                            mindict[contig.id] = depth

                        # Finally, deviationdict just has the depths at every position so we can calculate std dev later
                        deviationdict[contig.id].append(depth)

                bamfile.close()
            # Iterate through all the parsed alleles, and filter out ones that are either too short, or, if desired,
            # contain soft clipped sequences
            for allele in seqdict:
                # If the length of the match is greater or equal to the length of the gene/allele (multiplied by the
                # cutoff value) as determined using faidx indexing, then proceed
                if matchdict[allele] >= sample['faidict'][allele] * cutoff:
                    if has_clips_dict[allele] is False or allow_soft_clips:
                        # Calculate the average depth by dividing the total number of reads observed by the
                        # length of the gene
                        averagedepth = float(depthdict[allele]) / float(matchdict[allele])
                        percentidentity = \
                            float(matchdict[allele]) / float(sample['faidict'][allele]) * 100
                        # Only report a positive result if this average depth is greater than the desired average depth
                        # and if the percent identity is greater or equal to the cutoff
                        if averagedepth > desired_average_depth and percentidentity >= float(cutoff * 100):
                            # Populate resultsdict with the gene/allele name, the percent identity, and the avg depth
                            sample['results'][allele] = '{:.2f}'.format(percentidentity)
                            sample['avgdepth'][allele] = '{:.2f}'.format(averagedepth)
                            # Add the results to dictionaries
                            sample['resultssnp'][allele] = len(snplocationsdict[allele])
                            sample['snplocations'][allele] = snplocationsdict[allele]
                            sample['resultsgap'][allele] = len(gaplocationsdict[allele])
                            sample['gaplocations'][allele] = gaplocationsdict[allele]
                            sample['sequences'][allele] = seqdict[allele]
                            sample['maxcoverage'][allele] = maxdict[allele]
                            sample['mincoverage'][allele] = mindict[allele]
                            sample['standarddev'][allele] = '{:.2f}'.format(numpy.std(deviationdict[allele], ddof=1))
        return sample

    def parsebam(self):
        """
        Parse the dictionaries of the sorted bam files extracted using pysam
        """
        # Threading is actually the worst - need multiprocessing to make this work at all
        logging.info('Parsing BAM files')
        # The sample objects are too big to get pickled. To hack our way around this, try to dump the sample object to
        # json, and have the processing function turn the object into a dictionary.
        json_files = list()
        with tempfile.TemporaryDirectory() as tmpdir:
            best_assemblies = list()
            sample_names = list()
            for sample in self.runmetadata:
                json_name = os.path.join(tmpdir, '{sn}.json'.format(sn=sample.name))
                best_assemblies.append(sample.general.bestassemblyfile)
                sample_names.append(sample.name)
                with open(json_name, 'w') as f:
                    json.dump(sample[self.analysistype].dump(), f, sort_keys=True, indent=4)
                json_files.append(json_name)
            p = multiprocessing.Pool(processes=self.cpus)
            analysis_type_list = [self.analysistype] * len(self.runmetadata)
            iupac_list = [self.iupac] * len(self.runmetadata)
            cutoff_list = [self.cutoff] * len(self.runmetadata)
            depth_list = [self.averagedepth] * len(self.runmetadata)
            allow_soft_clip_list = [self.allow_soft_clips] * len(self.runmetadata)
            sample_results = p.starmap(Sippr.parse_one_sample,
                                       zip(json_files, sample_names, best_assemblies, analysis_type_list,
                                           iupac_list, cutoff_list, depth_list, allow_soft_clip_list))
            p.close()
            p.join()
        # Since we had to json-ize the sample objects, we now need to update the metadata for everything.
        for sample in self.runmetadata:
                sample[self.analysistype].faidict = dict()
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
                # Figure out which of the sample results to use.
                for sample_result in sample_results:
                    if sample_result['name'] == sample.name:
                        sample[self.analysistype].faidict = sample_result['faidict']
                        sample[self.analysistype].results = sample_result['results']
                        sample[self.analysistype].avgdepth = sample_result['avgdepth']
                        sample[self.analysistype].resultssnp = sample_result['resultssnp']
                        sample[self.analysistype].snplocations = sample_result['snplocations']
                        sample[self.analysistype].resultsgap = sample_result['resultsgap']
                        sample[self.analysistype].gaplocations = sample_result['gaplocations']
                        sample[self.analysistype].sequences = sample_result['sequences']
                        sample[self.analysistype].maxcoverage = sample_result['maxcoverage']
                        sample[self.analysistype].mincoverage = sample_result['mincoverage']
                        sample[self.analysistype].standarddev = sample_result['standarddev']
        logging.info('Done parsing BAM files')

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

    # noinspection PyDefaultArgument
    def __init__(self, inputobject, cutoff=0.98, averagedepth=2, k=None, allow_soft_clips=False):
        self.path = inputobject.path
        self.kmer_size = k
        self.sequencepath = inputobject.sequencepath
        self.targetpath = inputobject.targetpath
        self.reportpath = inputobject.reportpath
        try:
            self.runmetadata = inputobject.runmetadata.samples
        except AttributeError:
            self.runmetadata = inputobject.runmetadata
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
        self.allow_soft_clips = allow_soft_clips
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
