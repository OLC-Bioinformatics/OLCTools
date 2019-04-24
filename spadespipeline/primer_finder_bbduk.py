#!/usr/bin/env python3
from accessoryFunctions.accessoryFunctions import relative_symlink, filer, GenObject, make_path, MetadataObject, \
    run_subprocess, SetupLogging
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import SeqIO
from Bio import Seq
from argparse import ArgumentParser
from click import progressbar
from itertools import product
from threading import Thread
from csv import DictReader
from queue import Queue
import multiprocessing
from glob import glob
import itertools
import operator
import logging
import time
import os
import sys


__author__ = 'adamkoziol, duceppemo'


class PrimerFinder(object):

    def main(self):
        """
        Run the necessary methods
        """
        logging.info('Preparing metadata')
        # If this script is run as part of a pipeline, the metadata objects will already exist
        if not self.metadata:
            self.filer()
        else:
            self.objectprep()
        # Use the number of metadata objects to calculate the number of cores to use per sample in multi-threaded
        # methods with sequence calls to multi-threaded applications
        try:
            self.threads = int(self.cpus / len(self.metadata)) if self.cpus / len(
                self.metadata) > 1 else 1
        except (TypeError, ZeroDivisionError):
            self.threads = self.cpus
        logging.info('Reading and formatting primers')
        self.primers()
        logging.info('Baiting .fastq files against primers')
        self.bait()
        logging.info('Baiting .fastq files against previously baited .fastq files')
        self.doublebait()
        logging.info('Assembling contigs from double-baited .fastq files')
        self.assemble_amplicon_spades()
        logging.info('Creating BLAST database')
        self.make_blastdb()
        logging.info('Running BLAST analyses')
        self.blastnthreads()
        logging.info('Parsing BLAST results')
        self.parseblast()
        logging.info('Clearing amplicon files from previous iterations')
        self.ampliconclear()
        logging.info('Creating reports')
        self.reporter()

    def filer(self):
        """
        Find the sequence files, and create a metadata object for each sample
        """
        # Get a list of all the .fastq files in the sequence path
        fastq = sorted(glob(os.path.join(self.sequencepath, '*fastq*')))
        # Initialise a list to store .fasta formatted files
        fasta = list()
        # Get all the files in the sequence path
        for sequencefile in glob(os.path.join(self.sequencepath, '*')):
            # Split the extension from the file
            extension = os.path.splitext(sequencefile)[1]
            # If the extension is one of the acceptable FASTA extensions, add it to the list of FASTA files
            if extension in self.extensions:
                fasta.append(sequencefile)
        # Use filer from accessoryFunctions to find the base file name of each .fastq file/pair
        fastqnames = filer(fastq)
        # Get all the .fastq files into metadata objects
        for fastqname in sorted(fastqnames):
            # Create the metadata object
            sample = MetadataObject()
            sample.name = os.path.basename(fastqname)
            # Create an attribute for the current analysis
            setattr(sample, self.analysistype, GenObject())
            # Set the destination folder
            sample[self.analysistype].outputdir = os.path.join(self.path, fastqname)
            # Make the destination folder
            make_path(sample[self.analysistype].outputdir)
            # Get the fastq files specific to the fastqname
            specificfastq = glob(os.path.join(self.sequencepath, '{fq}*.fastq*'.format(fq=fastqname)))
            # Set the file type for the downstream analysis
            sample[self.analysistype].filetype = 'fastq'
            # Link the files to the output folder
            for fastq in specificfastq:
                relative_symlink(src_file=fastq,
                                 output_dir=sample[self.analysistype].outputdir)
            # Initialise the general and run categories
            sample.general = GenObject()
            # Populate the .fastqfiles category of :self.metadata
            sample.general.fastqfiles = [fastq for fastq in
                                         glob(os.path.join(sample[self.analysistype].outputdir, '{fq}*.fastq*'
                                                           .format(fq=fastqname))) if 'trimmed' not in fastq]
            # Populate certain attributes in order to be compatible with other local software
            sample.general.bestassemblyfile = sample.general.fastqfiles
            sample.general.trimmedcorrectedfastqfiles = sample.general.fastqfiles
            sample.general.outputdirectory = sample[self.analysistype].outputdir
            sample[self.analysistype].baitedfastq = os.path.join(
                sample[self.analysistype].outputdir,
                '{at}_targetMatches.fastq.gz'.format(at=self.analysistype))
            # Append the metadata to the list of samples
            self.metadata.append(sample)
        # Add any .fasta formatted files to the metadata object
        for fastafile in sorted(fasta):
            sample = MetadataObject()
            sample.name = os.path.splitext(os.path.split(fastafile)[-1])[0]
            setattr(sample, self.analysistype, GenObject())
            # Set the destination folder
            sample[self.analysistype].outputdir = os.path.join(self.path, 'detailed_reports', sample.name)
            # Make the destination folder
            make_path(sample[self.analysistype].outputdir)
            # Set the file type for the downstream analysis
            sample[self.analysistype].filetype = 'fasta'
            # Link the .fasta files to :self.path/:filename
            relative_symlink(src_file=fastafile,
                             output_dir=sample[self.analysistype].outputdir)
            # Initialise the general and run categories
            sample.general = GenObject()
            # Populate the .fastqfiles category of :self.metadata
            sample.general.fastqfiles = [fastq for fastq in glob(os.path.join(sample[self.analysistype].outputdir,
                                                                              '{sn}*{ext}*'
                                                                              .format(
                                                                                  sn=sample.name,
                                                                                  ext=os.path.splitext(fastafile)[1])))]
            sample.general.bestassemblyfile = sample.general.fastqfiles
            sample.general.trimmedcorrectedfastqfiles = sample.general.fastqfiles
            sample.general.outputdirectory = sample[self.analysistype].outputdir
            sample[self.analysistype].baitedfastq = os.path.join(
                sample[self.analysistype].outputdir,
                '{at}_targetMatches.fastq.gz'.format(at=self.analysistype))
            sample[self.analysistype].assemblyfile = fastafile
            # Append the metadata to the list of samples
            self.metadata.append(sample)

    def objectprep(self):
        """
        If the script is being run as part of a pipeline, create and populate the objects for the current analysis
        """
        for sample in self.metadata:
            setattr(sample, self.analysistype, GenObject())
            # Set the destination folder
            sample[self.analysistype].outputdir = os.path.join(self.path, self.analysistype)
            # Make the destination folder
            make_path(sample[self.analysistype].outputdir)
            sample[self.analysistype].baitedfastq = os.path.join(
                sample[self.analysistype].outputdir,
                '{at}_targetMatches.fastq.gz'.format(at=self.analysistype))
            # Set the file type for the downstream analysis
            sample[self.analysistype].filetype = self.filetype
            if self.filetype == 'fasta':
                sample[self.analysistype].assemblyfile = sample.general.bestassemblyfile

    def primers(self):
        """
        Read in the primer file, and create a properly formatted output file that takes any degenerate bases
        into account
        """
        with open(self.formattedprimers, 'w') as formatted:
            for record in SeqIO.parse(self.primerfile, 'fasta'):
                # from https://stackoverflow.com/a/27552377 - find any degenerate bases in the primer sequence, and
                # create all possibilities as a list
                degenerates = Seq.IUPAC.IUPACData.ambiguous_dna_values
                try:
                    primerlist = list(map("".join, product(*map(degenerates.get, str(record.seq)))))
                except TypeError:
                    print("Invalid Primer Sequence: {seq}".format(seq=str(record.seq)))
                    sys.exit()
                # As the record.id is being updated in the loop below, set the name of the primer here so that will
                # be able to be recalled when setting the new record.ids
                primername = record.id
                # Iterate through all the possible primers created from any degenerate bases
                for index, primer in enumerate(primerlist):
                    # Update the primer name with the position in the list to keep the name unique
                    record.id = primername + '_{index}'.format(index=index)
                    # Clear the description, as, otherwise, it will be added, and there will be duplicate information
                    record.description = ''
                    # Create a seqrecord from the primer sequence
                    record.seq = Seq.Seq(primer)
                    # Write the properly-formatted records to file
                    SeqIO.write(record, formatted, 'fasta')
                    # Populate a dictionary to store the length of the primers - will be used in determining whether
                    # BLAST hits are full-length
                    self.faidict[record.id] = len(str(record.seq))
                    # Ensure that the kmer length used in the initial baiting is no larger than the shorted primer
                    if len(str(record.seq)) < self.klength:
                        self.klength = len(str(record.seq))

    def bait(self):
        """
        Use bbduk to bait FASTQ reads from input files using the primer file as the target
        """
        with progressbar(self.metadata) as bar:
            for sample in bar:
                if sample.general.bestassemblyfile != 'NA':
                    # Only need to perform baiting on FASTQ files
                    if sample[self.analysistype].filetype == 'fastq':
                        # Make the system call - allow for single- or paired-end reads
                        if len(sample.general.fastqfiles) == 2:
                            # Create the command to run the baiting - ref: primer file, k: shortest primer length
                            # in1, in2: paired inputs, hdist: number of mismatches, interleaved: use interleaved output
                            # outm: single, zipped output file of reads that match the target file
                            sample[self.analysistype].bbdukcmd = \
                                'bbduk.sh ref={primerfile} k={klength} in1={forward} in2={reverse} ' \
                                'hdist={mismatches} threads={threads} interleaved=t outm={outfile}' \
                                .format(primerfile=self.formattedprimers,
                                        klength=self.klength,
                                        forward=sample.general.trimmedcorrectedfastqfiles[0],
                                        reverse=sample.general.trimmedcorrectedfastqfiles[1],
                                        mismatches=self.mismatches,
                                        threads=str(self.cpus),
                                        outfile=sample[self.analysistype].baitedfastq)
                        else:
                            sample[self.analysistype].bbdukcmd = \
                                'bbduk.sh ref={primerfile} k={klength} in={fastq} hdist={mismatches} ' \
                                'threads={threads} interleaved=t outm={outfile}' \
                                .format(primerfile=self.formattedprimers,
                                        klength=self.klength,
                                        fastq=sample.general.trimmedcorrectedfastqfiles[0],
                                        mismatches=self.mismatches,
                                        threads=str(self.cpus),
                                        outfile=sample[self.analysistype].baitedfastq)
                        # Run the system call (if necessary)
                        if not os.path.isfile(sample[self.analysistype].baitedfastq):
                            run_subprocess(sample[self.analysistype].bbdukcmd)

    def doublebait(self):
        """
        In order to ensure that there is enough sequence data to bridge the gap between the two primers, the paired
        .fastq files produced above will be used to bait the original input .fastq files
        """
        with progressbar(self.metadata) as bar:
            for sample in bar:
                if sample.general.bestassemblyfile != 'NA':
                    if sample[self.analysistype].filetype == 'fastq':
                        sample[self.analysistype].doublebaitedfastq = os.path.join(
                            sample[self.analysistype].outputdir,
                            '{at}_doubletargetMatches.fastq.gz'.format(at=self.analysistype))
                        # Make the system call
                        if len(sample.general.fastqfiles) == 2:
                            # Create the command to run the baiting
                            sample[self.analysistype].bbdukcmd2 = \
                                'bbduk.sh ref={baitfile} in1={forward} in2={reverse} hdist={mismatches} ' \
                                'threads={threads} interleaved=t outm={outfile}' \
                                .format(baitfile=sample[self.analysistype].baitedfastq,
                                        forward=sample.general.trimmedcorrectedfastqfiles[0],
                                        reverse=sample.general.trimmedcorrectedfastqfiles[1],
                                        mismatches=self.mismatches,
                                        threads=str(self.cpus),
                                        outfile=sample[self.analysistype].doublebaitedfastq)
                        else:
                            sample[self.analysistype].bbdukcmd2 = \
                                'bbduk.sh ref={baitfile} in={fastq} hdist={mismatches} threads={threads} ' \
                                'interleaved=t outm={outfile}' \
                                .format(baitfile=sample[self.analysistype].baitedfastq,
                                        fastq=sample.general.trimmedcorrectedfastqfiles[0],
                                        mismatches=self.mismatches,
                                        threads=str(self.cpus),
                                        outfile=sample[self.analysistype].doublebaitedfastq)
                        # Run the system call (if necessary)
                        if not os.path.isfile(sample[self.analysistype].doublebaitedfastq):
                            run_subprocess(sample[self.analysistype].bbdukcmd2)

    def assemble_amplicon_skesa(self):
        """
        Run skesa to assemble genomes
        """
        with progressbar(self.metadata) as bar:
            for sample in bar:
                # Initialise variables
                sample[self.analysistype].skesa_outdir = os.path.join(
                    sample[self.analysistype].outputdir, self.analysistype)
                make_path(sample[self.analysistype].skesa_outdir)
                sample[self.analysistype].assemblyfile = os.path.join(sample[self.analysistype].skesa_outdir,
                                                                      'contigs.fasta')
                # If there are two fastq files
                if len(sample.general.fastqfiles) == 2:
                    # Set the reverse fastq name https://github.com/ncbi/SKESA/issues/7
                    sample[self.analysistype].spadescommand = 'skesa --fastq {fastqfiles} --cores {threads} ' \
                                                              '--use_paired_ends --vector_percent 1 ' \
                                                              '--contigs_out {contigs}'\
                        .format(fastqfiles=','.join(sample.general.fastqfiles),
                                threads=self.cpus,
                                contigs=sample[self.analysistype].assemblyfile)
                # Same as above, but use single read settings for the assembler
                else:
                    sample[self.analysistype].spadescommand = 'skesa --fastq {fastqfiles} --cores {threads} ' \
                                                              '--vector_percent 1 --contigs_out {contigs}'\
                        .format(fastqfiles=','.join(sample.general.fastqfiles),
                                threads=self.cpus,
                                contigs=sample[self.analysistype].assemblyfile)
                if not os.path.isfile(sample[self.analysistype].assemblyfile):
                    # Create and run the sequence call as required
                    run_subprocess(sample[self.analysistype].spadescommand)
                    if not os.path.isfile(sample[self.analysistype].assemblyfile):
                        sample[self.analysistype].assemblyfile = 'NA'

    def assemble_amplicon_spades(self):
        """
        Use SPAdes to assemble the amplicons using the double-baited .fastq files
        """
        for _ in self.metadata:
            # Send the threads to the merge method. :args is empty as I'm using
            threads = Thread(target=self.assemble, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        for sample in self.metadata:
            if sample.general.bestassemblyfile != 'NA':
                if sample[self.analysistype].filetype == 'fastq':
                    sample[self.analysistype].spadesoutput = os.path.join(
                        sample[self.analysistype].outputdir, self.analysistype)
                    # Removed --careful, as there was an issue with the .fastq reads following baiting - something to
                    # do with the names, or the interleaving. Subsequent testing showed no real changes to assemblies
                    if len(sample.general.fastqfiles) == 2:
                        sample[self.analysistype].spadescommand = \
                            'spades.py -k {klength} --only-assembler --12 {fastq} -o {outfile} -t {threads}'\
                            .format(klength=self.kmers,
                                    fastq=sample[self.analysistype].doublebaitedfastq,
                                    outfile=sample[self.analysistype].spadesoutput,
                                    threads=self.threads)
                    else:
                        sample[self.analysistype].spadescommand = \
                            'spades.py -k {klength} --only-assembler -s {fastq} -o {outfile} -t {threads}' \
                            .format(klength=self.kmers,
                                    fastq=sample[self.analysistype].doublebaitedfastq,
                                    outfile=sample[self.analysistype].spadesoutput,
                                    threads=self.threads)
                    sample[self.analysistype].assemblyfile = os.path.join(sample[self.analysistype].spadesoutput,
                                                                          'contigs.fasta')
                    self.queue.put(sample)
        self.queue.join()

    def assemble(self):
        while True:
            sample = self.queue.get()
            if not os.path.isfile(sample[self.analysistype].assemblyfile):
                # Create and run the sequence call as required
                run_subprocess(sample[self.analysistype].spadescommand)
                if not os.path.isfile(sample[self.analysistype].assemblyfile):
                    sample[self.analysistype].assemblyfile = 'NA'
            self.queue.task_done()

    def make_blastdb(self):
        """
        Create a BLAST database of the primer file
        """
        # remove the path and the file extension for easier future globbing
        db = os.path.splitext(self.formattedprimers)[0]
        nhr = '{db}.nhr'.format(db=db)  # add nhr for searching
        if not os.path.isfile(str(nhr)):
            # Create the databases
            command = 'makeblastdb -in {primerfile} -parse_seqids -max_file_sz 2GB -dbtype nucl -out {outfile}'\
                .format(primerfile=self.formattedprimers,
                        outfile=db)
            run_subprocess(command)

    def blastnthreads(self):
        """
        Setup and create  threads for blastn and xml path
        """
        # Create the threads for the BLAST analysis
        for i in range(self.cpus):
            threads = Thread(target=self.runblast, args=())
            threads.setDaemon(True)
            threads.start()
        # Populate threads for each gene, genome combination
        for sample in self.metadata:
            if sample.general.bestassemblyfile != 'NA' and sample[self.analysistype].assemblyfile != 'NA':
                # Add each fasta file combination to the threads
                self.blastqueue.put(sample)
        # Join the threads
        self.blastqueue.join()

    def runblast(self):
        while True:
            sample = self.blastqueue.get()
            # Run the BioPython BLASTn module with the genome as query, fasta(target gene) as db.
            # Do not re-perform the BLAST search each time
            try:
                sample[self.analysistype].report = glob(os.path.join(
                    sample[self.analysistype].outputdir,
                    '{sn}*rawresults*'.format(sn=sample.name)))[0]
                size = os.path.getsize(sample[self.analysistype].report)
                # If a report was created, but no results entered - program crashed, or no sequences passed thresholds,
                # remove the report, and run the blast analyses again
                if size == 0:
                    os.remove(sample[self.analysistype].report)
                    sample[self.analysistype].report = os.path.join(sample[self.analysistype].outputdir,
                                                                    '{sn}_rawresults.csv'.format(sn=sample.name))
            except IndexError:
                sample[self.analysistype].report = os.path.join(
                    sample[self.analysistype].outputdir,
                    '{sn}_rawresults.csv'.format(sn=sample.name))
            db = os.path.splitext(self.formattedprimers)[0]
            # BLAST command line call. Note the high number of alignments.
            # Due to the fact that all the targets are combined into one database, this is to ensure that all potential
            # alignments are reported. Also note the custom outfmt: the doubled quotes are necessary to get it work
            blastn = NcbiblastnCommandline(query=sample[self.analysistype].assemblyfile,
                                           db=db,
                                           evalue=1e-1,
                                           task='blastn-short',
                                           num_alignments=1000000,
                                           num_threads=self.threads,
                                           outfmt="'6 qseqid sseqid positive mismatch gaps "
                                                  "evalue bitscore slen length qstart qend qseq sstart send sseq'",
                                           out=sample[self.analysistype].report)
            # Save the blast command in the metadata
            sample[self.analysistype].blastcommand = str(blastn)
            # Only run blast if the report doesn't exist
            # if not os.path.isfile(sample[self.analysistype].report):  # removed because don't want to
            # delete the blast files between each run when tweaking the settings and or primers
            try:
                blastn()
                # Add header to report
                blast_header = ['qseqid', 'sseqid', 'positive', 'mismatch', 'gaps', 'evalue', 'bitscore', 'slen',
                                'length', 'qstart', 'qend', 'qseq', 'sstart', 'send', 'sseq']
                with open(sample[self.analysistype].report, 'r+') as f:
                    content = f.read()
                    f.seek(0, 0)
                    f.write('\t'.join(blast_header) + '\n' + content)
            except:
                self.blastqueue.task_done()
                self.blastqueue.join()
                try:
                    os.remove(sample[self.analysistype].report)
                except IOError:
                    pass
                raise
            self.blastqueue.task_done()

    def parseblast(self):
        """
        Parse the BLAST results produced above. Find primer pairs with full-length hits with mismatches equal or
        lesser than the cutoff value
        """
        for sample in self.metadata:
            if sample.general.bestassemblyfile != 'NA' and sample[self.analysistype].assemblyfile != 'NA':
                # Initialise variables
                sample[self.analysistype].blastresults = dict()
                sample[self.analysistype].contigs = dict()
                sample[self.analysistype].hits = dict()
                sample[self.analysistype].mismatches = dict()
                sample[self.analysistype].blastrecords = list()
                sample[self.analysistype].range = dict()
                sample[self.analysistype].genespresent = dict()
                # Open blast output csv file
                csvfile = open(sample[self.analysistype].report)
                # Skip header
                csvfile.readline()
                # Open the sequence profile file as a dictionary
                blastdict = DictReader(csvfile, fieldnames=self.fieldnames, dialect='excel-tab')

                # Go through each BLAST result
                for row in blastdict:
                    # Ensure that the hit is full-length, and that the number of mismatches is equal to or lesser
                    # than the supplied cutoff value
                    if int(row['alignment_length']) == self.faidict[row['subject_id']] and \
                                    int(row['mismatches']) <= self.mismatches:
                        # Add the current row to the list for future work
                        sample[self.analysistype].blastrecords.append(row)
                        # Populate the dictionaries with the contig name (e.g. CA_CFIA-515_NODE_1_length_1791),
                        # the gene name (e.g. vtx2a), and the primer name (e.g. vtx2a-R3_1) as required
                        # accounts for primer names with "-" in addition to the terminal "-F" or "-R"
                        try:
                            sample[self.analysistype].blastresults[row['query_id']].add(row['subject_id'])
                            sample[self.analysistype].contigs[row['query_id']].add('-'.join(row['subject_id']
                                                                                            .split('-')[:-1]))
                        except KeyError:
                            sample[self.analysistype].blastresults[row['query_id']] = set()
                            sample[self.analysistype].blastresults[row['query_id']].add(row['subject_id'])
                            sample[self.analysistype].contigs[row['query_id']] = set()
                            sample[self.analysistype].contigs[row['query_id']].add('-'.join(row['subject_id']
                                                                                            .split('-')[:-1]))
                # Check to see if both forward and reverse primers are present for a particular gene within a contig
                for contig, genes in sample[self.analysistype].contigs.items():
                    # Split off the primer details (e.g. vtx2a-R3_1 -> vtx2a-R) from the blast results dictionary in
                    # order to create a searchable list of primers
                    # accounts for primer names with "-" in addition to the terminal "-F" or "-R"
                    reformatted = {'-'.join(['-'.join(x.split('-')[:-1]), x.split('-')[-1][0]])
                                   for x in sample[self.analysistype].blastresults[contig]}
                    # Iterate through the list of genes to check if primers are present
                    for gene in genes:
                        # Add -F and -R to the gene, and ensure that both options are in the reformatted list of genes
                        # e.g. vtx2a-F and vtx2a-R in [vtx1a-R, vtx2c-F ,vtx2a-F, vtx2a-R]
                        if gene + '-F' in reformatted and gene + '-R' in reformatted:
                            # Extract the precise primers from the dictionary e.g. vtx2a use to
                            # find vtx2a-F2_4 (forward) and vtx2a-R3_1 (reverse)
                            forwardprimers = [primer for primer in sample[self.analysistype].blastresults[contig]
                                              if gene == primer.split('-F')[0]]
                            reverseprimers = [primer for primer in sample[self.analysistype].blastresults[contig]
                                              if gene == primer.split('-R')[0]]
                            # Concatenate the lists
                            primers = forwardprimers + reverseprimers
                            # Populate the dictionary with the primers
                            try:
                                sample[self.analysistype].hits[contig].append(primers)
                            except KeyError:
                                sample[self.analysistype].hits[contig] = list()
                                sample[self.analysistype].hits[contig].append(primers)

                            for record in sample[self.analysistype].blastrecords:
                                for primer in primers:
                                    if record['query_id'] == contig and record['subject_id'] == primer:
                                        # Populate the dictionary with the primers
                                        try:
                                            sample[self.analysistype].mismatches[contig][gene]\
                                                .update({primer: int(record['mismatches'])})
                                        except KeyError:
                                            try:
                                                sample[self.analysistype].mismatches[contig][gene] = dict()
                                                sample[self.analysistype].mismatches[contig][gene] = \
                                                    {primer: int(record['mismatches'])}
                                            except KeyError:
                                                sample[self.analysistype].mismatches[contig] = dict()
                                                sample[self.analysistype].mismatches[contig][gene] = dict()
                                                sample[self.analysistype].mismatches[contig][gene] = \
                                                    {primer: int(record['mismatches'])}
                # Use query the stored blast dictionary to find the location of the hits
                for row in sample[self.analysistype].blastrecords:
                    try:
                        # Extract the primers corresponding to the contig
                        for primers in sample[self.analysistype].hits[row['query_id']]:
                            # Extract the name of the contig
                            contig = row['query_id']
                            # Iterate through the forward and reverse primers
                            for primer in primers:
                                # If the primer is present in the current row, then this is the row of interest
                                if row['subject_id'] == primer:
                                    # Split off the primer direction and numbering
                                    # accounts for primer names with "-" in addition to the terminal "-F" or "-R"
                                    gene = '-'.join(primer.split('-')[:-1])
                                    # Populate a dictionary for storing the genes present - will be used in creating
                                    # the report
                                    try:
                                        sample[self.analysistype].genespresent[contig].add(gene)
                                    except KeyError:
                                        sample[self.analysistype].genespresent[contig] = set()
                                        sample[self.analysistype].genespresent[contig].add(gene)
                                    # Populate the range of the hit - the forward primer will have a -F an the name
                                    if '-F' in primer:
                                        # Determine if the sequence is the reverse complement - based on the fact that
                                        # this is the forward primer, if the contig is reversed, then the primer
                                        # (subject) will be reversed.
                                        if int(row['subject_start']) > int(row['subject_end']):
                                            # For reversed sequences, take the larger value of the start and stop
                                            data = max(int(row['query_start']), int(row['query_end']))
                                        else:
                                            # Otherwise take the smaller value
                                            data = min(int(row['query_start']), int(row['query_end']))
                                        # Add the appropriately calculated value to the range dictionary
                                        try:
                                            sample[self.analysistype].range[contig][gene].add(data)
                                        except KeyError:
                                            try:
                                                sample[self.analysistype].range[contig][gene] = set()
                                                sample[self.analysistype].range[contig][gene].add(data)
                                            except KeyError:
                                                sample[self.analysistype].range[contig] = dict()
                                                sample[self.analysistype].range[contig][gene] = set()
                                                sample[self.analysistype].range[contig][gene].add(data)

                                    # Similar to the forward primer, except reverse the min() and max()
                                    elif '-R' in primer:
                                        if int(row['subject_start']) < int(row['subject_end']):
                                            data = min(int(row['query_start']), int(row['query_end']))
                                        else:
                                            data = max(int(row['query_start']), int(row['query_end']))
                                        # Add the appropriately calculated value to the range dictionary
                                        try:
                                            sample[self.analysistype].range[contig][gene].add(data)
                                        except KeyError:
                                            try:
                                                sample[self.analysistype].range[contig][gene] = set()
                                                sample[self.analysistype].range[contig][gene].add(data)
                                            except KeyError:
                                                sample[self.analysistype].range[contig] = dict()
                                                sample[self.analysistype].range[contig][gene] = set()
                                                sample[self.analysistype].range[contig][gene].add(data)
                    except KeyError:
                        pass

    def ampliconclear(self):
        """
        Clear previously created amplicon files to prepare for the appending of data to fresh files
        """
        for sample in self.metadata:
            # Set the name of the amplicon FASTA file
            sample[self.analysistype].ampliconfile = os.path.join(
                sample[self.analysistype].outputdir, '{sn}_amplicons.fa'.format(sn=sample.name))
            try:
                os.remove(sample[self.analysistype].ampliconfile)
            except IOError:
                pass

    def reporter(self):
        """
        Create reports of the analyses
        """
        # Create the report path if necessary
        make_path(self.reportpath)
        with open(self.report, 'w') as report:
            # Initialise the header
            data = 'Sample,Gene,GenomeLocation,AmpliconSize,Contig,ForwardPrimers,ReversePrimers,' \
                   'ForwardMismatches,ReverseMismatches\n'
            for sample in self.metadata:
                # Initialise variables to convert attributes from set to list
                sample[self.analysistype].genes = dict()
                sample[self.analysistype].ntrange = dict()
                try:
                    # Check to ensure that there were amplicons present
                    if sample[self.analysistype].range:
                        # Iterate through each contig with genes present, and the list of genes on the contig
                        for contig, genes in sorted(sample[self.analysistype].genespresent.items()):
                            sample[self.analysistype].genes[contig] = list(genes)
                            # Iterate through the set of genes
                            for gene in sorted(genes):
                                # Determine which primers had the best lowest number of mismatches
                                # Create variables to store primer names, and number of mismatches
                                forward = list()
                                reverse = list()
                                forwardmismatches = int()
                                reversemismatches = int()
                                # Sort the dictionary of mismatches based on the number of mismatches
                                sorteddict = sorted(sample[self.analysistype].mismatches[contig][gene].items(),
                                                    key=operator.itemgetter(1))
                                # For every primer in the sorted dictionary, check to see if it is a forward or
                                # reverse primer, and determine the number of mismatches; if it is the first primer
                                # encountered for a gene, add it to the list, as because the dictionary is sorted,
                                # the number of mismatches should be the best. Additionally, if subsequent primer
                                # hits have the same number of mismatches, also add them to the list
                                for entry in sorteddict:
                                    # Add primer to the list if it is the first primer encountered
                                    if '-F' in entry[0] and not forward:
                                        forward.append(entry[0])
                                        forwardmismatches = entry[1]
                                    # Add the primer to the list if it has the same number of mismatches as a previous
                                    # primer in the list
                                    elif '-F' in entry[0] and entry[1] == forwardmismatches:
                                        forward.append(entry[0])
                                        forwardmismatches = entry[1]
                                    elif '-R' in entry[0] and not reverse:
                                        reverse.append(entry[0])
                                        reversemismatches = entry[1]
                                    elif '-R' in entry[0] and entry[1] == reversemismatches:
                                        reverse.append(entry[0])
                                        reversemismatches = entry[1]

                                # Make a variable to prevent writing out this long attribute name multiple times
                                ntrange = list(sample[self.analysistype].range[contig][gene])
                                sample[self.analysistype].ntrange[gene] = ntrange
                                # Extract the amplicons from the sequence file
                                self.ampliconfile(sample, contig, sorted(ntrange), forward, reverse)
                                # This first gene for a sample gets the sample name printed
                                data += '{sn},'.format(sn=sample.name)
                                # Populate the string with the gene name, properly formatted range, the length of
                                # the amplicon, and the name of the contig on which the gene was found
                                data += '{},{},{},{},{},{},{},{}\n'\
                                    .format(gene,
                                            '-'.join(str(x) for x in sorted(ntrange)),
                                            max(ntrange) - min(ntrange) + 1,
                                            contig,
                                            ';'.join(sorted(forward)),
                                            ';'.join(sorted(reverse)),
                                            forwardmismatches,
                                            reversemismatches)
                    # If there were no amplicons, add the sample name and nothing else
                    else:
                        data += '{sn}\n'.format(sn=sample.name)
                # If there were no BLAST hits, add the sample name, and nothing else
                except AttributeError:
                    data += '{sn}\n'.format(sn=sample.name)
                # Remove attributes that either take up too much room in the .json output, or are not JSON serializable
                try:
                    delattr(sample[self.analysistype], "blastresults")
                except AttributeError:
                    pass
                try:
                    delattr(sample[self.analysistype], "genespresent")
                except AttributeError:
                    pass
                try:
                    delattr(sample[self.analysistype], "contigs")
                except AttributeError:
                    pass
                try:
                    delattr(sample[self.analysistype], "range")
                except AttributeError:
                    pass
            # Write the string to the report
            report.write(data)
        # Clean up the BLAST database files
        db = os.path.splitext(self.formattedprimers)[0]
        # A list of all the file extensions associated with the BLASTdb
        dbextensions = ['.nhr', '.nin', '.nog', '.nsd', '.nsi', '.nsq']
        # Iterate through all the files, and delete each one - pass on IO errors
        for dbfile in zip(itertools.repeat(db), dbextensions):
            try:
                os.remove(''.join(dbfile))
            except IOError:
                pass
        try:
            os.remove(self.formattedprimers)  # Maybe want to keep this file?
        except IOError:
            pass

    def ampliconfile(self, sample, contig, amplicon_range, forward_primer, reverse_primer):
        """
        Extracts amplicon sequence from contig file
        :param sample: sample metadata object
        :param contig: name of the contig hit by primers
        :param amplicon_range: range of the amplicon within the contig
        :param forward_primer: name of the forward primer
        :param reverse_primer: name of the reverse primer
        """
        # Open the file
        with open(sample[self.analysistype].ampliconfile, 'a') as ampliconfile:
            try:
                # Load the records from the assembly into the dictionary
                for record in SeqIO.parse(sample[self.analysistype].assemblyfile, 'fasta'):
                    if record.id == contig:
                        try:
                            # Extract the name of the gene from the primer name
                            genename = forward_primer[0].split('-')[0]
                            try:
                                # Sort the range calculated above
                                start = amplicon_range[0]
                                end = amplicon_range[1]
                                # Slice the gene sequence from the sequence record - remember to subtract one to
                                # allow for zero-based indexing
                                genesequence = str(record.seq)[int(start) - 1:int(end)]
                                # Set the record.id to be the sample name, the contig name,
                                # the range, and the primers
                                record.id = '{sn}_{contig}_{range}_{primers}' \
                                    .format(sn=sample.name,
                                            contig=contig,
                                            range='_'.join(str(x) for x in sorted(sample[self.analysistype]
                                                                                  .range[record.id][genename])),
                                            primers='_'.join(['_'.join(forward_primer), '_'.join(reverse_primer)]))
                                # Clear the record.description
                                record.description = ''
                                # Create a seq record from the sliced genome sequence
                                record.seq = Seq.Seq(genesequence)
                                # Write the amplicon to file
                                SeqIO.write(record, ampliconfile, 'fasta')
                            except IndexError:
                                pass
                        except AttributeError:
                            pass
            except FileNotFoundError:
                pass

    def __init__(self, args, analysistype, filetype='fastq'):
        if args.path.startswith('~'):
            self.path = os.path.abspath(os.path.expanduser(os.path.join(args.path)))
        else:
            self.path = os.path.abspath(os.path.join(args.path))
        if args.sequencepath.startswith('~'):
            self.sequencepath = os.path.abspath(os.path.expanduser(os.path.join(args.sequencepath)))
        else:
            self.sequencepath = os.path.abspath(os.path.join(args.sequencepath))
        self.start = args.start
        if args.primerfile.startswith('~'):
            self.primerfile = os.path.abspath(os.path.expanduser(os.path.join(args.primerfile)))
        else:
            self.primerfile = os.path.abspath(os.path.join(args.primerfile))
        self.mismatches = int(args.mismatches)
        try:
            self.metadata = args.runmetadata
        except AttributeError:
            self.metadata = list()
        # Use the argument for the number of threads to use, or default to the number of cpus in the system
        try:
            self.cpus = int(args.cpus)
        except (AttributeError, TypeError):
            self.cpus = multiprocessing.cpu_count()
        self.threads = int()
        try:
            self.analysistype = args.analysistype
        except AttributeError:
            self.analysistype = analysistype
        self.formattedprimers = os.path.join(self.path, 'formattedprimers.fa')
        self.faidict = dict()
        self.filetype = filetype
        # Use a long kmer for SPAdes assembly
        try:
            self.kmers = args.kmerlength
        except AttributeError:
            self.kmers = '99'
        self.queue = Queue()
        self.blastqueue = Queue()
        # Fields used for custom outfmt 6 BLAST output:
        self.fieldnames = ['query_id', 'subject_id', 'positives', 'mismatches', 'gaps',
                           'evalue', 'bit_score', 'subject_length', 'alignment_length',
                           'query_start', 'query_end', 'query_sequence',
                           'subject_start', 'subject_end', 'subject_sequence']
        # Set the report path
        self.reportpath = os.path.join(self.path, 'consolidated_report')
        self.report = os.path.join(self.reportpath, '{at}_report.csv'.format(at=self.analysistype))
        # The default length for the initial baiting - if there are primers shorter than this, then the shortest
        # value will be used
        self.klength = 20
        # A list of valid file extensions for FASTA formatted-files
        self.extensions = ['.fasta', '.fa', '.fas', '.fsa', '.fna', '.tfa', 'ffn']


if __name__ == '__main__':
    parser = ArgumentParser(description='Perform in silico PCR using bbduk and SPAdes')
    parser.add_argument('-p', '--path',
                        required=True,
                        help='Specify directory in which reports are to be created')
    parser.add_argument('-s', '--sequencepath',
                        required=True,
                        help='Path of folder containing .fasta/.fastq(.gz) files to process.')
    parser.add_argument('-n', '--cpus',
                        help='Number of threads. Default is the number of cores in the system')
    parser.add_argument('-pf', '--primerfile',
                        required=True,
                        help='Absolute path and name of the primer file (in FASTA format) to test. The file must have'
                             'every primer on a separate line AND -F/-R following the name e.g. '
                             '>primer1-F\n'
                             'ATCGACTGACAC....\n'
                             '>primer1-R\n'
                             'ATCGATCGATCGATG....\n'
                             '>primer2-F\n'
                             '.......\n')
    parser.add_argument('-m', '--mismatches',
                        default=1,
                        help='Number of mismatches allowed [0-3]. Default is 1')
    parser.add_argument('-k', '--kmerlength',
                        default='55,77,99,127',
                        help='The range of kmers used in SPAdes assembly. Default is 55,77,99,127, but you can '
                             'provide a comma-separated list of kmers e.g. 21,33,55,77,99,127 or a single kmer e.g. 33')
    parser.add_argument('-d', '--debug',
                        action='store_true',
                        help='Allow debug-level logging to be printed to the terminal')
    # Get the arguments into an object
    arguments = parser.parse_args()
    arguments.pipeline = False
    SetupLogging(debug=arguments.debug)
    # Define the start time
    arguments.start = time.time()

    # Run the script
    finder = PrimerFinder(args=arguments, 
                          analysistype='ePCR')
    # Run the script
    finder.main()
    logging.info('ePCR analyses complete')
