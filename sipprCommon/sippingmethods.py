#!/usr/bin/env python
from glob import glob
from subprocess import call
from threading import Thread
from Bio.Sequencing.Applications import *
from accessoryFunctions.accessoryFunctions import *
from accessoryFunctions.metadataprinter import *
from sipprCommon.bowtie import *
from io import StringIO
__author__ = 'adamkoziol'


class Sippr(object):
    def targets(self):
        printtime('Performing analysis with {} targets folder'.format(self.analysistype), self.start)
        if self.pipeline:
            for sample in self.runmetadata:
                if sample.general.bestassemblyfile != 'NA':
                    setattr(sample, self.analysistype, GenObject())
                    # Set attributes
                    try:
                        sample[self.analysistype].targetpath = \
                            os.path.join(self.targetpath, self.analysistype, sample.mash.closestrefseqgenus, '')
                    except KeyError:
                        sample[self.analysistype].targetpath = \
                            os.path.join(self.targetpath, self.analysistype, sample.general.closestrefseqgenus, '')

                    # Ignore any species that do not match the desired species e.g. Listeria monocytogenes is acceptable
                    # while Listeria grayi is not. Right now, this sets the best assembly file to 'NA' to get the script
                    # to ignore this isolate, but something more fleshed out may be required in the future
                    for genus, species in self.taxonomy.items():
                        try:
                            if genus == sample.mash.closestrefseqgenus and species != sample.mash.closestrefseqspecies:
                                sample.general.bestassemblyfile = 'NA'
                        except KeyError:
                            pass
                    # There is a relatively strict databasing scheme necessary for the custom targets. Eventually,
                    # there will be a helper script to combine individual files into a properly formatted combined file
                    try:
                        sample[self.analysistype].baitfile = glob('{}*.fasta'
                                                                  .format(sample[self.analysistype].targetpath))[0]
                    # If the fasta file is missing, raise a custom error
                    except IndexError as e:
                        # noinspection PyPropertyAccess
                        e.args = ['Cannot find the combined fasta file in {}. Please note that the file must have a '
                                  '.fasta extension'.format(sample[self.analysistype].targetpath)]
                        if os.path.isdir(sample[self.analysistype].targetpath):
                            raise
                        else:
                            sample.general.bestassemblyfile = 'NA'

            for sample in self.runmetadata:
                if sample.general.bestassemblyfile != 'NA':
                    # Create the hash file of the baitfile
                    targetbase = sample[self.analysistype].baitfile.split('.')[0]
                    sample[self.analysistype].hashfile = targetbase + '.mhs.gz'
                    sample[self.analysistype].hashcall = 'cd {} && mirabait -b {} -k 19 -K {}'\
                        .format(sample[self.analysistype].targetpath,
                                sample[self.analysistype].baitfile,
                                sample[self.analysistype].hashfile)
                    if not os.path.isfile(sample[self.analysistype].hashfile):
                        call(sample[self.analysistype].hashcall, shell=True, stdout=self.devnull, stderr=self.devnull)
                    # Ensure that the hash file was successfully created
                    assert os.path.isfile(sample[self.analysistype].hashfile), \
                        u'Hashfile could not be created for the target file {0!r:s}'.format(
                            sample[self.analysistype].baitfile)
                    sample[self.analysistype].outputdir = os.path.join(sample.run.outputdirectory, self.analysistype)
                    sample[self.analysistype].baitedfastq = \
                        '{}/{}_targetMatches.fastq'.format(sample[self.analysistype].outputdir, self.analysistype)
        else:
            # There is a relatively strict databasing scheme necessary for the custom targets. Eventually, there will
            # be a helper script to combine individual files into a properly formatted combined file
            try:
                self.baitfile = glob('{}*.fasta'.format(self.targetpath))[0]
            # If the fasta file is missing, raise a custom error
            except IndexError:
                # Combine any .tfa files in the directory into a combined targets .fasta file
                from Bio import SeqIO
                tfafiles = glob(os.path.join(self.targetpath, '*.tfa'))
                if tfafiles:
                    with open(os.path.join(self.targetpath, 'combinedtargets.fasta'), 'wb') as combined:
                        for tfafile in tfafiles:
                            for record in SeqIO.parse(tfafile, 'fasta'):
                                SeqIO.write(record, combined, 'fasta')
                try:
                    self.baitfile = glob('{}*.fasta'.format(self.targetpath))[0]
                except IndexError as e:
                    # noinspection PyPropertyAccess
                    e.args = ['Cannot find the combined fasta file in {}. Please note that the file must have a '
                              '.fasta extension'.format(self.targetpath)]
                    raise
            # Create the hash file of the baitfile
            targetbase = self.baitfile.split('.')[0]
            self.hashfile = targetbase + '.mhs.gz'
            self.hashcall = 'cd {} && mirabait -b {} -k 19 -K {}'.format(self.targetpath, self.baitfile, self.hashfile)
            if not os.path.isfile(self.hashfile):
                call(self.hashcall, shell=True, stdout=self.devnull, stderr=self.devnull)
            # Ensure that the hash file was successfully created
            assert os.path.isfile(self.hashfile), u'Hashfile could not be created for the target file {0!r:s}' \
                .format(self.baitfile)
            for sample in self.runmetadata:
                setattr(sample, self.analysistype, GenObject())
                # Set attributes
                sample[self.analysistype].baitfile = self.baitfile
                sample[self.analysistype].hashfile = self.hashfile
                sample[self.analysistype].hashcall = self.hashcall
                sample[self.analysistype].targetpath = self.targetpath
                sample[self.analysistype].outputdir = os.path.join(sample.run.outputdirectory, self.analysistype)
                sample[self.analysistype].baitedfastq = '{}/{}_targetMatches.fastq'.format(sample[self.analysistype]
                                                                                           .outputdir,
                                                                                           self.analysistype)
        # Bait
        self.baiting()

    def baiting(self):
        # Perform baiting
        printtime('Performing kmer baiting of fastq files with {} targets'.format(self.analysistype), self.start)
        # Create and start threads for each fasta file in the list
        for i in range(len(self.runmetadata)):
            # Send the threads to the bait method
            threads = Thread(target=self.bait, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        for sample in self.runmetadata:
            if sample.general.bestassemblyfile != 'NA':
                # Add the sample to the queue
                self.baitqueue.put(sample)
        self.baitqueue.join()
        # Run the bowtie2 read mapping module
        self.mapping()

    def bait(self):
        """
        Runs mirabait on the fastq files
        """
        while True:
            sample = self.baitqueue.get()
            # Create the folder (if necessary)
            make_path(sample[self.analysistype].outputdir)
            # Make the system call
            if len(sample.general.fastqfiles) == 2:
                sample[self.analysistype].mirabaitcall = 'mirabait -c -B {} -t 4 -o {} -p {} {}' \
                    .format(sample[self.analysistype].hashfile, sample[self.analysistype].baitedfastq,
                            sample.general.fastqfiles[0], sample.general.fastqfiles[1])
            else:
                sample[self.analysistype].mirabaitcall = 'mirabait -c -B {} -t 4 -o {} {}' \
                    .format(sample[self.analysistype].hashfile, sample[self.analysistype].baitedfastq,
                            sample.general.fastqfiles[0])
            # Run the system call (if necessary)
            if not os.path.isfile(sample[self.analysistype].baitedfastq):
                call(sample[self.analysistype].mirabaitcall, shell=True, stdout=self.devnull, stderr=self.devnull)
            self.baitqueue.task_done()

    def mapping(self):
        printtime('Performing reference mapping', self.start)
        for i in range(len(self.runmetadata)):
            # Send the threads to
            threads = Thread(target=self.map, args=())
            # Set the daemon to True - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        for sample in self.runmetadata:
            if sample.general.bestassemblyfile != 'NA':
                # Set the path/name for the sorted bam file to be created
                sample[self.analysistype].sortedbam = '{}/{}_sorted.bam'.format(sample[self.analysistype].outputdir,
                                                                                self.analysistype)
                # Remove the file extension of the bait file for use in the indexing command
                sample[self.analysistype].baitfilenoext = sample[self.analysistype].baitfile.split('.')[0]
                # Use bowtie2 wrapper to create index the target file
                bowtie2build = Bowtie2BuildCommandLine(reference=sample[self.analysistype].baitfile,
                                                       bt2=sample[self.analysistype].baitfilenoext,
                                                       **self.builddict)
                # Use samtools wrapper to set up the bam sorting command
                samsort = SamtoolsSortCommandline(input=sample[self.analysistype].sortedbam,
                                                  o=True,
                                                  out_prefix="-")
                # Determine the location of the SAM header editing script
                import sipprCommon.editsamheaders
                scriptlocation = sipprCommon.editsamheaders.__file__
                samtools = [
                    # When bowtie2 maps reads to all possible locations rather than choosing a 'best' placement, the
                    # SAM header for that read is set to 'secondary alignment', or 256. Please see:
                    # http://davetang.org/muse/2014/03/06/understanding-bam-flags/ The script below reads in the stdin
                    # and subtracts 256 from headers which include 256
                    'python3 {}'.format(scriptlocation),
                    # Use samtools wrapper to set up the samtools view
                    SamtoolsViewCommandline(b=True,
                                            S=True,
                                            h=True,
                                            input_file="-"),
                    samsort]
                # Add custom parameters to a dictionary to be used in the bowtie2 alignment wrapper
                indict = {'--very-sensitive-local': True,
                          # For short targets, the match bonus can be increased
                          '--ma': self.matchbonus,
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
                # Add the commands to the queue. Note that the commands would usually be set as attributes of the sample
                # but there was an issue with their serialization when printing out the metadata
                if not os.path.isfile(sample[self.analysistype].baitfilenoext + '.1' + self.bowtiebuildextension):
                    stdoutbowtieindex, stderrbowtieindex = map(StringIO,
                                                               bowtie2build(cwd=sample[self.analysistype].targetpath))
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
                self.mapqueue.put((sample, bowtie2build, bowtie2align, samindex))
        self.mapqueue.join()
        # Use samtools to index the sorted bam file
        self.indexing()

    def map(self):
        while True:
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
            self.mapqueue.task_done()

    def indexing(self):
        printtime('Indexing sorted bam files', self.start)
        for i in range(len(self.runmetadata)):
            # Send the threads to
            threads = Thread(target=self.index, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        for sample in self.runmetadata:
            if sample.general.bestassemblyfile != 'NA':
                bamindex = SamtoolsIndexCommandline(input=sample[self.analysistype].sortedbam)
                sample[self.analysistype].sortedbai = sample[self.analysistype].sortedbam + '.bai'
                sample[self.analysistype].bamindex = str(bamindex)
                self.indexqueue.put((sample, bamindex))
        self.indexqueue.join()
        # Parse the results
        self.parsing()

    def index(self):
        while True:
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
            self.indexqueue.task_done()

    def parsing(self):
        printtime('Parsing sorted bam files', self.start)
        for i in range(len(self.runmetadata)):
            # Send the threads to
            threads = Thread(target=self.parse, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        for sample in self.runmetadata:
            if sample.general.bestassemblyfile != 'NA':
                # Get the fai file into a dictionary to be used in parsing results
                with open(sample[self.analysistype].faifile, 'r') as faifile:
                    for line in faifile:
                        data = line.split('\t')
                        try:
                            sample[self.analysistype].faidict[data[0]] = int(data[1])
                        except KeyError:
                            sample[self.analysistype].faidict = dict()
                            sample[self.analysistype].faidict[data[0]] = int(data[1])
                self.parsequeue.put(sample)
        self.parsequeue.join()

    def parse(self):
        import pysamstats
        import operator
        import numpy
        while True:
            sample = self.parsequeue.get()
            # Initialise dictionaries to store parsed data
            matchdict = dict()
            depthdict = dict()
            seqdict = dict()
            snpdict = dict()
            gapdict = dict()
            maxdict = dict()
            mindict = dict()
            deviationdict = dict()
            sample[self.analysistype].results = dict()
            sample[self.analysistype].avgdepth = dict()
            sample[self.analysistype].resultssnp = dict()
            sample[self.analysistype].resultsgap = dict()
            sample[self.analysistype].sequences = dict()
            sample[self.analysistype].maxcoverage = dict()
            sample[self.analysistype].mincoverage = dict()
            sample[self.analysistype].standarddev = dict()
            # Variable to store the expected position in gene/allele
            pos = 0
            try:
                # Use the stat_variation function of pysam stats to return records parsed from sorted bam files
                # Values of interest can be retrieved using the appropriate keys
                for rec in pysamstats.stat_variation(alignmentfile=sample[self.analysistype].sortedbam,
                                                     fafile=sample[self.analysistype].baitfile,
                                                     max_depth=1000000):
                    # Initialise seqdict with the current gene/allele if necessary with an empty string
                    if rec['chrom'] not in seqdict:
                        seqdict[rec['chrom']] = str()
                        # Since this is the first position in a "new" gene/allele, reset the pos variable to 0
                        pos = 0
                    # Initialise gap dict with 0 gaps
                    if rec['chrom'] not in gapdict:
                        gapdict[rec['chrom']] = 0
                    # If there is a gap in the alignment, record the size of the gap in gapdict
                    if int(rec['pos']) > pos:
                        # Add the gap size to gap dict
                        gapdict[rec['chrom']] += rec['pos'] - pos
                        # Set the expected position to the current position
                        pos = int(rec['pos'])
                    # Increment pos in preparation for the next iteration
                    pos += 1
                    # Initialise snpdict if necessary
                    if rec['chrom'] not in snpdict:
                        snpdict[rec['chrom']] = 0
                    # Initialise the current gene/allele in depthdict with the depth (reads_all) if necessary,
                    # otherwise add the current depth to the running total
                    if rec['chrom'] not in depthdict:
                        depthdict[rec['chrom']] = int(rec['reads_all'])
                    else:
                        depthdict[rec['chrom']] += int(rec['reads_all'])
                    # Dictionary of bases and the number of times each base was observed per position
                    bases = {'A': rec['A'], 'C': rec['C'], 'G': rec['G'], 'T': rec['T']}
                    # If the most prevalent base (calculated with max() and operator.itemgetter()) does not match the
                    # reference base, add this prevalent base to seqdict
                    if max(bases.items(), key=operator.itemgetter(1))[0] != rec['ref']:
                        seqdict[rec['chrom']] += max(bases.items(), key=operator.itemgetter(1))[0]
                        # Increment the running total of the number of SNPs
                        snpdict[rec['chrom']] += 1
                    else:
                        # If the bases match, add the reference base to seqdict
                        seqdict[rec['chrom']] += (rec['ref'])
                        # Initialise posdict if necessary, otherwise, increment the running total of matches
                        if rec['chrom'] not in matchdict:
                            matchdict[rec['chrom']] = 1
                        else:
                            matchdict[rec['chrom']] += 1
                    # Find the max and min coverage for each strain/gene combo
                    try:
                        maxdict[rec['chrom']] = int(rec['reads_all']) if \
                            int(rec['reads_all']) >= maxdict[rec['chrom']] else maxdict[rec['chrom']]
                    except KeyError:
                        maxdict[rec['chrom']] = int(rec['reads_all'])
                    try:
                        mindict[rec['chrom']] = int(rec['reads_all']) if \
                            int(rec['reads_all']) <= mindict[rec['chrom']] else mindict[rec['chrom']]
                    except KeyError:
                        mindict[rec['chrom']] = int(rec['reads_all'])
                    # Create a list of all the depths in order to calculate the standard deviation
                    try:
                        deviationdict[rec['chrom']].append(int(rec['reads_all']))
                    except KeyError:
                        deviationdict[rec['chrom']] = list()
                        deviationdict[rec['chrom']].append(int(rec['reads_all']))
            # If there are no results in the bam file, then pass over the strain
            except ValueError:
                pass
            # Iterate through all the genes/alleles with results above
            for allele in sorted(matchdict):
                # If the length of the match is greater or equal to the length of the gene/allele (multiplied by the
                # cutoff value) as determined using faidx indexing, then proceed
                if matchdict[allele] >= sample[self.analysistype].faidict[allele] * self.cutoff:
                    # Calculate the average depth by dividing the total number of reads observed by the
                    # length of the gene
                    averagedepth = float(depthdict[allele]) / float(matchdict[allele])
                    percentidentity = float(matchdict[allele]) / float(sample[self.analysistype].faidict[allele]) * 100
                    # Only report a positive result if this average depth is greater than 10X
                    if averagedepth > 10:
                        # Populate resultsdict with the gene/allele name, the percent identity, and the average depth
                        sample[self.analysistype].results.update({allele: '{:.2f}'.format(percentidentity)})
                        sample[self.analysistype].avgdepth.update({allele: '{:.2f}'.format(averagedepth)})
                        # Add the SNP and gap results to dictionaries
                        sample[self.analysistype].resultssnp.update({allele: snpdict[allele]})
                        sample[self.analysistype].resultsgap.update({allele: gapdict[allele]})
                        sample[self.analysistype].sequences.update({allele: seqdict[allele]})
                        sample[self.analysistype].maxcoverage.update({allele: maxdict[allele]})
                        sample[self.analysistype].mincoverage.update({allele: mindict[allele]})
                        sample[self.analysistype]\
                            .standarddev.update({allele: '{:.2f}'.format(numpy.std(deviationdict[allele], ddof=1))})
            self.parsequeue.task_done()

    # noinspection PyDefaultArgument
    def __init__(self, inputobject, cutoff=0.98, matchbonus=2, builddict=dict(), extension='.bt2'):
        from queue import Queue
        self.path = inputobject.path
        self.sequencepath = inputobject.sequencepath
        self.targetpath = inputobject.targetpath
        self.reportpath = inputobject.reportpath
        self.runmetadata = inputobject.runmetadata.samples
        self.start = inputobject.starttime
        self.analysistype = inputobject.analysistype
        self.cpus = inputobject.cpus
        self.threads = inputobject.threads
        self.pipeline = inputobject.pipeline
        self.homepath = inputobject.homepath
        self.taxonomy = inputobject.taxonomy
        self.cutoff = cutoff
        self.matchbonus = matchbonus
        self.builddict = builddict
        self.bowtiebuildextension = extension
        self.baitfile = str()
        self.hashfile = str()
        self.hashcall = str()
        self.devnull = open(os.devnull, 'wb')  # define /dev/null
        self.baitqueue = Queue(maxsize=self.cpus)
        self.mapqueue = Queue(maxsize=self.cpus)
        self.indexqueue = Queue(maxsize=self.cpus)
        self.parsequeue = Queue(maxsize=self.cpus)
        # Run the analyses
        self.targets()
        # Print the metadata
        printer = MetadataPrinter(self)
        printer.printmetadata()
