#!/usr/bin/env python
from accessoryFunctions.accessoryFunctions import printtime, get_version, run_subprocess, write_to_logfile, dotter
import threading
from threading import Thread
from subprocess import call
from queue import Queue
import os
__author__ = 'adamkoziol'


class Spades(object):

    def spades(self):
        # Find the fastq files for each sample
        # Only make as many threads are there are samples with fastq files
        for i in range(len([sample.general for sample in self.metadata if type(sample.general.fastqfiles) is list])):
            # Send the threads to the merge method. :args is empty as I'm using
            threads = Thread(target=self.assemble, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        for sample in self.metadata:
            # Initialise the spades command
            spadescommand = ''
            # Split the string of the provided kmer argument
            kmerlist = self.kmers.split(',')
            # Regenerate the list of kmers to use if the kmer is less than the readlength
            # Had to type sample.run.forwardlength to an int, for some reason it was a string.
            sample.run.forwardlength = int(sample.run.forwardlength)
            sample.general.kmers = ','.join([kmer for kmer in kmerlist if int(kmer) <= sample.run.forwardlength])
            if sample.general.trimmedcorrectedfastqfiles:
                # Set the output directory
                sample.general.spadesoutput = os.path.join(sample.general.outputdirectory, 'spades_output')
                try:
                    fastqfiles = [sample.general.mergedreads,
                                  sample.general.unmergedforward,
                                  sample.general.unmergedreverse]
                    sample.general.assemblyfastq = fastqfiles
                    if os.path.isdir(sample.general.spadesoutput):
                        spadescommand = \
                            'spades.py -k {} --only-assembler --careful --continue -s {} -1 {} -2 {} -o {} -t {}'\
                            .format(sample.general.kmers,
                                    sample.general.mergedreads,
                                    sample.general.unmergedforward,
                                    sample.general.unmergedreverse,
                                    sample.general.spadesoutput,
                                    self.threads)
                    else:
                        spadescommand = 'spades.py -k {} --only-assembler --careful -s {} -1 {} -2 {} -o {} -t {}' \
                            .format(sample.general.kmers,
                                    sample.general.mergedreads,
                                    sample.general.unmergedforward,
                                    sample.general.unmergedreverse,
                                    sample.general.spadesoutput,
                                    self.threads)
                except AttributeError:
                    fastqfiles = sample.general.trimmedcorrectedfastqfiles
                    # Set the the forward fastq files
                    sample.general.assemblyfastq = fastqfiles
                    forward = fastqfiles[0]
                    # If there are two fastq files
                    if len(fastqfiles) == 2:
                        # Set the reverse fastq name
                        reverse = fastqfiles[1]
                        if sample.run.forwardlength < 50:
                            if os.path.isdir(sample.general.spadesoutput):
                                spadescommand = \
                                    'spades.py -k {} --only-assembler --careful --continue --s1 {} -o {} -t {}'\
                                    .format(sample.general.kmers, reverse, sample.general.spadesoutput,
                                            self.threads)
                            else:
                                spadescommand = 'spades.py -k {} --only-assembler --careful --s1 {} -o {} -t {}' \
                                    .format(sample.general.kmers, reverse, sample.general.spadesoutput,
                                            self.threads)
                        else:
                            # If a previous assembly was partially completed, continue from the most recent checkpoint
                            if os.path.isdir(sample.general.spadesoutput):
                                spadescommand = 'spades.py -k {} --only-assembler --careful --continue ' \
                                                '--pe1-1 {} --pe1-2 {} -o {} -t {}'\
                                                .format(sample.general.kmers, forward, reverse,
                                                        sample.general.spadesoutput,
                                                        self.threads)
                            else:
                                spadescommand = 'spades.py -k {} --only-assembler --careful ' \
                                                '--pe1-1 {} --pe1-2 {} -o {} -t {}'\
                                                .format(sample.general.kmers, forward, reverse,
                                                        sample.general.spadesoutput,
                                                        self.threads)
                    # Same as above, but use single read settings for spades
                    else:
                        if os.path.isdir(sample.general.spadesoutput):
                            spadescommand = 'spades.py -k {} --only-assembler --careful --continue --s1 {} -o {} -t {}'\
                                            .format(sample.general.kmers, forward, sample.general.spadesoutput,
                                                    self.threads)
                        else:
                            spadescommand = 'spades.py -k {} --only-assembler --careful --s1 {} -o {} -t {}'\
                                            .format(sample.general.kmers, forward, sample.general.spadesoutput,
                                                    self.threads)
            # If there are no fastq files, populate the metadata appropriately
            else:
                sample.general.spadesoutput = 'NA'
                sample.general.assemblyfastq = 'NA'
            # Put the arguments to pass to the assemble method into the queue
            self.assemblequeue.put((sample, spadescommand, sample.general.spadesoutput))
            # Add the command to the metadata
            sample.commands.spades = spadescommand
            # Add the version to the metadata
            sample.software.spades = get_version(['spades.py']).decode('utf-8').split('\n')[0].split()[3]
            # Set the name of the unfiltered spades assembly output file
            sample.general.bestassemblyfile = '{}/contigs.fasta'.format(sample.general.spadesoutput)
            # Set the name of the filtered assembly file
            filteredfile = os.path.join(sample.general.outputdirectory, '{}.fasta'.format(sample.name))
            # Add the name and path of the filtered file to the metadata
            sample.general.filteredfile = filteredfile
        # Join the threads
        self.assemblequeue.join()
        # Filter contigs shorter than 1000 bp, and rename remaining contigs with sample.name
        # self.filter()
        # self.insertsize()

    def assemble(self):
        """Run the assembly command in a multi-threaded fashion"""
        threadlock = threading.Lock()
        while True:
            (sample, command, output) = self.assemblequeue.get()
            if command and not os.path.isfile('{}/contigs.fasta'.format(output)):
                # execute(command)
                out, err = run_subprocess(command)
                threadlock.acquire()
                write_to_logfile(command, command, self.logfile, sample.general.logout, sample.general.logerr,
                                 None, None)
                write_to_logfile(out, err, self.logfile, sample.general.logout, sample.general.logerr, None, None)
                threadlock.release()
                #
                call(command, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
            dotter()
            # Signal to the queue that the job is done
            self.assemblequeue.task_done()

    def __init__(self, inputobject):
        self.metadata = inputobject.runmetadata.samples
        self.start = inputobject.starttime
        self.kmers = inputobject.kmers
        self.cpus = inputobject.cpus
        try:
            self.threads = int(self.cpus / len(self.metadata)) if self.cpus / len(self.metadata) > 1 else 1
        except TypeError:
            self.threads = self.cpus
        self.path = inputobject.path
        self.logfile = inputobject.logfile
        self.assemblequeue = Queue(maxsize=self.threads)
        printtime('Assembling sequences', self.start)
        self.spades()
