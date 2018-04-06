#!/usr/bin/env python3
from accessoryFunctions.accessoryFunctions import printtime, get_version, run_subprocess, write_to_logfile, make_path
import threading
from threading import Thread
from queue import Queue
import os
__author__ = 'adamkoziol'


class Skesa(object):

    def main(self):
        self.assemble_threads()
        self.best_assemblyfile()

    def assemble_threads(self):
        # Only make as many threads are there are samples with fastq files
        for i in range(self.cpus):
            # Send the threads to the merge method. :args is empty as I'm using
            threads = Thread(target=self.assemble, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        for sample in self.metadata:
            # Initialise the assembly command
            sample.commands.assemble = str()
            try:
                if sample.general.trimmedcorrectedfastqfiles:
                    # Set the output directory
                    sample.general.assembly_output = os.path.join(sample.general.outputdirectory, 'assembly_output')
                    make_path(sample.general.assembly_output)
                    sample.general.assemblyfile = os.path.join(sample.general.assembly_output, '{name}_unfiltered.fasta'
                                                               .format(name=sample.name))
                    sample.general.bestassemblyfile = os.path.join(sample.general.assembly_output, '{name}.fasta'
                                                                   .format(name=sample.name))
                    fastqfiles = sample.general.trimmedcorrectedfastqfiles

                    # Set the the forward fastq files
                    sample.general.assemblyfastq = fastqfiles
                    forward = fastqfiles[0]
                    gz = True if '.gz' in forward else False
                    # If there are two fastq files
                    if len(fastqfiles) == 2:
                        # Set the reverse fastq name
                        sample.commands.assemble = 'skesa --fastq {fastqfiles} --cores {threads} --gz {gz} ' \
                                                   '--use_paired_ends --contigs_out {contigs}'\
                            .format(fastqfiles=','.join(fastqfiles),
                                    threads=self.threads,
                                    gz=gz,
                                    contigs=sample.general.assemblyfile)
                    # Same as above, but use single read settings for the assembler
                    else:
                        sample.commands.assemble = 'skesa --fastq {fastqfiles} --cores {threads} --gz {gz} ' \
                                                   '--contigs_out {contigs}'\
                            .format(fastqfiles=','.join(fastqfiles),
                                    threads=self.threads,
                                    gz=gz,
                                    contigs=sample.general.assemblyfile)
                # If there are no fastq files, populate the metadata appropriately
                else:
                    sample.general.assembly_output = 'NA'
                    sample.general.assemblyfastq = 'NA'
                    sample.general.bestassemblyfile = 'NA'
            except KeyError:
                sample.general.assembly_output = 'NA'
                sample.general.assemblyfastq = 'NA'
                sample.general.trimmedcorrectedfastqfiles = 'NA'
                sample.general.bestassemblyfile = 'NA'
            if sample.commands.assemble:
                # Put the arguments to pass to the assemble method into the queue
                self.assemblequeue.put(sample)
        self.assemblequeue.join()

    def assemble(self):
        while True:
            sample = self.assemblequeue.get()
            if not os.path.isfile(sample.general.assemblyfile):
                # Run the assembly
                out, err = run_subprocess(sample.commands.assemble)
                self.threadlock.acquire()
                write_to_logfile(sample.commands.assemble,
                                 sample.commands.assemble,
                                 self.logfile,
                                 sample.general.logout,
                                 sample.general.logerr,
                                 None,
                                 None)
                write_to_logfile(out,
                                 err,
                                 self.logfile,
                                 sample.general.logout,
                                 sample.general.logerr,
                                 None,
                                 None)
                self.threadlock.release()
            self.assemblequeue.task_done()

    def best_assemblyfile(self):
        """
        Determine whether the contigs.fasta output file from SPAdes is present. If not, set the .bestassembly
        attribute to 'NA'
        """
        for sample in self.metadata:
            # Set the name of the unfiltered spades assembly output file
            if os.path.isfile(sample.general.assemblyfile):
                sample.general.bestassemblyfile = sample.general.assemblyfile
            else:
                sample.general.bestassemblyfile = 'NA'
            # Set the name of the filtered assembly file
            filteredfile = os.path.join(sample.general.outputdirectory, '{}.fasta'.format(sample.name))
            # Add the name and path of the filtered file to the metadata
            sample.general.filteredfile = filteredfile

    def __init__(self, inputobject):
        self.metadata = inputobject.runmetadata.samples
        self.start = inputobject.starttime
        self.cpus = inputobject.cpus
        try:
            self.threads = int(self.cpus / len(self.metadata)) if self.cpus / len(self.metadata) > 1 else 1
        except TypeError:
            self.threads = self.cpus
        self.path = inputobject.path
        self.logfile = inputobject.logfile
        self.assemblequeue = Queue(maxsize=self.threads)
        self.threadlock = threading.Lock()
        printtime('Assembling sequences', self.start)
