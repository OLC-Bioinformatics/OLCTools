#!/usr/bin/env python3
from accessoryFunctions.accessoryFunctions import run_subprocess, write_to_logfile, make_path
from biotools import bbtools
from subprocess import CalledProcessError
from click import progressbar
import logging
import shutil
import os
__author__ = 'adamkoziol'


class Skesa(object):

    def main(self):
        self.skesa_assemble()
        self.best_assemblyfile()

    def skesa_assemble(self):
        """
        Run skesa to assemble genomes
        """
        with progressbar(self.metadata) as bar:
            for sample in bar:
                # Initialise the assembly command
                sample.commands.assemble = str()
                try:
                    if sample.general.trimmedcorrectedfastqfiles:
                        # If the sample is a pure isolate, assemble it. Otherwise, run the pre-metagenome pipeline
                        try:
                            status = sample.run.Description
                        except AttributeError:
                            status = 'unknown'
                        if status == 'metagenome':
                            self.merge(sample)
                        else:
                            # Set the output directory
                            sample.general.assembly_output = os.path.join(sample.general.outputdirectory,
                                                                          'assembly_output')
                            make_path(sample.general.assembly_output)
                            sample.general.assemblyfile = os.path.join(sample.general.assembly_output,
                                                                       '{name}_unfiltered.fasta'
                                                                       .format(name=sample.name))
                            sample.general.bestassemblyfile = os.path.join(sample.general.assembly_output,
                                                                           '{name}.fasta'
                                                                           .format(name=sample.name))
                            fastqfiles = sample.general.trimmedcorrectedfastqfiles

                            # Set the the forward fastq files
                            sample.general.assemblyfastq = fastqfiles
                            forward = fastqfiles[0]
                            gz = True if '.gz' in forward else False
                            # If there are two fastq files
                            if len(fastqfiles) == 2:
                                # Set the reverse fastq name https://github.com/ncbi/SKESA/issues/7
                                sample.commands.assemble = 'skesa --fastq {fastqfiles} --cores {threads} ' \
                                                           '--use_paired_ends --vector_percent 1 ' \
                                                           '--contigs_out {contigs}'\
                                    .format(fastqfiles=','.join(fastqfiles),
                                            threads=self.cpus,
                                            contigs=sample.general.assemblyfile)
                            # Same as above, but use single read settings for the assembler
                            else:
                                sample.commands.assemble = 'skesa --fastq {fastqfiles} --cores {threads} ' \
                                                           '--vector_percent 1 --contigs_out {contigs}'\
                                    .format(fastqfiles=','.join(fastqfiles),
                                            threads=self.cpus,
                                            contigs=sample.general.assemblyfile)
                    # If there are no fastq files, populate the metadata appropriately
                    else:
                        sample.general.assembly_output = 'NA'
                        sample.general.assemblyfastq = 'NA'
                        sample.general.bestassemblyfile = 'NA'
                except AttributeError:
                    sample.general.assembly_output = 'NA'
                    sample.general.assemblyfastq = 'NA'
                    sample.general.trimmedcorrectedfastqfiles = 'NA'
                    sample.general.bestassemblyfile = 'NA'
                if sample.commands.assemble and not os.path.isfile(sample.general.assemblyfile):
                    # Run the assembly
                    out, err = run_subprocess(sample.commands.assemble)
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

    def merge(self, sample):
        """
        Use bbmerge to merge paired FASTQ files for use in metagenomics pipelines. Create a report with the
        total number of reads, and the number of reads that could be paired
        :param sample: metadata sample object flagged as a metagenome
        """
        # Set the assembly file to 'NA' as assembly is not desirable for metagenomes
        sample.general.assemblyfile = 'NA'
        # Can only merge paired-end
        if len(sample.general.fastqfiles) == 2:
            outpath = os.path.join(sample.general.outputdirectory, 'merged_reads')
            make_path(outpath)
            # Merge path - keep all the merged FASTQ files in one directory
            merge_path = os.path.join(self.path, 'merged_reads')
            make_path(merge_path)
            # Set the name of the merged, and unmerged files
            sample.general.mergedreads = \
                os.path.join(merge_path, '{}_paired.fastq.gz'.format(sample.name))
            log = os.path.join(outpath, 'log')
            error = os.path.join(outpath, 'err')
            try:
                if not os.path.isfile(sample.general.mergedreads):
                    # Run the merging command
                    out, err, cmd = bbtools.bbmerge(forward_in=sorted(sample.general.trimmedcorrectedfastqfiles)[0],
                                                    merged_reads=sample.general.mergedreads,
                                                    mix=True,
                                                    returncmd=True,
                                                    threads=self.cpus)
                    write_to_logfile(out, err, self.logfile, sample.general.logout, sample.general.logerr, None, None)
                    with open(log, 'w') as log_file:
                        log_file.write(out)
                    with open(error, 'w') as error_file:
                        error_file.write(err)
            except (CalledProcessError, IndexError):
                delattr(sample.general, 'mergedreads')
            # Set the name of the report to store the metagenome file merging results
            report = os.path.join(self.reportpath, 'merged_metagenomes.csv')
            # Extract the total number of reads, and the number of reads that could be paired from the bbmerge
            # err stream
            num_reads, num_pairs = self.reads(error)
            # If the report doesn't exist, create it with the header and the results from the first sample
            if not os.path.isfile(report):
                with open(report, 'w') as report_file:
                    report_file.write('Sample,TotalReads,PairedReads\n{sample},{total},{paired}\n'
                                      .format(sample=sample.name,
                                              total=num_reads,
                                              paired=num_pairs))
            # If the report exists, open it to determine which samples have already been added - useful if re-running
            # the analysis
            else:
                lines = list()
                with open(report, 'r') as report_file:
                    for line in report_file:
                        lines.append(line.split(',')[0])
                # Add the results to the report
                if sample.name not in lines:
                    with open(report, 'a+') as report_file:
                        report_file.write('{sample},{total},{paired}\n'
                                          .format(sample=sample.name,
                                                  total=num_reads,
                                                  paired=num_pairs))

    @staticmethod
    def reads(err_log):
        """
        Parse the outputs from bbmerge to extract the total number of reads, as well as the number of reads that
        could be paired
        :param err_log: bbmerge outputs the stats in the error file
        :return: num_reads, the total number of reads, paired_reads, number of paired readds
        """
        # Initialise variables
        num_reads = 0
        paired_reads = 0
        # Open the log file
        with open(err_log, 'r') as error_log:
            # Extract the necessary information
            for line in error_log:
                if 'Pairs:' in line:
                    num_reads = line.split('\t')[-1].rstrip()
                elif 'Joined:' in line:
                    paired_reads = line.split('\t')[-2].rstrip()
        return num_reads, paired_reads

    def best_assemblyfile(self):
        """
        Determine whether the contigs.fasta output file from the assembler is present. If not, set the .bestassembly
        attribute to 'NA'
        """
        for sample in self.metadata:
            try:
                # Set the name of the filtered assembly file
                filtered_outputfile = os.path.join(self.path, 'raw_assemblies', '{}.fasta'.format(sample.name))
                # Set the name of the unfiltered spades assembly output file
                if os.path.isfile(sample.general.assemblyfile):
                    size = os.path.getsize(sample.general.assemblyfile)
                    # Ensure that the assembly isn't just an empty file
                    if size == 0:
                        sample.general.bestassemblyfile = 'NA'
                    else:
                        sample.general.bestassemblyfile = sample.general.assemblyfile
                        shutil.copyfile(sample.general.bestassemblyfile, filtered_outputfile)
                else:
                    sample.general.bestassemblyfile = 'NA'
                # Add the name and path of the filtered file to the metadata
                sample.general.filteredfile = filtered_outputfile
            except AttributeError:
                sample.general.assemblyfile = 'NA'
                sample.general.bestassemblyfile = 'NA'

    def __init__(self, inputobject):
        self.metadata = inputobject.runmetadata.samples
        self.start = inputobject.starttime
        self.cpus = inputobject.cpus
        self.path = inputobject.path
        self.logfile = inputobject.logfile
        self.reportpath = inputobject.reportpath
        make_path(os.path.join(self.path, 'BestAssemblies'))
        make_path(os.path.join(self.path, 'raw_assemblies'))
        make_path(self.reportpath)
        logging.info('Assembling sequences')
