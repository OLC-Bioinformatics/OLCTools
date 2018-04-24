#!/usr/bin/env python
from accessoryFunctions.accessoryFunctions import make_path, printtime
from glob import glob
import os
__author__ = 'adamkoziol'


class FastqMover(object):

    def movefastq(self):
        """Find .fastq files for each sample and move them to an appropriately named folder"""
        printtime('Moving FASTQ files', self.start)
        # Iterate through each sample
        for sample in self.metadata.runmetadata.samples:
            # Retrieve the output directory
            outputdir = os.path.join(self.path, sample.name)
            # Find any fastq files with the sample name
            fastqfiles = sorted(glob(os.path.join(self.path, '{}_*.fastq*'.format(sample.name)))) \
                if sorted(glob(os.path.join(self.path, '{}_*.fastq*'.format(sample.name)))) \
                else sorted(glob(os.path.join(self.path, '{}.fastq*'.format(sample.name)))) \
                if sorted(glob(os.path.join(self.path, '{}.fastq*'.format(sample.name)))) \
                else sorted(glob(os.path.join(self.path, '{}*.fastq*'.format(sample.name))))
            # Only try and move the files if the files exist
            if fastqfiles:
                make_path(outputdir)
                # Symlink the fastq files to the directory
                try:
                    list(map(lambda x: os.symlink(os.path.join('..', os.path.basename(x)),
                                                  os.path.join(outputdir, os.path.basename(x))), fastqfiles))
                except OSError:
                    pass
                # Find any fastq files with the sample name
                fastqfiles = [fastq for fastq in sorted(glob(os.path.join(outputdir, '{}*.fastq*'.format(sample.name))))
                              if 'trimmed' not in fastq and 'normalised' not in fastq and 'corrected' not in fastq
                              and 'paired' not in fastq and 'unpaired' not in fastq]
            else:
                if outputdir:
                    # Find any fastq files with the sample name
                    fastqfiles = [fastq for fastq in sorted(glob(os.path.join(
                        outputdir, '{}*.fastq*'.format(outputdir, sample.name))))
                                  if 'trimmed' not in fastq and 'normalised' not in fastq and 'corrected' not in fastq
                                  and 'paired' not in fastq and 'unpaired' not in fastq]
            sample.general.fastqfiles = fastqfiles

    def __init__(self, inputobject):
        self.metadata = inputobject
        self.path = inputobject.path
        self.start = inputobject.starttime
        self.movefastq()
