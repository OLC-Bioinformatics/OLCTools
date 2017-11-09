#!/usr/bin/env python
from glob import glob
from accessoryFunctions.accessoryFunctions import make_path
import os
import time
__author__ = 'adamkoziol'


class FastqMover(object):

    def movefastq(self):
        """Find .fastq files for each sample and move them to an appropriately named folder"""
        # from accessoryFunctions import relativesymlink
        print("\r[{:}] Moving fastq files".format(time.strftime("%H:%M:%S")))
        # Iterate through each sample
        for sample in self.metadata.runmetadata.samples:
            # Retrieve the output directory
            outputdir = '{}{}'.format(self.path, sample.name)
            # Find any fastq files with the sample name
            fastqfiles = sorted(glob('{}{}_*.fastq*'.format(self.path, sample.name))) \
                if sorted(glob('{}{}_*.fastq*'.format(self.path, sample.name))) \
                else sorted(glob('{}{}.fastq*'.format(self.path, sample.name))) \
                if sorted(glob('{}{}.fastq*'.format(self.path, sample.name))) \
                else sorted(glob('{}{}*.fastq*'.format(self.path, sample.name)))
            # Only try and move the files if the files exist
            if fastqfiles:
                make_path(outputdir)
                # Symlink the fastq files to the directory
                try:
                    list(map(lambda x: os.symlink('../{}'.format(os.path.basename(x)),
                                                  '{}/{}'.format(outputdir, os.path.basename(x))), fastqfiles))
                except OSError:
                    pass
                # Find any fastq files with the sample name
                fastqfiles = [fastq for fastq in sorted(glob('{}/{}*.fastq*'.format(outputdir, sample.name)))
                              if 'trimmed' not in fastq and 'normalised' not in fastq and 'corrected' not in fastq
                              and 'paired' not in fastq and 'unpaired' not in fastq]
            else:
                if outputdir:
                    # Find any fastq files with the sample name
                    fastqfiles = [fastq for fastq in sorted(glob('{}/{}*.fastq*'.format(outputdir, sample.name)))
                                  if 'trimmed' not in fastq and 'normalised' not in fastq and 'corrected' not in fastq
                                  and 'paired' not in fastq and 'unpaired' not in fastq]
            sample.general.fastqfiles = fastqfiles

    def __init__(self, inputobject):
        self.metadata = inputobject
        self.path = inputobject.path
        self.movefastq()
