#!/usr/bin/env python
from accessoryFunctions.accessoryFunctions import make_path, filer, MetadataObject, GenObject
from glob import glob
import errno
import os
__author__ = 'adamkoziol'


class ObjectCreation(object):
    def createobject(self):
        # Grab any .fastq files in the path
        fastqfiles = glob(os.path.join(self.path, '*.fastq*'))
        # Extract the base name of the globbed name + path provided
        fastqnames = map(lambda x: os.path.split(x)[1], filer(fastqfiles))
        # Iterate through the names of the fastq files
        for fastqname in sorted(fastqnames):
            # Set the name
            metadata = MetadataObject()
            metadata.name = fastqname
            # Set the destination folder
            outputdir = os.path.join(self.path, fastqname)
            # Make the destination folder
            make_path(outputdir)
            # Get the fastq files specific to the fastqname
            specificfastq = glob(os.path.join(self.path,'{}*.fastq*'.format(fastqname)))
            # Make relative symlinks to the files in :self.path
            try:
                for fastq in specificfastq:
                    # Get the basename of the file
                    fastqfile = os.path.split(fastq)[-1]
                    # Set the destination fastq path as the base name plus the destination folder
                    destinationfastq = os.path.join(outputdir, fastqfile)
                    # Symlink the files
                    os.symlink('../{}'.format(fastqfile), destinationfastq)
            # Except os errors
            except OSError as exception:
                # If there is an exception other than the file exists, raise it
                if exception.errno != errno.EEXIST:
                    raise
            # Initialise the general and run categories
            metadata.general = GenObject()
            metadata.run = GenObject()
            # Populate the .fastqfiles category of :self.metadata
            metadata.general.fastqfiles = [fastq for fastq in glob(
                os.path.join(outputdir, '{}*.fastq*'.format(fastqname))) if 'trimmed' not in fastq]
            # Add the output directory to the metadata
            metadata.general.outputdirectory = outputdir
            metadata.run.outputdirectory = outputdir
            metadata.general.bestassemblyfile = True
            metadata.general.trimmedcorrectedfastqfiles = metadata.general.fastqfiles
            metadata.general.logout = os.path.join(metadata.general.outputdirectory, 'logout')
            metadata.general.logerr = os.path.join(metadata.general.outputdirectory, 'logerr')
            # Initialise an attribute to store commands
            metadata.commands = GenObject()
            # Append the metadata to the list of samples
            self.samples.append(metadata)

    def __init__(self, inputobject):
        self.samples = list()
        self.path = inputobject.sequencepath
        try:
            self.portallog = inputobject.portallog
        except AttributeError:
            self.portallog = ''
        try:
            self.debug = inputobject.debug
        except AttributeError:
            self.debug = False
        self.createobject()
