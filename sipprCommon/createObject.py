#!/usr/bin/env python
from accessoryFunctions.accessoryFunctions import make_path, filer, GenObject, MetadataObject, relative_symlink
from glob import glob
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
            specificfastq = glob(os.path.join(self.path, '{}*.fastq*'.format(fastqname)))
            # Initialise the general and run categories
            metadata.general = GenObject()
            metadata.run = GenObject()
            # Create the .fastqfiles category of :self.metadata
            metadata.general.fastqfiles = list()
            # Link the files to the output folder
            for fastq in sorted(specificfastq):
                relative_symlink(fastq,
                                 outputdir)
                metadata.general.fastqfiles.append(os.path.join(outputdir, os.path.basename(fastq)))
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
