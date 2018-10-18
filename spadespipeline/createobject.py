#!/usr/bin/env python
from accessoryFunctions.accessoryFunctions import filer, GenObject, make_path, MetadataObject, relative_symlink
from glob import glob
import logging
import os
__author__ = 'adamkoziol'


class ObjectCreation(object):

    def createobject(self):
        logging.info('Finding sequence files')
        # Find all the .fastq files in the sequence path
        filelist = glob(os.path.join(self.sequencepath, '*.fastq*'))
        if filelist:
            self.extension = 'fastq'
            self.filehandler(filelist)
        else:
            filelist = glob(os.path.join(self.sequencepath, '*.fa*'))
            self.extension = 'fasta'
            self.filehandler(filelist)

    def filehandler(self, filelist):
        # Extract the base name of the globbed name + path provided
        if self.extension == 'fastq':
            names = map(lambda x: os.path.split(x)[1], filer(filelist))
        else:
            names = map(lambda x: os.path.split(x)[1].split('.')[0], filelist)
        # Iterate through the names of the fastq files
        for name in sorted(names):
            # Set the name
            metadata = MetadataObject()
            metadata.name = name
            # Set the destination folder
            outputdir = os.path.join(self.sequencepath, name)
            # Make the destination folder
            make_path(outputdir)
            # Get the files specific to the sequence name
            specific = glob(os.path.join(self.sequencepath, '{}*{}*'.format(name, self.extension)))
            # Initialise the general and commands GenObjects
            metadata.general = GenObject()
            metadata.commands = GenObject()
            metadata.run = GenObject()
            # Create the .fastqfiles category of :self.metadata
            metadata.general.fastqfiles = list()
            # Link the files to the output folder
            for fastq in sorted(specific):
                relative_symlink(fastq,
                                 outputdir)
                metadata.general.fastqfiles.append(os.path.join(outputdir, os.path.basename(fastq)))
            # Add the output directory to the metadata
            metadata.general.outputdirectory = outputdir
            try:
                metadata.general.bestassemblyfile = metadata.general.fastqfiles[0]
            except IndexError:
                pass
            # Find the data files corresponding to the sample
            datafiles = glob(os.path.join(self.datapath, '{}*.csv'.format(metadata.name)))
            # Assign attributes to the files depending on whether they are abundance files or not
            for datafile in datafiles:
                if 'abundance' in datafile:
                    metadata.general.abundancefile = datafile
                else:
                    metadata.general.assignmentfile = datafile

            # Append the metadata to the list of samples
            self.samples.append(metadata)

    def __init__(self, inputobject):
        self.samples = list()
        self.path = inputobject.path
        self.sequencepath = inputobject.sequencepath
        try:
            self.datapath = inputobject.datapath
        except AttributeError:
            self.datapath = None
        if self.datapath:
            assert os.path.isdir(self.datapath), 'Data location supplied is not a valid directory {0!r:s}' \
                .format(self.datapath)
        self.start = inputobject.start
        self.extension = ''
        self.createobject()
