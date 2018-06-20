#!/usr/bin/env python
from accessoryFunctions.accessoryFunctions import MetadataObject
import sipprCommon.fastqCreator as fastqCreator
import sipprCommon.createObject as createObject
import os
__author__ = 'adamkoziol'


class Objectprep(object):

    def objectprep(self):
        """
        Creates fastq files from an in-progress Illumina MiSeq run or create an object and moves files appropriately
        """
        # Create .fastq files if necessary. Otherwise create the metadata object
        if self.bcltofastq:
            if self.customsamplesheet:
                assert os.path.isfile(self.customsamplesheet), 'Cannot find custom sample sheet as specified {}' \
                    .format(self.customsamplesheet)
            # Create the FASTQ files
            self.samples = fastqCreator.CreateFastq(self)
            # Create a dictionary of the object
            samples_dict = vars(self.samples)
            # Extract the required information from the dictionary
            self.index = samples_dict['index']
            self.index_length = samples_dict['indexlength']
            self.forward = samples_dict['forwardlength']
            self.reverse = samples_dict['reverselength']
            self.forwardlength = samples_dict['forward']
            self.reverselength = samples_dict['reverse']
            self.header = samples_dict['header']
        else:
            self.samples = createObject.ObjectCreation(self)

    def __init__(self, inputobject):
        self.path = inputobject.path
        self.starttime = inputobject.starttime
        self.sequencepath = inputobject.sequencepath
        try:
            self.customsamplesheet = inputobject.customsamplesheet
            self.bcltofastq = inputobject.bcltofastq
            self.miseqpath = inputobject.miseqpath
            self.miseqfolder = inputobject.miseqfolder
            self.fastqdestination = inputobject.fastqdestination
            self.forwardlength = inputobject.forwardlength
            self.reverselength = inputobject.reverselength
            self.numreads = 2 if self.reverselength != 0 else 1
            self.customsamplesheet = inputobject.customsamplesheet
            self.homepath = inputobject.homepath
            self.commit = inputobject.commit
            self.copy = inputobject.copy
            self.demultiplex = inputobject.demultiplex
        except AttributeError:
            self.bcltofastq = False
        try:
            self.debug = inputobject.debug
        except AttributeError:
            self.debug = False
        try:
            self.portallog = inputobject.portallog
        except AttributeError:
            self.portallog = ''
        self.samples = MetadataObject()
        self.forward = str()
        self.reverse = str()
        self.index = str()
        self.index_length = int()
        self.header = dict()
        self.run = dict()
