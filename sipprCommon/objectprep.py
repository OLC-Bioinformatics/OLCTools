#!/usr/bin/env python
from accessoryFunctions.accessoryFunctions import MetadataObject
import os
# import createFastq
import sipprCommon.fastqCreator as fastqCreator
import sipprCommon.createObject as createObject

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
            #
            self.samples = fastqCreator.CreateFastq(self)
            for key, value in vars(self.samples).items():
                if key == 'metadata':
                    self.forward = vars(value)['header'].forwardlength
                    self.reverse = vars(value)['header'].reverselength
                    self.header = vars(value)['header'].datastore
                elif key == 'samples':
                    self.index = value[0].run.modifiedindex
                    self.run = [x.run.datastore for x in value]

        else:
            self.samples = createObject.ObjectCreation(self)

    def __init__(self, inputobject):
        self.path = inputobject.path
        self.starttime = inputobject.starttime
        self.customsamplesheet = inputobject.customsamplesheet
        self.bcltofastq = inputobject.bcltofastq
        self.miseqpath = inputobject.miseqpath
        self.miseqfolder = inputobject.miseqfolder
        self.fastqdestination = inputobject.fastqdestination
        self.forwardlength = inputobject.forwardlength
        self.reverselength = inputobject.reverselength
        self.numreads = 2 if self.reverselength != 0 else 1
        self.customsamplesheet = inputobject.customsamplesheet
        self.sequencepath = inputobject.sequencepath
        self.homepath = inputobject.homepath
        self.commit = inputobject.commit
        self.copy = inputobject.copy
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
        self.header = dict()
        self.run = dict()
