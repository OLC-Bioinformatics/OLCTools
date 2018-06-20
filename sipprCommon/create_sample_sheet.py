#!/usr/bin/env python 3
from accessoryFunctions.accessoryFunctions import make_path
import os

__author__ = 'adamkoziol'


class SampleSheet(object):

    def samplesheet(self):
        """
        Create a custom sample sheet based on the original sample sheet for the run, but only including the samples
        that did not pass the quality threshold on the previous iteration
        """
        if self.demultiplex:
            make_path(self.samplesheetpath)
            self.customsamplesheet = os.path.join(self.samplesheetpath, 'SampleSheet.csv')
            header = ['Sample_ID', 'Sample_Name', 'Sample_Plate', 'Sample_Well', 'I7_Index_ID', 'index', 'I5_Index_ID',
                      'index2', 'Sample_Project', 'Description']
            with open(self.customsamplesheet, 'w') as samplesheet:
                lines = str()
                lines += '[Header]\n'
                lines += 'IEMFileVersion,{}\n'.format(self.header.IEMFileVersion)
                lines += 'Investigator Name,{}\n'.format(self.header.InvestigatorName)
                lines += 'Experiment Name,{}\n'.format(self.header.ExperimentName)
                lines += 'Date,{}\n'.format(self.header.Date)
                lines += 'Workflow,{}\n'.format(self.header.Workflow)
                lines += 'Application,{}\n'.format(self.header.Application)
                lines += 'Assay,{}\n'.format(self.header.Assay)
                lines += 'Description,{}\n'.format(self.header.Description)
                lines += 'Chemistry,{}\n'.format(self.header.Chemistry)
                lines += '\n'
                lines += '[Reads]\n'
                lines += str(self.forward) + '\n'
                lines += str(self.reverse) + '\n'
                lines += '\n'
                lines += '[Settings]\n'
                lines += 'ReverseComplement,{}\n'.format(self.header.ReverseComplement)
                lines += 'Adapter,{}\n'.format(self.header.Adapter)
                lines += '\n'
                lines += '[Data]\n'
                lines += ','.join(header)
                lines += '\n'
                # Correlate all the samples added to the list of incomplete samples with their metadata
                for incomplete in self.incomplete:
                    for sample in self.rundata:
                        if incomplete == sample['SampleID']:
                            # Use each entry in the header list as a key for the rundata dictionary
                            for data in header:
                                # Modify the key to be consistent with how the dictionary was populated
                                result = sample[data.replace('_', '')]
                                # Description is the final entry in the list, and shouldn't have a , following the value
                                if data != 'Description':
                                    lines += '{},'.format(result.replace('NA', ''))
                                # This entry should have a newline instead of a ,
                                else:
                                    lines += '{}\n'.format(result.replace('NA', ''))
                # Write the string to the sample sheet
                samplesheet.write(lines)

    def __init__(self, inputobject):
        self.forward = inputobject.forward
        self.reverse = inputobject.reverse
        self.incomplete = inputobject.incomplete
        self.header = inputobject.header
        self.rundata = inputobject.rundata
        self.samplesheetpath = inputobject.samplesheetpath
        self.customsamplesheet = str()
        self.demultiplex = inputobject.demultiplex
