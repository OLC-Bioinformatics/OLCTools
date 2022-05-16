#!/usr/bin/env python3
import logging
import json
import os
__author__ = 'adamkoziol'


class MetadataPrinter(object):

    def printmetadata(self):
        # Iterate through each sample in the analysis
        for sample in self.metadata:
            # Set the name of the json file
            jsonfile = os.path.join(sample.general.outputdirectory, '{}_metadata.json'.format(sample.name))
            try:
                # Open the metadata file to write
                with open(jsonfile, 'w') as metadatafile:
                    # Write the json dump of the object dump to the metadata file
                    json.dump(sample.dump(), metadatafile, sort_keys=True, indent=4, separators=(',', ': '))
            except IOError:
                # Print useful information in case of an error
                logging.warning('Error creating .json file for {sample}'.format(sample=sample.name))
                raise
            except TypeError as e:
                logging.debug(f'Encountered TypeError writing metadata to file with the following details: {e}')

    def __init__(self, inputobject):
        try:
            self.metadata = inputobject.runmetadata.samples
        except AttributeError:
            try:
                self.metadata = inputobject.metadata.samples
            except AttributeError:
                try:
                    self.metadata = inputobject.metadata
                except AttributeError:
                    self.metadata = inputobject.runmetadata
        self.printmetadata()
