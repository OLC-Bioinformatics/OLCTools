#!/usr/bin/env python
import json
import os
__author__ = 'adamkoziol'


class MetadataPrinter(object):

    def printmetadata(self):
        # Iterate through each sample in the analysis
        for sample in self.metadata:
            if type(sample.general.fastqfiles) is list:
                # Set the name of the json file
                jsonfile = os.path.join(sample.general.outputdirectory, '{}_metadata.json'.format(sample.name))
                try:
                    # Open the metadata file to write
                    with open(jsonfile, 'w') as metadatafile:
                        # Write the json dump of the object dump to the metadata file
                        json.dump(sample.dump(), metadatafile, sort_keys=True, indent=4, separators=(',', ': '))
                except IOError:
                    # Print useful information in case of an error
                    print(sample.name, sample.datastore)
                    raise

    def __init__(self, inputobject):
        try:
            self.metadata = inputobject.runmetadata.samples
        except AttributeError:
            try:
                self.metadata = inputobject.metadata.samples
            except AttributeError:
                self.metadata = inputobject.metadata
        self.printmetadata()
