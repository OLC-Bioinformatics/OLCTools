#!/usr/bin/env python
from accessoryFunctions.accessoryFunctions import MetadataObject, GenObject
import copy
import os
# Import ElementTree - try first to import the faster C version, if that doesn't
# work, try to import the regular version
try:
    import xml.etree.cElementTree as ElementTree
except ImportError:
    import xml.etree.ElementTree as ElementTree
__author__ = 'adamkoziol'


class Metadata(object):

    def parseruninfo(self):
        """Extracts the flowcell ID, as well as the instrument name from RunInfo.xml. If this file is not provided,
        NA values are substituted"""
        # Check if the RunInfo.xml file is provided, otherwise, yield N/A
        try:
            runinfo = ElementTree.ElementTree(file=self.runinfo)
            # Get the run id from the
            for elem in runinfo.iter():
                for run in elem:
                    try:
                        self.runid = run.attrib['Id']
                        self.runnumber = run.attrib['Number']
                    except KeyError:
                        break
            # pull the text from flowcell and instrument values using the .iter(tag="X") function
            for elem in runinfo.iter(tag="Flowcell"):
                self.flowcell = elem.text
            for elem in runinfo.iter(tag="Instrument"):
                self.instrument = elem.text
        except IOError:
            pass
        # Extract run statistics from either GenerateRunStatistics.xml or indexingQC.txt
        self.parserunstats()

    def parsesamplesheet(self):
        """Parses the sample sheet (SampleSheet.csv) to determine certain values
        important for the creation of the assembly report"""
        # Open the sample sheet
        with open(self.samplesheet, "r") as samplesheet:
            # Iterate through the sample sheet
            samples, prev, header = False, 0, []
            for count, line in enumerate(samplesheet):
                # Remove new lines, and split on commas
                # line = line.decode('utf-8')  # Turn from bytes to string, since python3 is finicky.
                data = line.rstrip().split(",")
                if any(data):
                    if "[Settings]" in line:
                        samples = False
                    if not line.startswith("[") and not samples and not data == ['']:
                        # Grab an data not in the [Data] Section
                        setattr(self.header, data[0].replace(" ", ""), "".join(data[1:]))
                    elif "[Data]" in line or "[Reads]" in line:
                        samples = True
                    elif samples and "Sample_ID" in line:
                        header.extend([x.replace("_", "").replace(' ', "") for x in data])
                        prev = count
                    elif header:
                        # Try and replicate the Illumina rules to create file names from "Sample_Name"
                        samplename = samplenamer(data)
                        # Create an object for storing nested static variables
                        strainmetadata = MetadataObject()
                        # Set the sample name in the object
                        strainmetadata.name = samplename
                        # Add the header object to strainmetadata
                        # strainmetadata.__setattr__("run", GenObject(dict(self.header)))
                        strainmetadata.run = GenObject(copy.copy(self.header.datastore))
                        # Create the run object, so it will be easier to populate the object (eg run.SampleName = ...
                        # instead of strainmetadata.run.SampleName = ...
                        run = strainmetadata.run
                        # Capture Sample_ID, Sample_Name, I7_Index_ID, index1, I5_Index_ID,	index2, Sample_Project
                        for idx, item in enumerate(data):
                            setattr(run, header[idx], item) if item else setattr(run, header[idx], "NA")
                        # Add the sample number
                        run.SampleNumber = count - prev
                        # Create the 'General' category for strainmetadata
                        strainmetadata.general = GenObject({'outputdirectory': os.path.join(self.path, samplename),
                                                            'pipelinecommit': self.commit})
                        strainmetadata.general.logout = os.path.join(self.path, samplename,
                                                                     '{}_log_out.txt'.format(samplename))
                        strainmetadata.general.logerr = os.path.join(self.path, samplename,
                                                                     '{}_log_err.txt'.format(samplename))
                        # Add the output directory to the general category
                        # Append the strainmetadata object to a list
                        self.samples.append(strainmetadata)
                    elif samples:
                        setattr(self.header, 'forwardlength', data[0]) \
                            if 'forwardlength' not in self.header.datastore else \
                            setattr(self.header, 'reverselength', data[0])
                        self.totalreads += int(data[0])
        self.date = self.header.Date if "Date" in self.header.datastore else self.date
        for sample in self.samples:
            if 'InvestigatorName' not in sample.run.datastore:
                sample.run.InvestigatorName = 'NA'

    def parserunstats(self):
        """Parses the XML run statistics file (GenerateFASTQRunStatistics.xml). In some cases, the file is not
        available. Equivalent data can be pulled from Basespace.Generate a text file  name indexingQC.txt containing
        the copied tables from the Indexing QC tab of the run on Basespace"""
        # metadata = GenObject()
        # If the default file GenerateFASTQRunStatistics.xml is present, parse it
        if os.path.isfile(os.path.join(self.path, "GenerateFASTQRunStatistics.xml")):
            # Create a list of keys for which values are to be extracted
            datalist = ["SampleNumber", "SampleID", "SampleName", "NumberOfClustersPF"]
            # Load the file as an xml ElementTree object
            runstatistics = ElementTree.ElementTree(file=os.path.join(self.path, "GenerateFASTQRunStatistics.xml"))
            # Iterate through all the elements in the object
            # .iterfind() allow for the matching and iterating though matches
            # This is stored as a float to allow subsequent calculations
            tclusterspf = [float(element.text) for element in runstatistics.iterfind("RunStats/NumberOfClustersPF")][0]
            # Iterate through all the elements (strains) in the OverallSamples/SummarizedSampleStatistics category
            for element in runstatistics.iterfind("OverallSamples/SummarizedSampleStatistics"):
                # List comprehension. Essentially iterate through each element for each category in datalist:
                # (element.iter(category) and pull out the value for nestedelement
                straindata = [nestedelement.text for category in datalist for nestedelement in element.iter(category)]
                # Try and replicate the Illumina rules to create file names from "Sample_Name"
                samplename = samplenamer(straindata, 1)
                # Calculate the percentage of clusters associated with each strain
                # noinspection PyTypeChecker
                percentperstrain = "{:.2f}".format((float(straindata[3]) / tclusterspf * 100))
                try:
                    # Use the sample number -1 as the index in the list of objects created in parsesamplesheet
                    strainindex = int(straindata[0]) - 1
                    # Set run to the .run object of self.samples[index]
                    run = self.samples[strainindex].run
                    # An assertion that compares the sample computed above to the previously entered sample name
                    # to ensure that the samples are the same
                    assert self.samples[strainindex].name == samplename, \
                        "Sample name does not match object name {0!r:s}".format(straindata[1])
                    # Add the appropriate values to the strain metadata object
                    run.SampleNumber = straindata[0]
                    run.NumberofClustersPF = straindata[3]
                    run.TotalClustersinRun = tclusterspf
                    run.PercentOfClusters = percentperstrain
                    run.flowcell = self.flowcell
                    run.instrument = self.instrument
                except IndexError:
                    pass
        else:
            strainindex = 0
            for i in range(len(self.samples)):
                # Set run to the .run object of self.samples[index]
                run = self.samples[strainindex].run
                # Update the object with the variables
                run.SampleNumber = strainindex + 1
                run.NumberofClustersPF = 'NA'
                run.TotalClustersinRun = 'NA'
                run.PercentOfClusters = 'NA'
                run.flowcell = self.flowcell
                run.instrument = self.instrument
                strainindex += 1

    def __init__(self, passed):
        """Initialise variables"""
        self.path = passed.path
        self.runinfo = passed.runinfo
        self.flowcell = "NA"
        self.instrument = "NA"
        self.samples = list()
        self.ids = list()
        self.date = str()
        self.totalreads = 0
        self.runid = str()
        self.runnumber = str()
        self.commit = passed.commit

        # Create and start to populate the header object
        self.header = GenObject()
        # If a custom sample sheet has been provided, use it
        if passed.customsamplesheet:
            self.samplesheet = passed.customsamplesheet
            assert os.path.isfile(self.samplesheet), u'Could not find CustomSampleSheet as entered: {0!r:s}'\
                .format(self.samplesheet)
        else:
            self.samplesheet = os.path.join(self.path, "SampleSheet.csv")
        # Extract data from SampleSheet.csv
        self.parsesamplesheet()


def samplenamer(listofdata, indexposition=0):
    """Tries to replicate the Illumina rules to create file names from 'Sample_Name'
    :param listofdata: a list of data extracted from a file
    :param indexposition:
    """
    samplename = listofdata[indexposition].rstrip().replace(" ", "-").replace(".", "-").replace("=", "-")\
        .replace("+", "").replace("/", "-").replace("#", "").replace("---", "-").replace("--", "-")
    return samplename
