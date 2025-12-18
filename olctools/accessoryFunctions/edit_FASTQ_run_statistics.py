#!/usr/bin/env python3
from olctools.accessoryFunctions.accessoryFunctions import filer, make_path, SetupLogging
from argparse import ArgumentParser
from Bio import SeqIO
from glob import glob
import logging
import gzip
import os

# Import ElementTree - try first to import the faster C version, if that doesn't
# work, try to import the regular version
try:
    import xml.etree.cElementTree as ElementTree
except ImportError:
    import xml.etree.ElementTree as ElementTree
__author__ = 'adamkoziol'


class EditFASTQRunStatistics(object):

    def main(self):
        if self.cycles:
            self.edit_config()
        self.parse_sample_sheet()
        self.parse_sample_statistics()
        if self.clusters:
            self.parse_run_statistics()
            self.locate_file()
            self.load_read_stats()
            self.recalculate_cluster_stats()
            self.update_cluster_stats()
            self.recalculate_sample_stats()
            self.update_read_stats()

    def edit_config(self):
        """
        Load the config.xml file, and replace the number of cycles as required. Note: the
        <Cycles First="1" Last="cycles" Number="cycles" />, as well as the RunFolder, Instrument, etc. must be
        replaced manually
        """
        # Load the file as an xml ElementTree object
        config = ElementTree.ElementTree(file=self.config_file)
        # Initialise variables to store the difference between the current number of cycles and the desired number of
        # cycles, as well as the the current number of cycles
        difference = 0
        initial = 0
        for element in config.iter():
            # The first read (of four) will be used to extract the necessary values
            if not difference:
                # Only need to worry about the 'LastCycle'
                if element.tag == 'LastCycle':
                    # Extract the original number of cycles in the file
                    initial = int(element.text)
                    # Calculate the difference between the current number of cycles and the desired number of cycles
                    difference = self.cycles - initial
                    # Replace the 'LastCycle' with the desired number of cycles
                    element.text = str(self.cycles)
            # All other reads can use the extracted values
            else:
                # Only edit the text in certain tags
                if element.tag in ['Cycle', 'FirstCycle', 'LastCycle']:
                    # Determine if the read is the fourth read
                    if int(element.text) == 2 * initial + 16:
                        # The fourth read will require twice the difference added to the text, as this is the second
                        # time the read length is different
                        element.text = str(2 * difference + int(element.text))
                    # Don't need to update the 'Cycle' or 'FirstCycle' if they are 1 or 2
                    elif element.tag in ['Cycle', 'FirstCycle'] and int(element.text) in [1, 2]:
                        pass
                    # All other cases, add the difference calculated above to the current cycle value
                    else:
                        element.text = str(difference + int(element.text))
        # Write the updated config.xml file
        config.write(self.replacement_config)

    def parse_sample_sheet(self):
        """
        Extract the values from the sample sheet
        """
        logging.info('Parsing SampleSheet: {ss}'.format(ss=self.sample_sheet))
        with open(self.sample_sheet, 'r') as sample_sheet:
            for line in sample_sheet:
                # Skip until the actual header line
                if 'Sample_ID' in line:
                    header = [entry for entry in line.rstrip().split(',')]
                    sample_count = 1
                    for subline in sample_sheet:
                        # Initialise the counter for the sample
                        self.sample_dict[sample_count] = dict()
                        # Capture Sample_ID, Sample_Name, I7_Index_ID, index1, I5_Index_ID,	index2, Sample_Project
                        for idx, item in enumerate(subline.rstrip().split(',')):
                            self.sample_dict[sample_count].update({header[idx].replace('_', ''): item})
                        sample_count += 1

    def parse_sample_statistics(self):
        """
        Parse the original GenerateFASTQRunStatistics.xml file to extract necessary values for downstream calculations
        """
        # Load the file as an xml ElementTree object
        runstatistics = ElementTree.ElementTree(file=self.stats_file)
        sample_number = 1
        # Iterate through all the elements (strains) in the OverallSamples/SummarizedSampleStatistics category
        for element in runstatistics.iterfind("OverallSamples/SummarizedSampleStatistics"):
            # Iterate through all the subelements in the element
            for subelement in element.iter():
                # Extract the values for the 'NumberOfClustersPF' and 'NumberOfClustersRaw'. Store the values in
                # the sample dictionary
                if subelement.tag == 'NumberOfClustersPF':
                    self.sample_dict[sample_number][subelement.tag] = subelement.text
                if subelement.tag == 'NumberOfClustersRaw':
                    self.sample_dict[sample_number][subelement.tag] = subelement.text
                # Update the 'SampleID' and 'SampleName' fields as required
                if subelement.tag == 'SampleID':
                    subelement.text = self.sample_dict[sample_number]['SampleName']
                if subelement.tag == 'SampleName':
                    subelement.text = self.sample_dict[sample_number]['SampleName']
                    sample_number += 1
        # Write the updated values to the edited file
        runstatistics.write(self.replacement_stats)

    def parse_run_statistics(self):
        """
        Extract the necessary values from the 'RunStats' block. These values include 'NumberOfClustersPF',
        'NumberOfClustersRaw', 'NumberOfUnalignedClusters', 'NumberOfUnalignedClustersPF', 'NumberOfUnindexedClusters',
        and 'NumberOfUnindexedClustersPF'
        """
        # Load the file as an xml ElementTree object
        runstatistics = ElementTree.ElementTree(file=self.stats_file)
        # Iterate through all the elements (strains) in the OverallSamples/SummarizedSampleStatistics category
        for element in runstatistics.iterfind("RunStats"):
            for subelement in element:
                # Only extract the text from the tags of interest
                if subelement.tag in self.cluster_stat_list:
                    self.cluster_stats[subelement.tag] = int(subelement.text)

    def locate_file(self):
        """
        Create a list of all the samples with FASTQ reads present in the supplied path
        """
        logging.info('Finding FASTQ files')
        files = glob(os.path.join(self.parent_path, '*.fastq.gz'))
        # Use the filer method to pair the files
        self.samples = sorted(filer(files))
        # Populate the sample_dict with the sample-specific files
        for count, sample in enumerate(self.samples):
            self.sample_dict[count + 1]['SampleFiles'] = sorted(glob(os.path.join('{sn}*.fastq.gz').format(sn=sample)))

    def load_read_stats(self):
        """
        Determine number of reads present in each pair. Create a summary file to avoid having to do this more than
        once
        """
        logging.info('Parsing read stats')
        # As this read counting takes a long time, only calculate it once - store the results in a file for subsequent
        # rapid extraction
        if not os.path.isfile(self.read_stats):
            for count, sample_dict in sorted(self.sample_dict.items()):
                # Initialise lists to store all the FASTQ read records
                self.sample_dict[count]['NumReads'] = int()
                # Open the compressed file with gzip - as the forward and reverse files have the same number of reads,
                # only need to process the forward reads
                with gzip.open(self.sample_dict[count]['SampleFiles'][0], 'rt') as reads:
                    # Append to the list of all the FASTQ reads using SeqIO
                    self.sample_dict[count]['NumReads'] += len(list(SeqIO.parse(reads, 'fastq')))
                # Open the read stats file to append the read information
                with open(self.read_stats, 'a+') as stats:
                    # Add the sample name and the number of reads to the file
                    stats.write('{sn},{reads}\n'.format(sn=self.sample_dict[count]['SampleName'],
                                                        reads=self.sample_dict[count]['NumReads']))
                # Update the total number of reads present - will be used downstream
                self.total_reads += int(self.sample_dict[count]['NumReads'])
        else:
            # Open the read stats file
            with open(self.read_stats, 'r') as stats:
                # Start the count at one to match the dictionary keys defined above
                count = 1
                # Iterate through the entire file
                for line in stats:
                    # Split the sample name from the number of reads
                    sample_name, num_reads = line.rstrip().split(',')
                    # Populate the dictionary with the sample-specific number of reads
                    self.sample_dict[count]['NumReads'] = int(num_reads)
                    # Update the total read count
                    self.total_reads += int(self.sample_dict[count]['NumReads'])
                    count += 1

    def recalculate_cluster_stats(self):
        """
        Use the run-specific number of reads to adjust the values extracted from the GenerateFASTQRunStatistics.xml
        file, so that the values match what is present in the files
        """
        # Determine the ratio of the number of reads in the simulated files against the total number of reads in the
        # run from which the file was generated - this ratio will be used to adjust all the other values
        total_ratio = self.total_reads / self.cluster_stats['NumberOfClustersPF']
        # Set the updated value for the NumberOfClustersPF as the number of simulated reads
        self.updated_cluster_stats['NumberOfClustersPF'] = self.total_reads
        # Multiply the other previous values by the ratio to calculate the run-specific values, and populate the
        # dictionary with these updated values
        self.updated_cluster_stats['NumberOfClustersRaw'] = \
            int(total_ratio * self.cluster_stats['NumberOfClustersRaw'])
        self.updated_cluster_stats['NumberOfUnalignedClusters'] = \
            int(total_ratio * self.cluster_stats['NumberOfUnalignedClusters'])
        self.updated_cluster_stats['NumberOfUnalignedClustersPF'] = \
            int(total_ratio * self.cluster_stats['NumberOfUnalignedClustersPF'])
        self.updated_cluster_stats['NumberOfUnindexedClusters'] = \
            int(total_ratio * self.cluster_stats['NumberOfUnindexedClusters'])
        self.updated_cluster_stats['NumberOfUnindexedClustersPF'] = \
            int(total_ratio * self.cluster_stats['NumberOfUnindexedClustersPF'])

    def update_cluster_stats(self):
        """
        Load the XML file, and update the necessary values within the 'RunStats' block
        """
        # Load the file as an xml ElementTree object
        runstatistics = ElementTree.ElementTree(file=self.replacement_stats)
        # Iterate through all the elements (strains) in the OverallSamples/SummarizedSampleStatistics category
        for element in runstatistics.iterfind("RunStats"):
            for subelement in element:
                if subelement.tag in self.cluster_stat_list:
                    # Replace the existing value with the string of the new values calculated above
                    subelement.text = str(self.updated_cluster_stats[subelement.tag])
        # Write the updated file to disk
        runstatistics.write(self.replacement_stats)

    def recalculate_sample_stats(self):
        """
        Recalculate the sample-specific values, 'NumberOfClustersRaw', and 'NumberOfClustersPF'
        """
        for count, sample_dict in self.sample_dict.items():
            # The new 'NumberOfClustersRaw' value is the sample-specific number of reads divided by the previous ratio
            # of 'NumberOfClustersRaw':'NumberOfClustersPF'
            sample_dict['NumberOfClustersRaw'] = \
                int(sample_dict['NumReads']
                    * (int(sample_dict['NumberOfClustersRaw'])
                       / int(sample_dict['NumberOfClustersPF'])))
            # Update the sample dictionary with the number of reads
            sample_dict['NumberOfClustersPF'] = sample_dict['NumReads']

    def update_read_stats(self):
        """
        Update the values for 'NumberOfClustersRaw', and 'NumberOfClustersPF' in the XML file
        """
        # Load the file as an xml ElementTree object
        runstatistics = ElementTree.ElementTree(file=self.replacement_stats)
        sample_number = 1
        # Iterate through all the elements (strains) in the OverallSamples/SummarizedSampleStatistics category
        for element in runstatistics.iterfind("OverallSamples/SummarizedSampleStatistics"):
            for subelement in element.iter():
                # Find the appropriate tags, and replace the text with the values calculated above
                if subelement.tag == 'NumberOfClustersPF':
                    subelement.text = str(self.sample_dict[sample_number][subelement.tag])
                if subelement.tag == 'NumberOfClustersRaw':
                    subelement.text = str(self.sample_dict[sample_number][subelement.tag])
                    sample_number += 1
        runstatistics.write(self.replacement_stats)

    def __init__(self, path, stats_file, sample_sheet, clusters, cycles):
        if path.startswith('~'):
            self.path = os.path.abspath(os.path.expanduser(os.path.join(path)))
        else:
            self.path = os.path.abspath(os.path.join(path))
        assert os.path.isdir(self.path), 'Supplied path {sp} does not exist!'.format(sp=self.path)
        self.stats_file = os.path.join(self.path, stats_file)
        self.edit_path = os.path.join(self.path, 'edited')
        make_path(self.edit_path)
        self.replacement_stats = os.path.join(self.edit_path, 'GenerateFASTQRunStatistics.xml')
        assert os.path.isfile(self.stats_file), 'Cannot locate stats file, {sf}, in supplied path {sp}'\
            .format(sf=stats_file,
                    sp=self.path)
        self.sample_sheet = os.path.join(self.path, sample_sheet)
        assert os.path.isfile(self.sample_sheet), 'Cannot locate sample sheet, {ss} in supplied path {sp}'\
            .format(ss=sample_sheet,
                    sp=self.path)
        self.clusters = clusters
        if self.clusters:
            self.parent_path = os.path.dirname(self.path)
            self.read_stats = os.path.join(self.edit_path, 'reads')
            self.total_reads = int()
        self.cycles = cycles
        if self.cycles:
            self.config_file = os.path.join(self.path, 'config.xml')
            self.replacement_config = os.path.join(self.edit_path, 'config.xml')
            assert os.path.isfile(self.config_file), 'Cannot locate {cf}, which is necessary when the "cycles" ' \
                                                     'argument is provided'.format(cf=self.config_file)
        self.cluster_stats = dict()
        self.updated_cluster_stats = dict()
        self.cluster_stat_list = ['NumberOfClustersPF', 'NumberOfClustersRaw', 'NumberOfUnalignedClusters',
                                  'NumberOfUnalignedClustersPF', 'NumberOfUnindexedClusters',
                                  'NumberOfUnindexedClustersPF']
        self.sample_dict = dict()
        self.samples = list()


def cli():
    # Parser for arguments
    parser = ArgumentParser(description='Edits a GenerateFASTQRunStatistics.xml based on a SampleSheet.csv file')
    parser.add_argument('-p', '--path',
                        required=True,
                        help='Path to folder containing sequencing files')
    parser.add_argument('-g', '--generatefastqrunstatisticsfile',
                        default='GenerateFASTQRunStatistics.xml',
                        help='Name of the GenerateFASTQRunStatistics.xml file to process. Must be in the supplied '
                             'sequence path. Default is GenerateFASTQRunStatistics.xml')
    parser.add_argument('-s', '--samplesheet',
                        default='SampleSheet.csv',
                        help='Name of the SampleSheet.csv file to use. Must be in the supplied sequence path. '
                             'Default is SampleSheet.csv')
    parser.add_argument('-cl', '--clusters',
                        action='store_true',
                        help='Use the FASTQ files present in the parent directory to update any "cluster"-related '
                             'entry')
    parser.add_argument('-cy', '--cycles',
                        type=int,
                        default=0,
                        help="Number of cycles to use to adjust the config.xml file")
    SetupLogging()
    # Get the arguments into an object
    arguments = parser.parse_args()
    rename = EditFASTQRunStatistics(path=arguments.path,
                                    stats_file=arguments.generatefastqrunstatisticsfile,
                                    sample_sheet=arguments.samplesheet,
                                    clusters=arguments.clusters,
                                    cycles=arguments.cycles)
    rename.main()


if __name__ == '__main__':
    cli()
