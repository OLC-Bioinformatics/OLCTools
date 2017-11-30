#!/usr/bin/env python
from glob import glob
import sipprCommon.runMetadata as runMetadata
from sipprCommon.offhours import Offhours
from accessoryFunctions.accessoryFunctions import make_path, printtime, GenObject, run_subprocess, write_to_logfile
import multiprocessing
import os
# Import ElementTree - try first to import the faster C version, if that doesn't
# work, try to import the regular version
try:
    import xml.etree.cElementTree as ElementTree
except ImportError:
    import xml.etree.ElementTree as ElementTree

__author__ = 'adamkoziol'


class CreateFastq(object):

    def createfastq(self):
        """Uses bcl2fastq to create .fastq files from a MiSeqRun"""
        from time import sleep
        from subprocess import call
        # Initialise samplecount
        samplecount = 0
        # If the fastq destination folder is not provided, make the default value of :path/:miseqfoldername
        self.fastqdestination = self.fastqdestination if self.fastqdestination else self.path + self.miseqfoldername
        # Make the path
        make_path(self.fastqdestination)
        # Initialise variables for storing index information
        index = ''
        indexlength = int()
        # A list of all projects supplied in the Sample_Project field in the SampleSheet.csv - these will be replaced
        projectlist = list()
        for strain in self.samples:
            projectlist.append(strain.run.SampleProject)
            samplecount += 1
            # Create a combined index of index1-index2
            try:
                strain.run.modifiedindex = '{}-{}'.format(strain.run.index, strain.run.index2)
                indexlength = 16
                index = 'I8,I8'
            except KeyError:
                strain.run.modifiedindex = strain.run.index
                indexlength = 6
                index = 'I6'
        # Create a new sample sheet using self.project name instead of the provided Sample_Project. This ensures
        # that all the FASTQ files are stored in the same output folder
        projectsamplesheet = os.path.join(self.fastqdestination, 'SampleSheet_modified.csv')
        with open(projectsamplesheet, "w") as modifiedsamplesheet:
            # Use the 'original' sample sheet as the template for the new sheet
            with open(self.customsamplesheet) as samplesheet:
                # Iterate through the template sheet, and write lines until the header for the data portion of the sheet
                for line in samplesheet:
                    modifiedsamplesheet.write(line)
                    if 'Sample_ID' in line:
                        for subline in samplesheet:
                            # Split the line on commas
                            data = subline.split(',')
                            updateddata = list()
                            # Find the Sample_Project value, and update it to be self.projectname
                            for value in data:
                                if value in projectlist:
                                    value = self.projectname
                                updateddata.append(value)
                            # Write the updated string to the new sheet
                            modifiedsamplesheet.write(','.join(updateddata))
        # Set :forward/reverse length to :header.forward/reverse length if the argument is not provided, or it's 'full',
        # otherwise  use the supplied argument
        self.forwardlength = self.metadata.header.forwardlength if self.forwardlength.lower()\
            == 'full' else self.forwardlength
        # Set :reverselength to :header.reverselength
        self.reverselength = self.metadata.header.reverselength if self.reverselength.lower() \
            == 'full' else self.reverselength
        # As the number of cycles required is the number of forward reads + the index(8) + the second index(8)
        # Also set the basemask variable as required
        if self.reverselength != '0':
            self.readsneeded = int(self.forwardlength) + int(self.reverselength) + indexlength
            basemask = "Y{}n*,{},Y{}n*".format(self.forwardlength, index, self.reverselength)
        else:
            self.readsneeded = int(self.forwardlength) + indexlength
            basemask = "Y{}n*,{},n*".format(self.forwardlength, index)
        # Handle plurality appropriately
        samples = 'samples' if samplecount != 1 else 'sample'
        number = 'are' if samplecount != 1 else 'is'
        printtime('There {} {} {} in this run.\n'
                  'MiSeqPath: {},\n'
                  'MiSeqFolder: {},\n'
                  'Fastq destination: {},\n'
                  'SampleSheet: {}'
                  .format(number, samplecount, samples, self.miseqpath, self.miseqfolder,
                          self.fastqdestination, projectsamplesheet),
                  self.start,
                  output=self.portallog)
        # Count the number of completed cycles in the run of interest
        cycles = glob(os.path.join(self.miseqpath, self.miseqfolder, 'Data', 'Intensities', 'BaseCalls', 'L001', 'C*'))
        while len(cycles) < self.readsneeded:
            printtime('Currently at {} cycles. Waiting until the MiSeq reaches cycle {}'.format(len(cycles),
                      self.readsneeded), self.start, output=self.portallog)
            sleep(300)
            cycles = glob(os.path.join(self.miseqpath, self.miseqfolder,
                                       'Data', 'Intensities', 'BaseCalls', 'L001', 'C*'))
        # configureBClToFastq requires :self.miseqfolder//Data/Intensities/BaseCalls/config.xml in order to work
        # When you download runs from BaseSpace, this file is not provided. There is an empty config.xml file that
        # can be populated with run-specific values and moved to the appropriate folder
        if not os.path.isfile(os.path.join(self.miseqfolder, 'Data', 'Intensities', 'BaseCalls', 'config.xml')):
            self.configfilepopulator()
        if self.debug:
            # Define the bcl2fastq system call for the unit test
            bclcall = "bcl2fastq --input-dir {basecalls} " \
                      "--output-dir {outdir} --sample-sheet {samplesheet} " \
                      "--barcode-mismatches 0 -r 4 -p {cpus} -w 4 -R {runfolder} --use-bases-mask {mask} " \
                      "--tiles s_1_1101 --minimum-trimmed-read-length 1" \
                .format(basecalls=os.path.join(self.miseqfolder, 'Data', 'Intensities', 'BaseCalls'),
                        outdir=self.fastqdestination,
                        samplesheet=projectsamplesheet,
                        cpus=self.cpus,
                        runfolder=self.miseqfolder,
                        mask=basemask)
        else:
            bclcall = "bcl2fastq --input-dir {basecalls} " \
                      "--output-dir {outdir} --sample-sheet {samplesheet} " \
                      "--barcode-mismatches 1 -r 4 -p {cpus} -w 4 -R {runfolder} --use-bases-mask {mask}"\
                .format(basecalls=os.path.join(self.miseqfolder, 'Data', 'Intensities', 'BaseCalls'),
                        outdir=self.fastqdestination,
                        samplesheet=projectsamplesheet,
                        cpus=self.cpus,
                        runfolder=self.miseqfolder,
                        mask=basemask)
        if not os.path.isdir(self.projectpath):
            # Call configureBclToFastq.pl
            printtime('Running bcl2fastq', self.start, output=self.portallog)
            # Run the command
            out, err = run_subprocess(bclcall)
            write_to_logfile(bclcall,
                             bclcall,
                             self.logfile)
            write_to_logfile(out,
                             err,
                             self.logfile)
        # Populate the metadata
        for sample in self.metadata.samples:
            sample.commands = GenObject()
            sample.commands.bcl = bclcall
            sample.run.forwardlength = self.forwardlength
            sample.run.reverselength = self.reverselength
        # Copy the fastq files to a central folder so they can be processed
        self.fastqmover()

    def configfilepopulator(self):
        """Populates an unpopulated config.xml file with run-specific values and creates
        the file in the appropriate location"""
        # Set the number of cycles for each read and index using the number of reads specified in the sample sheet
        self.forwardlength = int(self.metadata.header.forwardlength)
        self.reverselength = int(self.metadata.header.reverselength)
        # Create a list of lists containing [cycle start, cycle end, and :runid] for each of forward reads, index 1
        # index 2, and reverse reads
        cycles = [[1, self.forwardlength, self.runid],
                  [self.forwardlength + 1, self.forwardlength + 8, self.runid],
                  [self.forwardlength + 9, self.forwardlength + 16, self.runid],
                  [self.forwardlength + 17, self.forwardlength + 16 + self.reverselength, self.runid]]
        # A dictionary of parameters (keys) and the values to use when repopulating the config file
        parameters = {'RunFolder': self.runid, 'RunFolderDate': self.metadata.date.replace("-", ""),
                      'RunFolderId': self.metadata.runnumber, 'RunFlowcellId': self.metadata.flowcell}
        # Load the xml file using element tree
        config = ElementTree.parse(os.path.join(self.miseqpath, self.miseqfolder, 'Data', 'Intensities', 'BaseCalls',
                                                'config.xml'))
        # Get the root of the tree
        configroot = config.getroot()
        # The run node is the only child node of the root
        for run in configroot:
            # Iterate through the child nodes. There are three nodes sections that must be populated
            for child in run:
                # Find the cycles tag
                if child.tag == 'Cycles':
                    # Set the attributes with a dictionary containing the total reads
                    child.attrib = {'Last': '{}'.format(self.forwardlength + 16 + self.reverselength),
                                    'Number': '{}'.format(self.totalreads), 'First': '1'}
                elif child.tag == 'RunParameters':
                    # Name the child as runparameter for easier coding
                    runparameters = child
                    for runparameter in runparameters:
                        # This replaces data in both 'ImagingReads' and 'Reads' nodes
                        if 'Reads' in runparameter.tag:
                            # Enumerate through the run parameters
                            for indexcount, reads in enumerate(runparameter):
                                # The values for the index are 1, 2, 3, 4. Subtract one to get the index of the first
                                # list in cycles
                                index = int(runparameter.attrib['Index']) - 1
                                # Set the text value as the appropriate value from cycles
                                reads.text = str(cycles[index][indexcount])
                        # Populate the instrument value
                        if runparameter.tag == 'Instrument':
                            runparameter.text = self.instrument
                        # Iterate through the parameters in the parameter dictionary
                        for parameter in parameters:
                            # If the key is encountered
                            if runparameter.tag == parameter:
                                # Replace the text with the value
                                runparameter.text = parameters[parameter]
                        if 'Barcode' in runparameter.tag:
                            for cycle, barcode in enumerate(runparameter):
                                # Add the barcode cycles. These are the number of forward reads (+ 1 as the barcode
                                # starts 1 cycle after the first run) plus the current iterator
                                barcode.text = str(self.forwardlength + 1 + cycle)
        # Write the modified config file to the desired location
        config.write(os.path.join(self.miseqfolder, 'Data', 'Intensities', 'BaseCalls', 'config.xml'))

    def fastqmover(self):
        """Links .fastq files created above to :sequencepath"""
        import errno
        # Create the sequence path if necessary
        make_path(self.sequencepath)
        # Iterate through all the sample names
        for sample in self.metadata.samples:
            # Make directory variables
            outputdir = os.path.join(self.sequencepath, sample.name)
            # Glob all the .gz files in the subfolders - projectpath/Sample_:sample.name/*.gz
            for fastq in sorted(glob(os.path.join(self.projectpath, '*.gz'))):
                fastqname = os.path.basename(fastq)
                # Set the name of the destination file
                outputfile = os.path.join(self.sequencepath, fastqname)
                if not self.copy:
                    # Try/except loop link .gz files to self.path
                    try:
                        # Symlink fastq file to the seq path
                        relativepath = os.path.relpath(self.projectpath, self.sequencepath)
                        os.symlink(
                            os.path.join(relativepath, fastqname),
                            os.path.join(outputfile)
                        )
                    # Except os errors
                    except OSError as exception:
                        # If there is an exception other than the file exists, raise it
                        if exception.errno != errno.EEXIST:
                            raise
                else:
                    import shutil
                    # Copy the file if it doesn't already exist
                    if not os.path.isfile(outputfile):
                        shutil.copyfile(fastq, outputfile)
            # Repopulate .strainfastqfiles with the freshly-linked/copied files
            fastqfiles = glob(os.path.join(self.sequencepath, '{}*.fastq*'.format(sample.name)))
            fastqfiles = [fastq for fastq in fastqfiles if 'trimmed' not in fastq]
            # Populate the metadata object with the name/path of the fastq files
            sample.general.fastqfiles = fastqfiles
            # Save the outputdir to the metadata object
            sample.run.outputdirectory = outputdir
            sample.general.outputdirectory = outputdir
            sample.general.bestassemblyfile = True
            sample.general.trimmedcorrectedfastqfiles = sample.general.fastqfiles
            sample.general.logout = os.path.join(sample.general.outputdirectory, 'logout')
            sample.general.logerr = os.path.join(sample.general.outputdirectory, 'logerr')
            sample.commands = GenObject()

    def __init__(self, inputobject):
        """Initialise variables"""
        self.path = inputobject.path
        self.sequencepath = inputobject.sequencepath
        self.start = inputobject.starttime
        self.fastqdestination = inputobject.fastqdestination
        self.homepath = inputobject.homepath
        self.miseqout = str()
        self.projectname = 'fastqCreation'
        self.projectpath = os.path.join(self.fastqdestination, self.projectname)
        self.numreads = inputobject.numreads
        self.forwardlength = inputobject.forwardlength
        self.reverselength = inputobject.reverselength if self.numreads > 1 else '0'
        self.cpus = multiprocessing.cpu_count()
        try:
            self.logfile = inputobject.logfile
        except AttributeError:
            self.logfile = os.path.join(self.path, 'log')
        try:
            self.debug = inputobject.debug
        except AttributeError:
            self.debug = False
        self.readsneeded = 0
        self.commit = inputobject.commit
        self.copy = inputobject.copy
        try:
            self.portallog = inputobject.portallog
        except AttributeError:
            self.portallog = ''
        if inputobject.miseqpath:
            self.miseqpath = os.path.join(inputobject.miseqpath, '')
        else:
            print('MiSeqPath argument is required in order to use the fastq creation module. Please provide this '
                  'argument and run the script again.')
            quit()
        self.customsamplesheet = inputobject.customsamplesheet
        if self.customsamplesheet:
            assert os.path.isfile(self.customsamplesheet), 'Cannot find custom sample sheet as specified {}' \
                .format(self.customsamplesheet)
        # Use the assertions module from offhours to validate whether provided arguments are valid
        self.assertions = Offhours(inputobject)
        self.assertions.assertpathsandfiles()
        # Populate variables from this object
        self.miseqfolder = self.assertions.miseqfolder
        self.miseqfoldername = self.assertions.miseqfoldername
        self.customsamplesheet = self.assertions.customsamplesheet if self.assertions.customsamplesheet \
            else os.path.join(self.miseqfolder, 'SampleSheet.csv')
        self.runinfo = os.path.join(self.miseqfolder, 'RunInfo.xml')
        # Parse the sample sheet and other metadata files here
        self.metadata = runMetadata.Metadata(self)
        self.metadata.parseruninfo()
        # Create variables from this method
        self.flowcell = self.metadata.flowcell
        self.instrument = self.metadata.instrument
        self.samples = self.metadata.samples
        self.runid = self.metadata.runid
        # self.header = self.metadata.
        self.ids = self.metadata.ids
        self.date = self.metadata.date
        self.totalreads = self.metadata.totalreads
        # Create fastq files
        self.createfastq()


# If the script is called from the command line, then call the argument parser
if __name__ == '__main__':
    import subprocess
    from time import time
    # Get the current commit of the pipeline from git
    # Extract the path of the current script from the full path + file name
    homepath = os.path.split(os.path.abspath(__file__))[0]
    # Find the commit of the script by running a command to change to the directory containing the script and run
    # a git command to return the short version of the commit hash
    commit = subprocess.Popen('cd {} && git tag | tail -n 1'.format(homepath),
                              shell=True, stdout=subprocess.PIPE).communicate()[0].rstrip()
    from argparse import ArgumentParser
    # Parser for arguments
    parser = ArgumentParser(description='Assemble genomes from Illumina fastq files')
    parser.add_argument('-v', '--version',
                        action='version', version='%(prog)s commit {}'.format(commit))
    parser.add_argument('path',
                        help='Specify path')
    parser.add_argument('-n', '--numreads',
                        default=2,
                        type=int,
                        help='Specify the number of reads. Paired-reads:'
                        ' 2, unpaired-reads: 1. Default is paired-end')
    parser.add_argument('-t', '--threads',
                        help='Number of threads. Default is the number of cores in the system')
    parser.add_argument('-d', '--fastqdestination',
                        help='Optional folder path to store .fastq files created using the fastqCreation module. '
                             'Defaults to path/miseqfolder')
    parser.add_argument('-m', '--miseqpath',
                        required=True,
                        help='Path of the folder containing MiSeq run data folder e.g. /mnt/MiSeq')
    parser.add_argument('-f', '--miseqfolder',
                        required=True,
                        help='Name of the folder containing MiSeq run data e.g. 161129_M02466_0007_000000000-AW5L5')
    parser.add_argument('-r1', '--forwardlength',
                        default='full',
                        help='Length of forward reads to use. Can specify "full" to take the full length of forward '
                             'reads specified on the SampleSheet. Defaults to "full"')
    parser.add_argument('-r2', '--reverselength',
                        default='full',
                        help='Length of reverse reads to use. Can specify "full" to take the full length of reverse '
                             'reads specified on the SampleSheet. Defaults to "full"')
    parser.add_argument('-c', '--customsamplesheet',
                        help='Path of folder containing a custom sample sheet and name of sample sheet file '
                             'e.g. /home/name/folder/BackupSampleSheet.csv. Note that this sheet must still have the '
                             'same format of Illumina SampleSheet.csv files')
    parser.add_argument('-C', '--copy',
                        action='store_true',
                        help='Normally, the program will create symbolic links of the files into the sequence path, '
                             'however, the are occasions when it is necessary to copy the files instead')
    # Get the arguments into an object
    arguments = parser.parse_args()
    arguments.starttime = time()
    arguments.commit = commit
    arguments.homepath = homepath
    # Run the pipeline
    CreateFastq(arguments)
    # Print a bold, green exit statement
    print('\033[92m' + '\033[1m' + "\nElapsed Time: %0.2f seconds" % (time() - arguments.starttime) + '\033[0m')
