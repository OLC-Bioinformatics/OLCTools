#!/usr/bin/env python
# noinspection PyProtectedMember
from Bio.Application import _Option, AbstractCommandline, _Switch
from Bio.Alphabet import generic_dna, generic_protein
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from subprocess import Popen, PIPE, STDOUT
from collections import defaultdict
import subprocess
import datetime
import logging
import shutil
import errno
import shlex
import time
import glob
import os
import re
import sys

__author__ = 'adamkoziol', 'andrewlow'


def dependency_check(dependency):
    """
    Checks a program to see if it's installed (or at least, checks whether or not some sort of executable
    for it is on your path).
    :param dependency: Name of program you want to check, as a string.
    :return: True if dependency is present, False if it isn't.
    """
    check = shutil.which(dependency)
    if not check:
        return False
    else:
        return True


def find_paired_reads(fastq_directory, forward_id='_R1', reverse_id='_R2'):
    """
    Looks at a directory to try to find paired fastq files. Should be able to find anything fastq.
    :param fastq_directory: Complete path to directory containing fastq files.
    :param forward_id: Identifier for forward reads. Default R1.
    :param reverse_id: Identifier for reverse reads. Default R2.
    :return: List containing pairs of fastq files, in format [[forward_1, reverse_1], [forward_2, reverse_2]], etc.
    """
    pair_list = list()
    fastq_files = glob.glob(os.path.join(fastq_directory, '*.f*q*'))
    for name in fastq_files:
        if forward_id in name and os.path.isfile(name.replace(forward_id, reverse_id)):
            pair_list.append([name, name.replace(forward_id, reverse_id)])
    return pair_list


def find_unpaired_reads(fastq_directory, forward_id='_R1', reverse_id='_R2'):
    """
    Looks at a directory to try to find unpaired fastq files.
    :param fastq_directory: Complete path to directory containing fastq files.
    :param forward_id: Identifier for paired reads, forward.
    :param reverse_id: Identifier for paired reads, reverse.
    :return: List of paths to unpaired files.
    """
    unpaired_list = list()
    fastq_files = glob.glob(os.path.join(fastq_directory, '*.f*q*'))
    for name in fastq_files:
        if forward_id not in name and reverse_id not in name:
            unpaired_list.append(name)
        elif forward_id in name and not os.path.isfile(name.replace(forward_id, reverse_id)):
            unpaired_list.append(name)
        elif reverse_id in name and not os.path.isfile(name.replace(reverse_id, forward_id)):
            unpaired_list.append(name)
    return unpaired_list


def download_file(address, output_name, hour_start=18, hour_end=6, day_start=5, day_end=6, timeout=600):
    """
    Downloads a file, between specified hours. (Hour start has to be greater than hour end for this to work in current
    iteration).
    :param address: Address of file that you want to download.
    :param output_name: Where you want to save the file to.
    :param hour_start: Start of window where downloading is acceptable. Default 6PM (1800h)
    :param hour_end: End of window where downloading is acceptable. Default 6AM (600h)
    :param day_start: Start of window where it's always OK to download. Default Saturday (day 5).
    :param day_end: End of window where it's always OK to download. Default Sunday (day 6).
    :param timeout: How often to check if you're outside the acceptable download window (default 600 seconds).
    :return:
    """
    out = open(os.devnull, 'w')
    returncode = 28  # While loop is based on returncode given by curl, so need to initialize it to something.
    while returncode != 0:  # 0 means that the file has already been downloaded completely, so stop looping then.
        # Figure out what hour it is. If not in acceptable download window, wait a while before checking again.
        hour = datetime.datetime.now().time().hour
        minute = datetime.datetime.now().time().minute
        day = datetime.datetime.today().weekday()
        acceptable_hour = not(hour_end < hour < hour_start)  # True if current hour is between start and end.
        acceptable_day = day_start <= day <= day_end  # True if current day is a weekend day.
        if not(acceptable_hour or acceptable_day):
            print('Current time is {hour}:{minute}. I am not allowed to start downloading until'
                  ' {start_hour}:00.'.format(hour=hour, minute=minute, start_hour=hour_start))
            time.sleep(timeout)
        # If the file doesn't already exist, start downloading it.
        elif not os.path.exists(output_name):
            cmd = 'curl -o {outname} --max-time {timeout} {address}'.format(timeout=timeout,
                                                                            address=address,
                                                                            outname=output_name)
            returncode = subprocess.call(cmd, shell=True, stdout=out, stderr=out)
        # If the file does already exist, download it starting from filesize offset.
        else:
            file_size = os.path.getsize(output_name)
            cmd = 'curl -o {outname} --max-time {timeout} -C {file_size} {address}'.format(timeout=timeout,
                                                                                           address=address,
                                                                                           outname=output_name,
                                                                                           file_size=file_size)
            returncode = subprocess.call(cmd, shell=True, stdout=out, stderr=out)


def write_to_logfile(out, err, logfile, samplelog=None, sampleerr=None, analysislog=None, analysiserr=None):
    """
    Writes out and err (both should be strings) to logfile.
    """
    # Run log
    with open(logfile + '_out.txt', 'a+') as outfile:
        outfile.write(out + '\n')
    with open(logfile + '_err.txt', 'a+') as outfile:
        outfile.write(err + '\n')
    # Sample log
    if samplelog:
        with open(samplelog, 'a+') as outfile:
            outfile.write(out + '\n')
        with open(sampleerr, 'a+') as outfile:
            outfile.write(err + '\n')
    # Analysis log
    if analysislog:
        with open(analysislog, 'a+') as outfile:
            outfile.write(out + '\n')
        with open(analysiserr, 'a+') as outfile:
            outfile.write(err + '\n')


def clear_logfile(logfile):
    """
    As logfiles are appended to each time the same data are processed, sometimes it is desirable to clear out
    logsfiles from previous iterations
    :param logfile: Base name of logfile
    """
    try:
        os.remove(logfile + '_out.txt')
    except IOError:
        pass
    try:
        os.remove(logfile + '_err.txt')
    except IOError:
        pass


def run_subprocess(command):
    """
    command is the command to run, as a string.
    runs a subprocess, returns stdout and stderr from the subprocess as strings.
    """
    x = Popen(command, shell=True, stdout=PIPE, stderr=PIPE)
    out, err = x.communicate()
    out = out.decode('utf-8')
    err = err.decode('utf-8')
    return out, err


def get_version(exe):
    """
    :param exe: :type list required
    """
    assert isinstance(exe, list)
    return Popen(exe, stdout=PIPE, stderr=STDOUT).stdout.read()


def logstr(*args):
    yield "{}\n".__add__("-".__mul__(60)).__mul__(len(args)).format(*args)


def make_path(inpath):
    """
    from: http://stackoverflow.com/questions/273192/check-if-a-directory-exists-and-create-it-if-necessary \
    does what is indicated by the URL
    :param inpath: string of the supplied path
    """
    if not os.path.isfile(inpath):
        try:
            os.makedirs(inpath)
        except FileExistsError:
            pass
    else:
        raise OSError


def make_dict():
    """Makes Perl-style dictionaries"""
    return defaultdict(make_dict)


class CustomLogs(logging.StreamHandler):
    """
    Uses the logging module to create custom-coloured logs. The colours correspond to the level
    Modified from:
    http://uran198.github.io/en/python/2016/07/12/colorful-python-logging.html
    https://plumberjack.blogspot.com/2010/12/colorizing-logging-output-in-terminals.html
    """
    # Dictionary mapping logging level to colours
    level_map = {
        logging.DEBUG: '\033[1;92m',
        logging.INFO: '\033[1;94m',
        logging.ERROR: '\033[1;91m',
        logging.WARNING: '\033[1;93m',
        logging.CRITICAL: '\033[1;95m'
    }

    def emit(self, record):
        try:
            # Write the formatted record to the stream
            self.stream.write(self.format(record))
            self.stream.write(getattr(self, 'terminator', '\n'))
            # Flush the output to terminal
            self.flush()
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            self.handleError(record)

    def colorize(self, message, record):
        if record.levelno in self.level_map:
            # Extract the colour corresponding to the current level
            color = self.level_map[record.levelno]
            # Add the colour to the message. Reset the formatting with '\x1b[0m'
            message = ''.join((color, message, '\x1b[0m'))
        return message

    def format(self, record):

        message = logging.StreamHandler.format(self, record)
        parts = message.split('\n', 1)
        # Add the custom formatted date to the message
        parts[0] = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' ' + parts[0]
        parts[0] = self.colorize(parts[0], record)
        # Reconstitute the message from its updated parts
        message = '\n'.join(parts)
        return message


class SetupLogging(object):
    """
    Runs the CustomLogs class
    """

    def __init__(self, debug=False):
        # Create a logging object
        logger = logging.getLogger()
        # Set whether debug level messages should be displayed
        if debug:
            logger.setLevel(logging.DEBUG)
        else:
            logger.setLevel(logging.INFO)
        # Use CustomLogs to modify the handler
        # Only add a handler if it hasn't been added by another script in the pipeline
        if not logger.handlers:
            logger.addHandler(CustomLogs())


def printtime(string, start, option=None, output=None):
    """Prints a string with colour options with the elapsed time
    # Reset
    Color_Off='\033[0m'       # Text Reset

    # Regular Colors
    Black='\033[0;30m'        # Black
    Red='\033[0;31m'          # Red
    Green='\033[0;32m'        # Green
    Yellow='\033[0;33m'       # Yellow
    Blue='\033[0;34m'         # Blue
    Purple='\033[0;35m'       # Purple
    Cyan='\033[0;36m'         # Cyan
    White='\033[0;37m'        # White

    # Bold
    BBlack='\033[1;30m'       # Black
    BRed='\033[1;31m'         # Red
    BGreen='\033[1;32m'       # Green
    BYellow='\033[1;33m'      # Yellow
    BBlue='\033[1;34m'        # Blue
    BPurple='\033[1;35m'      # Purple
    BCyan='\033[1;36m'        # Cyan
    BWhite='\033[1;37m'       # White

    # Underline
    UBlack='\033[4;30m'       # Black
    URed='\033[4;31m'         # Red
    UGreen='\033[4;32m'       # Green
    UYellow='\033[4;33m'      # Yellow
    UBlue='\033[4;34m'        # Blue
    UPurple='\033[4;35m'      # Purple
    UCyan='\033[4;36m'        # Cyan
    UWhite='\033[4;37m'       # White

    # Background
    On_Black='\033[40m'       # Black
    On_Red='\033[41m'         # Red
    On_Green='\033[42m'       # Green
    On_Yellow='\033[43m'      # Yellow
    On_Blue='\033[44m'        # Blue
    On_Purple='\033[45m'      # Purple
    On_Cyan='\033[46m'        # Cyan
    On_White='\033[47m'       # White

    # High Intensity
    IBlack='\033[0;90m'       # Black
    IRed='\033[0;91m'         # Red
    IGreen='\033[0;92m'       # Green
    IYellow='\033[0;93m'      # Yellow
    IBlue='\033[0;94m'        # Blue
    IPurple='\033[0;95m'      # Purple
    ICyan='\033[0;96m'        # Cyan
    IWhite='\033[0;97m'       # White

    # Bold High Intensity
    BIBlack='\033[1;90m'      # Black
    BIRed='\033[1;91m'        # Red
    BIGreen='\033[1;92m'      # Green
    BIYellow='\033[1;93m'     # Yellow
    BIBlue='\033[1;94m'       # Blue
    BIPurple='\033[1;95m'     # Purple
    BICyan='\033[1;96m'       # Cyan
    BIWhite='\033[1;97m'      # White

    # High Intensity backgrounds
    On_IBlack='\033[0;100m'   # Black
    On_IRed='\033[0;101m'     # Red
    On_IGreen='\033[0;102m'   # Green
    On_IYellow='\033[0;103m'  # Yellow
    On_IBlue='\033[0;104m'    # Blue
    On_IPurple='\033[0;105m'  # Purple
    On_ICyan='\033[0;106m'    # Cyan
    On_IWhite='\033[0;107m'   # White
    :param string: a string to be printed
    :param start: integer of the starting time
    :param option: Additional option for the text style
    :param output: name and path of the logfile to store the message
    """
    # If not option is provided, default to bold high-intensity white
    if not option:
        # option = '\033[1;97m'
        option = '\033[1;94m'
    # Add the string formatting option to the message. Reset the format back to normal at the end with \033[0m
    print('{} [Elapsed Time: {:.2f} seconds] {} \033[0m'.format(option, time.time() - start, string))
    if output:
        try:
            with open(output, 'a') as log:
                log.write('[Elapsed Time: {:.2f} seconds] {}\n'.format(time.time() - start, string))
        except FileNotFoundError:
            pass


class Dotter(object):

    def globalcounter(self):
        """Resets the globalcount to 0"""
        self.globalcount = 0

    def dotter(self):
        """Prints formatted time to stdout at the start of a line, as well as a "."
        whenever the length of the line is equal or lesser than 80 "." long"""
        if self.globalcount <= 80:
            sys.stdout.write('.')
            self.globalcount += 1
        else:
            sys.stdout.write('\n.')
            self.globalcount = 1

    def __init__(self):
        self.globalcount = 0


# Initialise globalcount
globalcount = 0


def globalcounter():
    """Resets the globalcount to 0"""
    global globalcount
    globalcount = 0


def dotter():
    """Prints formatted time to stdout at the start of a line, as well as a "."
    whenever the length of the line is equal or lesser than 80 "." long"""
    # Use a global variable
    global globalcount
    if globalcount <= 80:
        sys.stdout.write('.')
        globalcount += 1
    else:
        sys.stdout.write('\n.')
        globalcount = 1


def execute(command, outfile=""):
    """
    Allows for dots to be printed to the terminal while waiting for a long system call to run
    :param command: the command to be executed
    :param outfile: optional string of an output file
    from https://stackoverflow.com/questions/4417546/constantly-print-subprocess-output-while-process-is-running
    """
    # Initialise count
    count = 0
    # Initialise the starting time
    start = int(time.time())
    maxtime = 0
    # Removing Shell=True to prevent excess memory use thus shlex split if needed
    if type(command) is not list:
        command = shlex.split(command)
    # Run the commands - direct stdout to PIPE and stderr to stdout
    # DO NOT USE subprocess.PIPE if not writing it!
    if outfile:
        process = Popen(command, stdout=PIPE, stderr=STDOUT)
    else:
        devnull = open(os.devnull, 'wb')
        process = Popen(command, stdout=devnull, stderr=STDOUT)
    # Write the initial time
    sys.stdout.write('[{:}] '.format(time.strftime('%H:%M:%S')))
    # Create the output file - if not provided, then nothing should happen
    writeout = open(outfile, "ab+") if outfile else ""
    # Poll process for new output until finished
    while True:
        # If an output file name is provided
        if outfile:
            # Get stdout into a variable
            nextline = process.stdout.readline()
            # Print stdout to the file
            writeout.write(nextline)
        # Break from the loop if the command is finished
        if process.poll() is not None:
            break
        # Adding sleep commands slowed down this method when there was lots of output. Difference between the start time
        # of the analysis and the current time. Action on each second passed
        currenttime = int(time.time())
        if currenttime - start > maxtime:
            # Set the max time for each iteration
            maxtime = currenttime - start
            # Print up to 80 dots on a line, with a one second delay between each dot
            if count <= 80:
                sys.stdout.write('.')
                count += 1
            # Once there are 80 dots on a line, start a new line with the the time
            else:
                sys.stdout.write('\n[{:}] .'.format(time.strftime('%H:%M:%S')))
                count = 1
    # Close the output file
    writeout.close() if outfile else ""
    sys.stdout.write('\n')


def filer(filelist, extension='fastq'):
    """
    Helper script that creates a set of the stain names created by stripping off parts of the filename.
    Hopefully handles different naming conventions (e.g. 2015-SEQ-001_S1_L001_R1_001.fastq(.gz),
    2015-SEQ-001_R1_001.fastq.gz, 2015-SEQ-001_R1.fastq.gz, 2015-SEQ-001_1.fastq.gz, and 2015-SEQ-001_1.fastq.gz
    all become 2015-SEQ-001)
    :param filelist: List of files to parse
    :param extension: the file extension to use. Default value is 'fastq
    """
    # Initialise the set
    fileset = set()
    for seqfile in filelist:
        # Search for the conventional motifs present following strain names
        # _S\d+_L001_R\d_001.fastq(.gz) is a typical unprocessed Illumina fastq file
        if re.search("_S\d+_L001", seqfile):
            fileset.add(re.split("_S\d+_L001", seqfile)[0])
        # Files with _R\d_001.fastq(.gz) are created in the SPAdes assembly pipeline
        elif re.search("_R\d_001", seqfile):
            fileset.add(re.split("_R\d_001", seqfile)[0])
        # _R\d.fastq(.gz) represents a simple naming scheme for paired end reads
        elif re.search("R\d.{}".format(extension), seqfile):
            fileset.add(re.split("_R\d.{}".format(extension), seqfile)[0])
        # _\d.fastq is always possible
        elif re.search("[-_]\d.{}".format(extension), seqfile):
            fileset.add(re.split("[-_]\d.{}".format(extension), seqfile)[0])
        # .fastq is the last option
        else:
            fileset.add(re.split(".{}".format(extension), seqfile)[0])
    return fileset


def relativesymlink(src_file, dest_file):
    """
    https://stackoverflow.com/questions/9793631/creating-a-relative-symlink-in-python-without-using-os-chdir
    :param src_file: the file to be linked
    :param dest_file: the path and filename to which the file is to be linked
    """
    # Perform relative symlinking
    try:
        print(os.path.relpath(src_file), os.path.relpath(dest_file))
        os.symlink(
            # Find the relative path for the source file and the destination file
            os.path.relpath(src_file),
            os.path.relpath(dest_file)
        )
    # Except os errors
    except OSError as exception:
        # If the os error is anything but directory exists, then raise
        if exception.errno != errno.EEXIST:
            raise


def relative_symlink(src_file, output_dir, output_name=None):
    """
    Create relative symlinks files - use the relative path from the desired output directory to the storage path
    e.g. ../../2013-SEQ-0072/simulated/40/50_150/simulated_trimmed/2013-SEQ-0072_simulated_40_50_150_R1.fastq.gz
    is the relative path to the output_dir. The link name is the base name of the source file joined to the desired
    output directory e.g. output_dir/2013-SEQ-0072/2013-SEQ-0072_simulated_40_50_150_R1.fastq.gz
    https://stackoverflow.com/questions/9793631/creating-a-relative-symlink-in-python-without-using-os-chdir
    :param src_file: Source file to be symbolically linked
    :param output_dir: Destination folder for the link
    :param output_name: Optionally allow for the link to have a different name
    """
    if output_name:
        file_name = output_name
    else:
        file_name = os.path.basename(src_file)
    try:
        os.symlink(
            os.path.relpath(
                src_file,
                output_dir),
            os.path.join(
                output_dir,
                file_name
            )
        )
    # Ignore FileExistsErrors
    except FileExistsError:
        pass


class GenObject(object):
    """Object to store static variables"""
    def __init__(self, x=None):
        start = x if x else {}
        super(GenObject, self).__setattr__('datastore', start)

    def __getattr__(self, key):
        if key in self.datastore:
            return self.datastore[key]
        else:
            raise AttributeError('The GenObject has not been initialised with the following key: {key}'
                                 .format(key=key))

    def __setattr__(self, key, value):
        try:
            self.datastore[key] = value
        except TypeError:
            raise AttributeError('The GenObject cannot accept the following key:value pair provided {key}:{value}'
                                 .format(key=key,
                                         value=value))

    def __delattr__(self, key):
        try:
            del self.datastore[key]
        except KeyError:
            raise AttributeError('The GenObject does not contain the following key: {key}'
                                 .format(key=key))

    def returnattr(self, key):
        """
        Returns a string of either datastore[key], or 'ND' if datastore[key] doesn't exist formatted for a CSV report
        Replace any commas with semicolons.
        :param key: Dictionary key to be used to return the value from datastore[key]
        """
        try:
            if key in self.datastore:
                # Return the string of the value with any commas replaced by semicolons. Append a comma to the
                # end of the string for the CSV format
                return '{},'.format(str(self.datastore[key]).replace(',', ';'))
            else:
                return 'ND,'
        except AttributeError:
            return 'ND,'

    def isattr(self, key):
        """
        Checks to see if an attribute exists. If it does, returns True, otherwise returns False
        :param key: Dictionary key to be checked for presence in the datastore
        :return: True/False depending on whether an attribute exists
        """
        try:
            if key in self.datastore:
                return True
            else:
                return False
        except AttributeError:
            return False

    def __getitem__(self, key):
        """
        Make GenObjects subscriptable in order to allow for nested GenObject
        """
        return getattr(self, key)

    def __setitem__(self, key, value):
        """
        Allow item assignment for GenObjects
        """
        try:
            self.datastore[key] = value
        except TypeError:
            raise AttributeError('The GenObject cannot accept the following key:value pair provided {key}:{value}'
                                 .format(key=key,
                                         value=value))


class MetadataObject(object):
    """Object to store static variables"""
    def __init__(self):
        """Create datastore attr with empty dict"""
        super(MetadataObject, self).__setattr__('datastore', {})
        self.unwanted_keys = ['allelenames', 'alleles', 'faidict', 'gaplocations', 'maxcoverage',
                              'mincoverage', 'profiledata', 'resultsgap', 'averagedepth', 'avgdepth',
                              'resultssnp', 'sequences', 'sequence', 'snplocations', 'standarddev',
                              'totaldepth']

    def __getattr__(self, key):
        """:key is retrieved from datastore if exists, for nested attr recursively :self.__setattr__"""
        try:
            return self.datastore[key]
        except KeyError:
            raise AttributeError('The MetadataObject has not been initialised with the following key: {key}'
                                 .format(key=key))

    def __setattr__(self, key, value=GenObject(), **args):
        """Add :value to :key in datastore or create GenObject for nested attr"""
        if args:
            self.datastore[key].value = args
        else:
            try:
                self.datastore[key] = value
            except TypeError:
                raise AttributeError('The MetadataObject cannot accept the following key:value pair '
                                     'provided {key}:{value}'.format(key=key,
                                                                     value=value))

    def __getitem__(self, key):
        try:
            return self.datastore[key]
        except KeyError:
            raise AttributeError('The MetadataObject has not been initialised with the following key: {key}'
                                 .format(key=key))

    def dump(self):
        """Prints only the nested dictionary values; removes __methods__ and __members__ attributes"""
        metadata = dict()
        for attr in sorted(self.datastore):
            # Initialise the attribute (e.g. sample.general) key in the metadata dictionary
            metadata[attr] = dict()
            # Ignore attributes that begin with '__'
            if not attr.startswith('__'):
                # If self.datastore[attribute] is a primitive datatype, populate the metadata dictionary with
                # the attr: self.datastore[attr] pair
                # e.g. attr: name,  self.datastore[attr]: 2013-SEQ-0072
                if isinstance(self.datastore[attr], str) or \
                        isinstance(self.datastore[attr], list) or \
                        isinstance(self.datastore[attr], dict) or \
                        isinstance(self.datastore[attr], int):
                    metadata[attr] = self.datastore[attr]
                else:
                    # Otherwise, recursively convert GenObjects to nested dictionaries
                    metadata.update(self.nested_genobject(metadata, attr, self.datastore))
        return metadata

    def nested_genobject(self, metadata, attr, datastore):
        """
        Allow for the printing of nested GenObjects
        :param metadata: Nested dictionary containing the metadata. Will be further populated by this method
        :param attr: Current attribute being evaluated. Must be a GenObject e.g. sample.general
        :param datastore: The dictionary of the current attribute. Will be converted to nested dictionaries
        :return: Updated nested metadata dictionary with all GenObjects safely converted to dictionaries
        """
        # Iterate through all the key: value pairs of the current datastore[attr] datastore
        # e.g. reverse_reads <accessoryFunctions.accessoryFunctions.GenObject object at 0x7fe153b725f8>
        for key, value in sorted(datastore[attr].datastore.items()):
            # If the type(value) is a GenObject, then JSON serialization will not work
            if 'GenObject' in str(type(value)):
                # Initialise the nested attribute: key nested dictionary within the metadata dictionary
                # e.g. attr: 100_100, key: reverse_reads
                metadata[attr][key] = dict()
                # Iterate through the nested keys and nested values within the value datastore
                # e.g. nested_key: length, nested_value: 100
                for nested_key, nested_datastore in sorted(value.datastore.items()):
                    # Create an additional dictionary layer within the metadata dictionary
                    metadata[attr][key][nested_key] = dict()
                    # If the type(nested_datastore) is a GenObject, recursively run this method to update the
                    # metadata dictionary, supply the newly created nested dictionary: metadata[attr][key] as
                    # the input metadata dictionary, the nested key as the input attribute, and the datastore of
                    # value as the input datastore
                    # e.g. key: 100_100,
                    # datastore: <accessoryFunctions.accessoryFunctions.GenObject object at 0x7fc526001e80>
                    if 'GenObject' in str(type(nested_datastore)):
                        metadata[attr][key].update(
                            self.nested_genobject(metadata[attr][key], nested_key, value.datastore))
                    # If the nested datastore is not a GenObject, populate the nested metadata dictionary with
                    # the attribute, key, nested key, and nested datastore
                    # e.g. attr: 100_100, key: reverse_reads, nested_key: length, nested_datastore: 100
                    else:
                        metadata[attr][key][nested_key] = nested_datastore
            # Non-GenObjects can (usually) be added to the metadata dictionary without issues
            else:
                try:
                    if key not in self.unwanted_keys:
                        metadata[attr][key] = value
                except AttributeError:
                    print('dumperror', attr)
        # Return the metadata
        return metadata


class MakeBlastDB(AbstractCommandline):
    """Base makeblastdb wrapper"""
    def __init__(self, cmd='makeblastdb', **kwargs):
        assert cmd is not None
        extra_parameters = [
            # Core:
            _Switch(["-h", "h"],
                    "Print USAGE and DESCRIPTION;  ignore other arguments."),
            _Switch(["-help", "help"],
                    "Print USAGE, DESCRIPTION and ARGUMENTS description; "
                    "ignore other arguments."),
            _Switch(["-version", "version"],
                    "Print version number;  ignore other arguments."),
            # Output configuration options
            _Option(["-out", "out"],
                    "Output file prefix for db.",
                    filename=True,
                    equate=False),
            _Option(["-in", "db"],
                    "The sequence create db with.",
                    filename=True,
                    equate=False),  # Should this be required?
            _Option(["-dbtype", "dbtype"],
                    "Molecule type of target db (string, 'nucl' or 'prot').",
                    equate=False)]
        try:
            # Insert extra parameters - at the start just in case there
            # are any arguments which must come last:
            self.parameters = extra_parameters + self.parameters
        except AttributeError:
            # Should we raise an error?  The subclass should have set this up!
            self.parameters = extra_parameters
        AbstractCommandline.__init__(self, cmd, **kwargs)


def combinetargets(targets, targetpath, mol_type='nt'):
    """
    Creates a set of all unique sequences in a list of supplied FASTA files. Properly formats headers and sequences
    to be compatible with local pipelines. Splits hybrid entries. Removes illegal characters.
    :param targets: fasta gene targets to combine
    :param targetpath: folder containing the targets
    """
    make_path(targetpath)
    with open(os.path.join(targetpath, 'combinedtargets.fasta'), 'w') as combined:
        idset = set()
        for target in targets:
            # Remove non-unicode characters present in the FASTA files
            cleanedstring = str()
            # Read in the file as binary
            with open(target, 'rb') as fasta:
                # Import all the text
                text = fasta.read()
                # Convert the binary variable to a string, ignoring non-UTF-8 characters
                cleanedstring += text.decode('utf-8', 'ignore')
            # Overwrite the file with the clean string
            with open(target, 'w') as fasta:
                fasta.write(cleanedstring)
            # Clean up each record
            for record in SeqIO.parse(target, 'fasta'):
                # In case FASTA records have been spliced together, allow for the splitting of
                # these records
                if '>' in record.seq:
                    # Split the two records apart on '>' symbols
                    record.seq, hybrid = record.seq.split('>')
                    # Split the header from the sequence e.g. sspC:6:CP003808.1ATGGAAAGTACATTAGA...
                    # will be split into sspC:6:CP003808.1 and ATGGAAAGTACATTAGA
                    hybridid, seq = re.findall('(.+\d+\.\d)(.+)', str(hybrid))[0]
                    # Replace and dashes in the record.id with underscores
                    hybridid = hybridid.replace('-', '_')
                    # Convert the string to a seq object
                    if mol_type == 'nt':
                        hybridseq = Seq(seq, generic_dna)
                    else:
                        hybridseq = Seq(seq, generic_protein)
                    # Create a SeqRecord of the sequence - use the sequence object and id
                    hybridrecord = SeqRecord(hybridseq,
                                             description='',
                                             id=hybridid)

                    # Remove and dashes or 'N's from the sequence data - makeblastdb can't handle sequences
                    # with gaps
                    # noinspection PyProtectedMember
                    hybridrecord.seq._data = hybridrecord.seq._data.replace('-', '').replace('N', '')
                    # Write the original record to the file
                    # Extract the sequence record from each entry in the multifasta
                    # Replace and dashes in the record.id with underscores
                    record.id = record.id.replace('-', '_')
                    # Remove and dashes or 'N's from the sequence data - makeblastdb can't handle sequences
                    # with gaps
                    # noinspection PyProtectedMember
                    record.seq._data = record.seq._data.replace('-', '').replace('N', '')
                    # Clear the name and description attributes of the record
                    record.name = ''
                    record.description = ''
                    if record.id not in idset:
                        SeqIO.write(record, combined, 'fasta')
                    if hybridrecord.id not in idset:
                        # Write the second record to file
                        SeqIO.write(hybridrecord, combined, 'fasta')
                        idset.add(hybridrecord.id)

                else:
                    # Extract the sequence record from each entry in the multifasta
                    # Replace and dashes in the record.id with underscores
                    record.id = record.id.replace('-', '_')
                    # Remove and dashes or 'N's from the sequence data - makeblastdb can't handle sequences
                    # with gaps
                    # noinspection PyProtectedMember
                    record.seq._data = record.seq._data.replace('-', '').replace('N', '')
                    # Clear the name and description attributes of the record
                    record.name = ''
                    record.description = ''
                    if record.id not in idset:
                        SeqIO.write(record, combined, 'fasta')
                        idset.add(record.id)


class KeyboardInterruptError(Exception):
    pass
