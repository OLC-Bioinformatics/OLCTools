#!/usr/bin/env python 3
from accessoryFunctions.accessoryFunctions import combinetargets, MetadataObject, make_path, \
    run_subprocess, setup_logging
from databasesetup import get_mlst, get_rmlst
from argparse import ArgumentParser
import urllib.request
from time import time
from glob import glob
import tarfile
import shutil
import click
import os
__author__ = 'adamkoziol'


class DatabaseSetup(object):

    def main(self):
        """
        Run the methods
        """
        self.custom_databases(databasepath=self.databasepath,
                              database_name='sipprverse',
                              download_id='13143830')
        self.custom_databases(databasepath=self.databasepath,
                              database_name='COWBAT',
                              download_id='13152896')
        self.custom_databases(databasepath=self.databasepath,
                              database_name='ConFindr',
                              download_id='11864267',
                              nested=True)
        self.custom_databases(databasepath=self.databasepath,
                              database_name='plasmidextractor',
                              download_id='9827323',
                              nested=True,
                              complete=True)
        self.mash(databasepath=self.databasepath)
        self.mlst(databasepath=self.databasepath)
        self.rmlst(databasepath=self.databasepath)
        self.univec(databasepath=self.databasepath)
        self.cge_db_downloader(databasepath=self.databasepath,
                               analysistype='resfinder',
                               dbname='resfinder_db')
        self.cge_db_downloader(databasepath=self.databasepath,
                               analysistype='virulence',
                               dbname='virulencefinder_db')
        self.cge_db_downloader(databasepath=self.databasepath,
                               analysistype='serosippr',
                               dbname='serotypefinder_db')
        self.cge_db_downloader(databasepath=self.databasepath,
                               analysistype='plasmidfinder',
                               dbname='plasmidfinder_db')
        self.cge_db_downloader(databasepath=self.databasepath,
                               analysistype='virulence',
                               dbname='virulencefinder_db')
        self.clark(databasepath=self.databasepath)

    def mlst(self, databasepath, genera=('Escherichia', 'Vibrio', 'Campylobacter', 'Listeria',
                                         'Bacillus', 'Staphylococcus', 'Salmonella')):
        """
        Download the necessary up-to-date MLST profiles and alleles
        :param databasepath: path to use to save the database
        :param genera: default genera for which alleles and profiles should be downloaded
        """
        self.logging.info('Downloading MLST databases')
        for genus in genera:
            # Create an object to pass to the get_mlst script
            args = MetadataObject()
            # Populate the object with the necessary attributes
            args.species = genus
            args.repository_url = 'http://pubmlst.org/data/dbases.xml'
            args.force_scheme_name = False
            args.path = os.path.join(databasepath, 'MLST', genus)
            # Create the name of the file to be used to determine if the database download and setup was successful
            completefile = os.path.join(args.path, 'complete')
            # Only download the files if the download was not previously successful
            if not os.path.isfile(completefile):
                # Run the download
                get_mlst.main(args)
                # Create and populate the complete.txt file
                with open(completefile, 'w') as complete:
                    complete.write('\n'.join(glob(os.path.join(args.path, '*'))))

    def rmlst(self, databasepath):
        """
        Get the most up-to-date profiles and alleles from pubmlst. Note that you will need the necessary access token
        and secret for this to work
        :param databasepath: path to use to save the database
        """
        self.logging.info('Downloading rMLST database')
        # Set the name of the file to be used to determine if the database download and set-up was successful
        completefile = os.path.join(self.databasepath, 'rMLST', 'complete')
        if not os.path.isfile(completefile):
            # Create an object to send to the rMLST download script
            args = MetadataObject()
            # Add the path and start time attributes
            args.path = databasepath
            args.logging = self.logging
            # Run the rMLST download
            get_rmlst.Get(args)
            # Create and populate the complete.txt file
            with open(completefile, 'w') as complete:
                complete.write('\n'.join(glob(os.path.join(databasepath, 'rMLST', '*'))))

    def clark(self, databasepath):
        """
        Download and set-up the CLARK database using the set_targets.sh script. Use defaults of bacteria for database
        type, and species for taxonomic level
        :param databasepath: path to use to save the database
        """
        if self.clarkpath:
            self.logging.info('Downloading CLARK database')
            # Create the folder in which the database is to be stored
            databasepath = self.create_database_folder(databasepath, 'clark')
            # Set the call to create the database - use the --light option, as we don't require the full database
            targetcall = 'cd {clarkpath} && ../opt/clark/set_targets.sh {dbpath} bacteria --species --light'\
                .format(clarkpath=self.clarkpath,
                        dbpath=databasepath)
            # Download the database
            self.database_clone(targetcall, databasepath)
        else:
            self.logging.warning('No CLARK scripts detected in $PATH. Cannot download database.')

    def mash(self, databasepath):
        """
        Download the pre-computed sketch of the RefSeq database, and compress it with gzip
        :param databasepath: path to use to save the database
        """
        self.logging.info('Downloading pre-computed RefSeq MASH sketches')
        # Create the folder in which the database is to be stored
        databasepath = self.create_database_folder(databasepath=databasepath,
                                                   database='mash')
        output_file = os.path.join(databasepath, 'assembly_summary_refseq.txt')
        # Download the assembly summary RefSeq document
        if not os.path.isfile(output_file):
            self.database_download(output_file=output_file,
                                   database_path=databasepath,
                                   target_url='ftp://ftp.ncbi.nih.gov/genomes/ASSEMBLY_REPORTS/'
                                              'assembly_summary_refseq.txt'
                                   )
        # Set the call to create the database
        output_file = os.path.join(databasepath, 'RefSeqSketchesDefaults.msh')
        # Download the database
        if not os.path.isfile(output_file):
            self.database_download(output_file=output_file,
                                   database_path=databasepath,
                                   target_url='https://gembox.cbcb.umd.edu/mash/refseq.genomes.k21s1000.msh',
                                   complete=True)

    def univec(self, databasepath):
        """
        Download the UniVec core database
        :param databasepath: path to use to save the database
        """
        self.logging.info('Downloading univec database')
        databasepath = self.create_database_folder(databasepath, 'univec')
        # Set the name of the output file
        outputfile = os.path.join(databasepath, 'UniVec_core.tfa')
        target_url = 'ftp://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec_Core'
        self.database_download(output_file=outputfile,
                               target_url=target_url,
                               database_path=databasepath)
        # Create a copy of the file with a .fasta extension
        if os.path.isfile(outputfile):
            renamed = os.path.splitext(outputfile)[0] + '.fasta'
            shutil.copy(outputfile, renamed)

    @staticmethod
    def url_request(target_url, output_file):
        """
        Use urllib to download the requested file from the target URL. Use the click progress bar to print download
        progress
        :param target_url: URL from which the file is to be downloaded
        :param output_file: Name and path of local copy of file
        """
        # Create the request
        request = urllib.request.urlopen(target_url)
        # Open the destination file to write
        with open(output_file, 'wb') as targets:
            # Calculate the total file size - will be used by the progress bar
            total_length = int(request.headers.get('content-length'))
            # Create a click progress bar using the total length calculated above
            with click.progressbar(length=total_length,
                                   label='Downloading files') as bar:
                while True:
                    # Break up the download into chunks of 4096 bytes
                    data = request.read(4096)
                    # Break the loop when the download finishes/errors
                    if not data:
                        break
                    # Write the chunk to file
                    targets.write(data)
                    # Update the progress bar
                    bar.update(len(data))

    def custom_databases(self, databasepath, database_name, download_id, nested=False, complete=False):
        """
        Download and extract a .tar.gz file from figshare
        :param databasepath: Name and path of where the database files are to be downloaded
        :param database_name: Name of the database e.g. sipprverse
        :param download_id: Figshare ID of the targets file
        :param nested: Boolean of whether the targets file has nested folders that must be accounted for
        :param complete: Boolean of whether the completefile should be completed
        """
        self.logging.info('Downloading {} databases'.format(database_name))
        # ConFindr has a nested 'databases' folder
        if nested:
            databasepath = os.path.join(databasepath, database_name)
        # Set the name and path of the file that is created when the download is successful
        completefile = os.path.join(databasepath, 'complete')
        # Create the database folder if necessary
        make_path(databasepath)
        # Set the name of the targets file
        tar_file = os.path.join(databasepath, download_id)
        # Create the target download call
        target_url = 'https://ndownloader.figshare.com/files/{id}'.format(id=download_id)
        if not os.path.isfile(completefile):
            self.url_request(target_url=target_url,
                             output_file=tar_file)
        # Extract the databases from the archives
        self.logging.info('Extracting {dbname} database from archives'.format(dbname=database_name))
        if os.path.isfile(tar_file):
            with tarfile.open(tar_file, 'r') as tar:
                # Decompress the archive
                tar.extractall(path=databasepath)
            # Delete the archive file
            os.remove(tar_file)
        # Create the completefile
        if complete:
            with open(completefile, 'w') as complete:
                complete.write('')

    def cge_db_downloader(self,  databasepath, analysistype, dbname, extension_in='fsa', extension_out='tfa'):
        """
        Clones CGE databases into appropriate folder. Creates properly formatted file with non-redundant sequences
        :param databasepath: path to use to save the database
        :param analysistype: The name of the database folder to create
        :param dbname: The name of the database repository on bitbucket
        :param extension_in: The file extension of the FASTA files in the database
        :param extension_out: The desired extension for the FASTA files
        """
        self.logging.info('Downloading {} database'.format(analysistype))
        if analysistype == 'serosippr':
            databasepath = os.path.join(databasepath, analysistype, 'Escherichia')
        else:
            databasepath = os.path.join(databasepath, analysistype)
        targetcall = 'git clone https://bitbucket.org/genomicepidemiology/{db}.git {atype}'\
            .format(db=dbname,
                    atype=databasepath)
        # Download the database
        self.database_clone(targetcall, databasepath)
        # Create a variable to use in creating the combined targets file
        extension = extension_in
        # If the extension_out is different than extension_in, rename the files to have the appropriate extension
        if extension_in != extension_out:
            # Create a list of all the FASTA files with the input extension
            fastafiles = glob(os.path.join(databasepath, '*.{ext}'.format(ext=extension_in)))
            for fasta in fastafiles:
                # Split the extension
                filename = os.path.splitext(fasta)[0]
                # Rename the files
                os.rename(fasta, '{fn}.{ex}'.format(fn=filename,
                                                    ex=extension_out))
            # Update the variable to use when creating the combined targets file
            extension = extension_out
        # Create the combined targets file to use in the OLC pipelines
        if not os.path.isfile(os.path.join(databasepath, 'combinedtargets.fasta')):
            # Create the combinedtargets.fasta file - this will combine all the FASTA files in the downloaded database
            # into a properly-formatted, non-redundant FASTA database
            databasefiles = glob(os.path.join(databasepath, '*.{ext}'.format(ext=extension)))
            combinetargets(databasefiles, databasepath)

    def create_database_folder(self, databasepath, database):
        """
        Create an appropriately named folder in which the database is to be stored
        :param databasepath: path to use to save the database
        :param database: the name of the database folder to create
        :return: the absolute path of the folder
        """
        self.logging.info('Setting up {} database'.format(database))
        # Define the path to store the database files
        databasepath = os.path.join(databasepath, database)
        # Create the path as required
        make_path(databasepath)
        return databasepath

    def database_download(self, output_file, target_url, database_path, complete=False):
        """
        Check to see if the download has previously been completed. Run the download if necessary. Create the
        completefile if required
        :param output_file: Name and path of local copy of downloaded target
        :param target_url: URL of the target to download
        :param database_path: Path on the local filesystem in which the file is to be downloaded
        :param complete: Boolean to determine whether a completefile should be created
        """
        # Create a file to store the logs; it will be used to determine if the database was downloaded and set-up
        completefile = os.path.join(database_path, 'complete')
        if not os.path.isfile(completefile):
            self.url_request(target_url=target_url,
                             output_file=output_file)
            if complete:
                # Create the completefile
                with open(completefile, 'w') as complete:
                    complete.write('')

    @staticmethod
    def database_clone(targetcall, databasepath, complete=False):
        """
        Checks to see if the database has already been downloaded. If not, runs the system call to
        download the database, and writes stdout and stderr to the logfile
        :param targetcall: system call to download, and possibly set-up the database
        :param databasepath: absolute path of the database
        :param complete: boolean variable to determine whether the complete file should be created
        """
        # Create a file to store the logs; it will be used to determine if the database was downloaded and set-up
        completefile = os.path.join(databasepath, 'complete')
        # Run the system call if the database is not already downloaded
        if not os.path.isfile(completefile):
            out, err = run_subprocess(targetcall)
            if complete:
                # Create the database completeness assessment file and populate it with the out and err streams
                with open(completefile, 'w') as complete:
                    complete.write(out)
                    complete.write(err)

    def __init__(self, args=None, debug=False):
        # Allow args to be absent
        if args:
            self.databasepath = os.path.join(args.databasepath)
            make_path(self.databasepath)
        # Create the logging object
        try:
            self.logging = setup_logging(debug=args.debug)
        except AttributeError:
            self.logging = setup_logging(debug=debug)
        try:
            # Determine the location of the CLARK scripts
            self.clarkpath = os.path.dirname(shutil.which('CLARK'))
        except AttributeError:
            self.clarkpath = None


# If the script is called from the command line, then call the argument parser
if __name__ == '__main__':
    # Parser for arguments
    parser = ArgumentParser(description='Downloads and sets up required databases')
    parser.add_argument('-d', '--databasepath',
                        required=True,
                        help='Absolute path to location to store database files. Include any version numbers if '
                             'required.')
    parser.add_argument('--debug',
                        default=False,
                        action='store_true',
                        help='Option to include debug level logging messages. Default is false')
    # Get the arguments into an object
    arguments = parser.parse_args()
    arguments.start = time()
    # Run the pipeline
    pipeline = DatabaseSetup(arguments)
    pipeline.main()
