#!/usr/bin/env python3
from olctools.accessoryFunctions.accessoryFunctions import combinetargets, MetadataObject, make_path, \
    run_subprocess, SetupLogging
from olctools.databasesetup import enterobase_api_download, get_mlst, get_rmlst
from datetime import datetime
from argparse import ArgumentParser
from subprocess import call
import urllib.request
from glob import glob
import logging
import tarfile
import zipfile
import shutil
import click
import gzip
import ssl
import os
__author__ = 'adamkoziol'


class DatabaseSetup(object):

    def cowbat(self):
        """
        Run all the methods
        """
        logging.info('Beginning COWBAT database downloads')
        if self.overwrite or not os.path.isdir(os.path.join(self.databasepath, 'genesippr')):
            self.sipprverse_targets(databasepath=self.databasepath)
        if self.overwrite or not os.path.isdir(os.path.join(self.databasepath, 'coregenome')):
            self.cowbat_targets(databasepath=self.databasepath)
        if self.overwrite or not os.path.isdir(os.path.join(self.databasepath, 'ConFindr')):
            self.confindr_targets()
        if self.overwrite or not os.path.isdir(os.path.join(self.databasepath, 'mash')):
            self.mash(databasepath=self.databasepath)
        if self.overwrite or not os.path.isdir(os.path.join(self.databasepath, 'MLST')):
            if self.enterobase:
                self.enterobase_mlst(databasepath=self.databasepath)
                self.mlst(databasepath=self.databasepath)
            else:
                self.mlst(databasepath=self.databasepath,
                          genera=('Bacillus', 'Campylobacter', 'Cronobacter', 'Escherichia', 'Listeria',
                                  'Salmonella', 'Staphylococcus', 'Vibrio', 'Yersinia'))
        if self.overwrite or not os.path.isdir(os.path.join(self.databasepath, 'rMLST')):
            self.rmlst(databasepath=self.databasepath,
                       credentials=self.credentials)
        if self.overwrite or not os.path.isdir(os.path.join(self.databasepath, 'cgMLST')):
            self.enterobase_cgmlst(databasepath=self.databasepath)
        if self.overwrite or not os.path.isdir(os.path.join(self.databasepath, 'univec')):
            self.univec(databasepath=self.databasepath)
        if self.overwrite or not os.path.isdir(os.path.join(self.databasepath, 'resfinder')):
            self.cge_db_downloader(databasepath=self.databasepath,
                                   analysistype='resfinder',
                                   dbname='resfinder_db')
        if self.overwrite or not os.path.isdir(os.path.join(self.databasepath, 'virulence')):
            self.cge_db_downloader(databasepath=self.databasepath,
                                   analysistype='virulence',
                                   dbname='virulencefinder_db')
        if self.overwrite or not os.path.isdir(os.path.join(self.databasepath, 'serosippr')):
            self.cge_db_downloader(databasepath=self.databasepath,
                                   analysistype='serosippr',
                                   dbname='serotypefinder_db')
        if self.overwrite or not os.path.isdir(os.path.join(self.databasepath, 'pointfinder')):
            self.cge_db_downloader(databasepath=self.databasepath,
                                   analysistype='pointfinder',
                                   dbname='pointfinder_db')
        if self.overwrite or not os.path.isdir(os.path.join(self.databasepath, 'clark')):
            self.clark(databasepath=self.databasepath)
        self.download_date()

    def sipprverse_full(self):
        """
        Run a subset of the methods - only the targets used in the sipprverse are required here
        """
        logging.info('Beginning sipprverse full database downloads')
        if self.overwrite or not os.path.isdir(os.path.join(self.databasepath, 'genesippr')):
            self.sipprverse_targets(databasepath=self.databasepath)
        if self.overwrite or not os.path.isdir(os.path.join(self.databasepath, 'ConFindr')):
            self.confindr_targets()
        if self.overwrite or not os.path.isdir(os.path.join(self.databasepath, 'mash')):
            self.mash(databasepath=self.databasepath)
        if self.overwrite or not os.path.isdir(os.path.join(self.databasepath, 'MLST')):
            self.mlst(databasepath=self.databasepath)
        if self.overwrite or not os.path.isdir(os.path.join(self.databasepath, 'rMLST')):
            self.rmlst(databasepath=self.databasepath,
                       credentials=self.credentials)
        if self.overwrite or not os.path.isdir(os.path.join(self.databasepath, 'resfinder')):
            self.cge_db_downloader(databasepath=self.databasepath,
                                   analysistype='resfinder',
                                   dbname='resfinder_db')
        if self.overwrite or not os.path.isdir(os.path.join(self.databasepath, 'virulence')):
            self.cge_db_downloader(databasepath=self.databasepath,
                                   analysistype='virulence',
                                   dbname='virulencefinder_db')
        if self.overwrite or not os.path.isdir(os.path.join(self.databasepath, 'serosippr')):
            self.cge_db_downloader(databasepath=self.databasepath,
                                   analysistype='serosippr',
                                   dbname='serotypefinder_db')

    def sipprverse_method(self):
        """
        Reduced subset again. Only sipprverse, MASH, and confindr targets are required
        """
        logging.info('Beginning sipprverse method database downloads')
        if self.overwrite or not os.path.isdir(os.path.join(self.databasepath, 'genesippr')):
            self.sipprverse_targets(databasepath=self.databasepath)
        if self.overwrite or not os.path.isdir(os.path.join(self.databasepath, 'ConFindr')):
            self.confindr_targets()
        if self.overwrite or not os.path.isdir(os.path.join(self.databasepath, 'mash')):
            self.mash(databasepath=self.databasepath)

    def rmlst_method(self):
        """
        Run only the rMLST download
        """
        if self.overwrite or not os.path.isdir(os.path.join(self.databasepath, 'rMLST')):
            self.rmlst(databasepath=self.databasepath,
                       credentials=self.credentials)

    def sipprverse_targets(self, databasepath, database_name='sipprverse', download_id='18130808'):
        """
        Download OLC-specific sipprverse targets
        :param databasepath: path to use to save the database
        :param database_name: name of current database
        :param download_id: figshare identifier of .tar.gz file
        """
        self.custom_databases(databasepath=databasepath,
                              database_name=database_name,
                              download_id=download_id)

    def cowbat_targets(self, databasepath, database_name='COWBAT', download_id='25319129'):
        """
        Download OLC-specific COWBAT targets
        :param databasepath: path to use to save the database
        :param database_name: name of current database
        :param download_id: figshare identifier of .tar.gz file
        """
        self.custom_databases(databasepath=databasepath,
                              database_name=database_name,
                              download_id=download_id)

    def confindr_targets(self, database_name='ConFindr'):
        """
        Download OLC-specific ConFindr targets
        :param database_name: name of current database
        """
        logging.info('Downloading ConFindr databases.')
        # NOTE: Need ConFindr >= 0.5.0 for this to work.
        secret_file = os.path.join(self.credentials, 'secret.txt')
        condfindr_download = 'confindr_database_setup -s {secret} -o {output}'\
            .format(secret=secret_file,
                    output=os.path.join(self.databasepath, database_name))
        # Use subprocess.call rather than run_subprocess, as there is a prompt for input from the user
        call(condfindr_download, shell=True)

    @staticmethod
    def mlst(databasepath, genera=('Escherichia', 'Salmonella', 'Yersinia',
                                   'Campylobacter', 'Cronobacter', 'Listeria',
                                   'Bacillus', 'Staphylococcus', 'Vibrio')):
        """
        Download the necessary up-to-date MLST profiles and alleles from pubmlst
        :param databasepath: path to use to save the database
        :param genera: default genera for which alleles and profiles should be downloaded
        """
        logging.info('Downloading MLST databases from PubMLST')
        for genus in genera:
            # Create an object to pass to the get_mlst script
            args = MetadataObject()
            # Populate the object with the necessary attributes
            args.genus = genus
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

    def enterobase_mlst(self, databasepath, genera=('Escherichia', 'Salmonella', 'Yersinia')):
        """
        Download the necessary up-to-date MLST profiles and alleles from Enterobase
        :param databasepath: path to use to save the database
        :param genera: default genera for which alleles and profiles should be downloaded
        """
        logging.info('Downloading MLST databases from Enterobase')
        for genus in genera:
            self.enterobase_download(scheme='MLST_Achtman',
                                     databasepath=databasepath,
                                     analysis='MLST',
                                     genus=genus)

    def enterobase_download(self, scheme, databasepath, analysis, genus):
        """
        Download the appropriate scheme (MLST_Achtman/cgMLST) for the requested genus
        :param scheme: Scheme to download
        :param databasepath: path to use to save the database
        :param analysis: MLST or cgMLST
        :param genus: Genus for which alleles and profiles should be downloaded
        """
        enterobase_genus = self.genus_dict[genus]
        enterobase_cmd = 'python -m olctools.databasesetup.enterobase_api_download -o {eg} -s {scheme} -d {db}' \
            .format(eg=enterobase_genus,
                    scheme=scheme,
                    db=databasepath)
        if not os.path.isdir(os.path.join(databasepath, analysis, genus)):
            #call(enterobase_cmd, shell=True)
            out, err = run_subprocess(enterobase_cmd)
            print(out, err)

    @staticmethod
    def rmlst(databasepath, credentials):
        """
        Get the most up-to-date profiles and alleles from pubmlst. Note that you will need the necessary access token
        and secret for this to work
        :param databasepath: path to use to save the database
        :param credentials: path to folder containing accessory token and secret.txt files
        """
        logging.info('Downloading rMLST database')
        # Set the name of the file to be used to determine if the database download and set-up was successful
        completefile = os.path.join(databasepath, 'rMLST', 'complete')
        if not os.path.isfile(completefile):
            # Create an object to send to the rMLST download script
            args = MetadataObject()
            # Add the path and start time attributes
            args.path = databasepath
            args.logging = logging
            args.credentials = credentials
            # Run the rMLST download
            get_rmlst.Get(args)
            # Create and populate the complete.txt file
            with open(completefile, 'w') as complete:
                complete.write('\n'.join(glob(os.path.join(databasepath, 'rMLST', '*'))))

    def enterobase_cgmlst(self, databasepath, genera=('Escherichia', 'Yersinia')):
        """
        Download the necessary up-to-date cgMLST profiles and alleles from Enterobase
        :param databasepath: path to use to save the database
        :param genera: default genera for which alleles and profiles should be downloaded
        """
        logging.info('Downloading MLST databases from Enterobase')
        for genus in genera:
            self.enterobase_download(scheme='cgMLST',
                                     databasepath=databasepath,
                                     analysis='cgMLST',
                                     genus=genus)

    def clark(self, databasepath):
        """
        Download and set-up the CLARK database using the set_targets.sh script. Use defaults of bacteria for database
        type, and species for taxonomic level
        :param databasepath: path to use to save the database
        """
        # Determine the location of the CLARK scripts
        self.clarkpath = os.path.dirname(shutil.which('CLARK'))
        if self.clarkpath is not None:
            self.clarkpath = os.path.join(self.clarkpath)
        if self.clarkpath:
            logging.info('Downloading CLARK database')
            # Create the folder in which the database is to be stored
            databasepath = self.create_database_folder(databasepath, 'clark')
            # Set the call to create the database - use the --light option, as we don't require the full database
            targetcall = 'cd {clarkpath} && ../opt/clark/set_targets.sh {dbpath} bacteria --species --light'\
                .format(clarkpath=self.clarkpath,
                        dbpath=databasepath)
            # Download the database
            self.database_clone(targetcall, databasepath)
        else:
            logging.warning('No CLARK scripts detected in $PATH. Cannot download database.')

    def mash(self, databasepath):
        """
        Download the pre-computed sketch of the RefSeq database, and compress it with gzip
        :param databasepath: path to use to save the database
        """
        logging.info('Downloading pre-computed RefSeq MASH sketches')
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
                                   target_url='https://obj.umiacs.umd.edu/marbl_publications/mash/refseq.genomes.k21s1000.msh',
                                   complete=True)

    def univec(self, databasepath):
        """
        Download the UniVec core database
        :param databasepath: path to use to save the database
        """
        logging.info('Downloading univec database')
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
    def url_request(target_url, output_file, context):
        """
        Use urllib to download the requested file from the target URL. Use the click progress bar to print download
        progress
        :param target_url: URL from which the file is to be downloaded
        :param output_file: Name and path of local copy of file
        :param context: ssl.create_default_context() with check_hostname set to False, and verify_mode 
        set to ssl.CERT_NONE
        """
        # Create the request
        request = urllib.request.urlopen(target_url, context=context)
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

    def custom_databases(self, databasepath, database_name, download_id, f_type='files', post_id=None,
                         compression='tar', nested=False, complete=False):
        """
        Download and extract a .tar.gz file from figshare
        :param databasepath: Name and path of where the database files are to be downloaded
        :param database_name: Name of the database e.g. sipprverse
        :param download_id: Figshare ID of the targets file
        :param f_type: STR MOB-suite databases have the 'articles' keyword in the figshare URL, while OLC databases
        all have the 'files' keyword
        :param post_id: STR MOB-suite databases have 'versions/1' appended at the end of the figshare URL.
        :param compression: STR MOB-suite databases are .zip files, while OLC databases are .tar.gz
        :param nested: Boolean of whether the targets file has nested folders that must be accounted for
        :param complete: Boolean of whether the completefile should be completed
        """
        logging.info('Downloading {} databases'.format(database_name))
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
        target_url = 'https://ndownloader.figshare.com/{type}/{id}'.format(type=f_type,
                                                                           id=download_id)
        if post_id:
            target_url += '/{post}'.format(post=post_id)
        logging.debug(target_url)
        if not os.path.isfile(completefile):
            self.url_request(target_url=target_url,
                             output_file=tar_file,
                             context=self.context)
        # Decompress the file
        self.decompress(databasepath=databasepath,
                        database_name=database_name,
                        compression=compression,
                        compressed_file=tar_file)
        # Create the completefile
        if complete:
            with open(completefile, 'w') as complete:
                complete.write('')

    @staticmethod
    def decompress(databasepath, database_name, compression, compressed_file):
        """
        Decompress the provided file using the appropriate library
        :param databasepath: Name and path of where the database files are to be downloaded
        :param database_name: Name of the database e.g. sipprverse
        :param compression: STR MOB-suite databases are .zip files, while OLC databases are .tar.gz
        :param compressed_file: Compressed file to process
        """
        # Extract the databases from the archives
        if os.path.isfile(compressed_file):
            if compression == 'tar':
                logging.info('Extracting {dbname} from archives'.format(dbname=database_name))
                with tarfile.open(compressed_file, 'r') as tar:
                    # Decompress the archive
                    tar.extractall(path=databasepath)
            elif compression == 'gz':
                with gzip.open(compressed_file, 'rb') as gz:
                    file_name = os.path.basename(os.path.splitext(compressed_file)[0])
                    output_file = os.path.join(databasepath,
                                               database_name,
                                               file_name)
                    logging.info('Extracting {file_name} from archives'.format(file_name=file_name))
                    with open(output_file, 'wb') as output:
                        shutil.copyfileobj(gz, output)
            else:
                logging.info('Extracting {dbname} from archives'.format(dbname=database_name))
                with zipfile.ZipFile(compressed_file, 'r') as zip_file:
                    zip_file.extractall(path=databasepath)
            # Delete the archive file
            os.remove(compressed_file)

    def cge_db_downloader(self,  databasepath, analysistype, dbname, extension_in='fsa', extension_out='tfa'):
        """
        Clones CGE databases into appropriate folder. Creates properly formatted file with non-redundant sequences
        :param databasepath: path to use to save the database
        :param analysistype: The name of the database folder to create
        :param dbname: The name of the database repository on bitbucket
        :param extension_in: The file extension of the FASTA files in the database
        :param extension_out: The desired extension for the FASTA files
        """
        logging.info('Downloading {} database'.format(analysistype))
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

    @staticmethod
    def create_database_folder(databasepath, database):
        """
        Create an appropriately named folder in which the database is to be stored
        :param databasepath: path to use to save the database
        :param database: the name of the database folder to create
        :return: the absolute path of the folder
        """
        logging.info('Setting up {} database'.format(database))
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
                             output_file=output_file,
                             context=self.context)
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

    def download_date(self):
        """
        Write the current date to file.
        """
        with open(os.path.join(self.databasepath, 'download_date'), 'w') as download:
            download.write('{:%Y-%m-%d}'.format(datetime.today()))

    def __init__(self, databasepath=None, debug=False, credentials=None, overwrite=False, enterobase=False):
        # Initialise the custom logging handler
        SetupLogging(debug)
        # Create class variables from arguments
        if databasepath.startswith('~'):
            self.databasepath = os.path.abspath(os.path.expanduser(os.path.join(databasepath)))
        else:
            self.databasepath = os.path.abspath(os.path.join(databasepath))
        make_path(self.databasepath)
        assert os.path.isdir(self.databasepath)
        if credentials:
            if credentials.startswith('~'):
                self.credentials = os.path.abspath(os.path.expanduser(os.path.join(credentials)))
            else:
                self.credentials = os.path.abspath(os.path.join(credentials))
        self.overwrite = overwrite
        assert type(self.overwrite) is bool, 'Overwrite variable must be a Boolean. You provided "{var}" with ' \
                                             'type {type}'.format(var=self.overwrite,
                                                                  type=type(self.overwrite))
        self.clarkpath = str()
        # Enterobase
        self.enterobase = enterobase
        self.genus_dict = {
            'Escherichia': 'ecoli',
            'Salmonella': 'senterica',
            'Yersinia': 'yersinia'
        }
        # Create a context to all for disabling SSL verifying
        self.context = ssl.create_default_context()
        self.context.check_hostname = False
        self.context.verify_mode = ssl.CERT_NONE


# If the script is called from the command line, then call the argument parser
if __name__ == '__main__':
    # Parser for arguments
    parser = ArgumentParser(description='Downloads and sets up required databases')
    parser.add_argument('-d', '--databasepath',
                        required=True,
                        help='Absolute path to location to store database files. Include any version numbers if '
                             'required.')
    parser.add_argument('-c,', '--credentials',
                        help='Name and path of folder containing required rMLST credentials.')
    parser.add_argument('-o', '--overwrite',
                        action='store_true',
                        help='Optionally allow for the overwriting of database files in the databasepath. Defaults to '
                             'False, so if the output folder exists, that part of the download will be skipped.')
    parser.add_argument('-s', '--sipprverse_full',
                        action='store_true',
                        help='Optionally only download the databases used in the sipprverse. These include: genesippr, '
                             'GDCS, sixteenS, ConFindr, MASH, MLST, rMLST, ResFindr, VirulenceFinder, '
                             'and SerotypeFinder')
    parser.add_argument('-m', '--sipprverse_method',
                        action='store_true',
                        help='Optionally only download the databases used by the sipprverse method: genesippr, '
                             'sixteenS, GDCS, MASH, and ConFindr')
    parser.add_argument('-e', '--enterobase',
                        action='store_true',
                        help='Use Enterobase to download MLST definitions for Escherichia, Salmonella, and Yersinia, '
                             'as well as cgMLST schemes for Escherichia and Yersinia. Disabled by default')
    parser.add_argument('-v', '--verbose',
                        action='store_true',
                        help='Option to include debug level logging messages. Default is false')
    parser.add_argument('-r', '--rmlst',
                        action='store_true',
                        help='Optionally only download the rMLST database')
    parser.add_argument('-res', '--resfinder',
                        action='store_true',
                        help='Only download the Resfinder database')
    # Get the arguments into an object
    arguments = parser.parse_args()
    # Create an object
    pipeline = DatabaseSetup(databasepath=arguments.databasepath,
                             debug=arguments.verbose,
                             credentials=arguments.credentials,
                             overwrite=arguments.overwrite,
                             enterobase=arguments.enterobase)
    # Run the appropriate analyses
    if arguments.resfinder:
        pipeline.cge_db_downloader(
            databasepath=pipeline.databasepath,
            analysistype='resfinder',
            dbname='resfinder_db'
        )
        raise SystemExit
    if not arguments.credentials:
        logging.error('Please provide the name and path of the folder containing your rMLST credentials with the -c '
                      'argument.')
        raise SystemExit
    if arguments.sipprverse_full:
        pipeline.sipprverse_full()
    elif arguments.sipprverse_method:
        pipeline.sipprverse_method()
    elif arguments.rmlst:
        pipeline.rmlst_method()
    else:
        pipeline.cowbat()
