#!/usr/bin/env python3

"""
Sets up the required databases for the COWBAT pipeline
"""

# Standard imports
from argparse import ArgumentParser
from datetime import datetime
from glob import glob
import gzip
import logging
import os
import shutil
import ssl
from subprocess import call
import tarfile
from typing import (
    Any,
    Optional
)
import urllib.request
import zipfile

# Third-party imports
from olctools.accessoryFunctions.accessoryFunctions import (
    combinetargets,
    run_subprocess,
)
from olctools.accessoryFunctions.metadata import CustomBox
from olctools.database_setup import get_mlst, get_rmlst

__author__ = 'adamkoziol'


class DatabaseSetup:
    """
    Run the necessary database setup methods
    """

    def __init__(
        self,
        database_path=None,
        credentials=None,
        overwrite=False,
        enterobase=False
    ):
        # Create class variables from arguments
        self.database_path = self._expand_path(database_path)

        # Create the database_path if required
        os.makedirs(self.database_path, exist_ok=True)

        # Ensure that the database path was created successfully
        assert os.path.isdir(self.database_path)

        # Set the path of the credentials file
        self.credentials = self._expand_path(
            credentials
        ) if credentials else None

        # Determine whether files should be overwritten
        self.overwrite = overwrite
        assert isinstance(
            self.overwrite, bool), (
                f'Overwrite variable must be a Boolean. You provided '
                f'"{self.overwrite}" ' f'with type {type(self.overwrite)}'
            )

        # Determine whether Enterobase is to be used for the downloads
        self.enterobase = enterobase

        # Set the MLST database names for desired genera
        self.genus_dict = {
            'Escherichia': 'ecoli',
            'Salmonella': 'senterica',
            'Yersinia': 'yersinia'
        }

        # Create a context to allow for disabling SSL verifying
        self.context = ssl.create_default_context()
        self.context.check_hostname = False
        self.context.verify_mode = ssl.CERT_NONE

    @staticmethod
    def _expand_path(path):
        """
        Expand user and relative paths to absolute paths.
        """
        return os.path.abspath(os.path.expanduser(path))

    def cowbat(self):
        """
        Run all the methods to set up the COWBAT pipeline databases.
        """
        logging.info('Beginning COWBAT database downloads')
        self._run_method_if_needed(self.sipprverse_targets, 'genesippr')
        self._run_method_if_needed(self.cowbat_targets, 'coregenome')
        self._run_method_if_needed(self.confindr_targets, 'ConFindr')
        self._run_method_if_needed(self.rmlst, 'rMLST', self.credentials)
        self._run_method_if_needed(self.mash, 'mash')
        self._run_method_if_needed(self._mlst_wrapper, 'MLST')
        self._run_method_if_needed(self.univec, 'univec')
        self._run_method_if_needed(
            self._cge_db_downloader_wrapper, 'resfinder')
        self._run_method_if_needed(
            self._cge_db_downloader_wrapper, 'virulence')
        self._run_method_if_needed(
            self._cge_db_downloader_wrapper, 'serosippr')
        self._run_method_if_needed(
            self._cge_db_downloader_wrapper, 'pointfinder')
        self.download_date()

    def _run_method_if_needed(self, method, folder, *args):
        """
        Run the provided method if the folder does not exist or overwrite is
        True.
        """
        if self.overwrite or not os.path.isdir(
                os.path.join(self.database_path, folder)):
            method(self.database_path, *args)

    def sipprverse_full(self):
        """
        Run a subset of the methods for the sipprverse full database setup.
        """
        logging.info('Beginning sipprverse full database downloads')
        self._run_method_if_needed(self.sipprverse_targets, 'genesippr')
        self._run_method_if_needed(self.confindr_targets, 'ConFindr')
        self._run_method_if_needed(self.mash, 'mash')
        self._run_method_if_needed(self.mlst, 'MLST')
        self._run_method_if_needed(self.rmlst, 'rMLST', self.credentials)
        self._run_method_if_needed(
            self._cge_db_downloader_wrapper, 'resfinder')
        self._run_method_if_needed(
            self._cge_db_downloader_wrapper, 'virulence')
        self._run_method_if_needed(
            self._cge_db_downloader_wrapper, 'serosippr')

    def sipprverse_method(self):
        """
        Run a reduced subset of methods for the sipprverse database setup.
        """
        logging.info('Beginning sipprverse method database downloads')
        self._run_method_if_needed(self.sipprverse_targets, 'genesippr')
        self._run_method_if_needed(self.confindr_targets, 'ConFindr')
        self._run_method_if_needed(self.mash, 'mash')

    def rmlst_method(self):
        """
        Run only the rMLST download.
        """
        self._run_method_if_needed(self.rmlst, 'rMLST', self.credentials)

    def sipprverse_targets(
        self,
        database_path,
        database_name='sipprverse',
        download_id='18130808'
    ):
        """
        Download OLC-specific sipprverse targets.
        """
        self.custom_databases(
            database_path=database_path,
            database_name=database_name,
            download_id=download_id,
            compression='gz'
        )

    def cowbat_targets(
        self,
        database_path,
        database_name='COWBAT',
        download_id='25319129',
    ):
        """
        Download OLC-specific COWBAT targets.
        """
        self.custom_databases(
            database_path=database_path,
            database_name=database_name,
            download_id=download_id,
            compression='gz'
        )

    def confindr_targets(self, database_name='ConFindr'):
        """
        Download OLC-specific ConFindr targets.

        Args:
            database_name (str): Name of the database directory to store
            ConFindr targets.
        """
        logging.info('Downloading ConFindr databases.')

        # Path to the secret file containing credentials
        secret_file = os.path.join(self.credentials, 'secret.txt')

        # Construct the ConFindr database setup command
        confindr_download = (
            f'confindr_database_setup -s {secret_file} -o '
            f'{os.path.join(self.database_path, database_name)}'
        )

        # Execute the ConFindr database setup command
        call(confindr_download, shell=True)

    @staticmethod
    def mlst(
        databasepath: str,
        genera: tuple = (
            'Escherichia', 'Salmonella', 'Yersinia', 'Campylobacter',
            'Cronobacter', 'Listeria', 'Bacillus', 'Staphylococcus', 'Vibrio'
        )
    ) -> None:
        """
        Download the up-to-date MLST profiles and alleles from PubMLST.

        Args:
            databasepath (str): Path to the directory to store MLST databases.
            genera (tuple): Tuple of genera to for which MLST profiles and
                alleles will be downloaded
        """
        logging.info('Downloading MLST databases from PubMLST')

        # Iterate through each genus in the provided genera tuple
        for genus in genera:
            args = CustomBox()
            args.genus = genus
            args.repository_url = 'http://pubmlst.org/data/dbases.xml'
            args.force_scheme_name = False
            args.path = os.path.join(databasepath, 'MLST', genus)

            # Path to the file indicating the download is complete
            complete_file = os.path.join(args.path, 'complete')

            # Only download if the complete file does not exist
            if not os.path.isfile(complete_file):
                # Download the MLST profiles and alleles for the genus
                get_mlst.main(args)

                # Create the complete file and write the list of downloaded
                # files
                with open(complete_file, 'w', encoding='utf-8') as complete:
                    complete.write(
                        '\n'.join(
                            glob(
                                os.path.join(args.path, '*')
                            )
                        )
                    )

    def enterobase_mlst(
        self,
        databasepath,
        genera=(
            'Escherichia',
            'Salmonella',
            'Yersinia'
        )
    ):
        """
        Download the necessary up-to-date MLST profiles and alleles from
        Enterobase.
        """
        logging.info('Downloading MLST databases from Enterobase')
        for genus in genera:
            self.enterobase_download(
                'MLST_Achtman', databasepath, 'MLST', genus)

    def enterobase_download(
            self,
            scheme: str,
            databasepath: str,
            analysis: str,
            genus: str) -> None:
        """
        Download the appropriate scheme (MLST_Achtman/cgMLST) for the
        requested genus.

        Args:
            scheme (str): The scheme to download (e.g., MLST_Achtman, cgMLST).
            databasepath (str): Path to the directory to store the downloaded
                scheme.
            analysis (str): Type of analysis (e.g., MLST, cgMLST).
            genus (str): Genus for which the scheme is to be downloaded.
        """
        # Get the genus name as recognized by EnteroBase
        enterobase_genus = self.genus_dict[genus]

        # Construct the EnteroBase API download command
        enterobase_cmd = (
            f'python -m olctools.database_setup.enterobase_api_download -o '
            f'{enterobase_genus} -s {scheme} -d {databasepath}'
        )

        # Check if the directory for the analysis and genus exists
        if not os.path.isdir(os.path.join(databasepath, analysis, genus)):
            # Run the EnteroBase API download command
            out, err = run_subprocess(enterobase_cmd)
            print(out, err)

    @staticmethod
    def rmlst(databasepath: str, credentials: str) -> None:
        """
        Get the most up-to-date profiles and alleles from PubMLST.

        Args:
            databasepath (str): Path to the directory to store rMLST databases.
            credentials (str): Path to the credentials file for PubMLST access.
        """
        logging.info('Downloading rMLST database')

        # Path to the file indicating the download is complete
        complete_file = os.path.join(databasepath, 'rMLST', 'complete')

        # Only download if the complete file does not exist
        if not os.path.isfile(complete_file):
            args = CustomBox()
            args.path = databasepath
            args.logging = logging
            args.credentials = credentials

            # Download the rMLST profiles and alleles
            get_rmlst.Get(args)

            # Create the complete file and write the list of downloaded files
            with open(complete_file, 'w', encoding='utf-8') as complete:
                complete.write(
                    '\n'.join(
                        glob(
                            os.path.join(databasepath, 'rMLST', '*')
                        )
                    )
                )

    def enterobase_cgmlst(
        self,
        databasepath,
        genera=(
            'Escherichia',
            'Yersinia')):
        """
        Download the necessary up-to-date cgMLST profiles and alleles from
        Enterobase.
        """
        logging.info('Downloading cMLST databases from Enterobase')
        for genus in genera:
            self.enterobase_download('cgMLST', databasepath, 'cgMLST', genus)

    def mash(self, database_path: str) -> None:
        """
        Download the pre-computed sketch of the RefSeq database, and compress
        it with gzip.

        Args:
            databasepath (str): Path to the directory to store MASH sketches.
        """
        logging.info('Downloading pre-computed RefSeq MASH sketches')

        # Create the MASH database folder
        database_path = self.create_database_folder(database_path, 'mash')

        # Path to the RefSeq assembly summary file
        output_file = os.path.join(
            database_path, 'assembly_summary_refseq.txt')

        # Download the assembly summary file if it doesn't exist
        if not os.path.isfile(output_file):
            self.database_download(
                output_file,
                'ftp://ftp.ncbi.nih.gov/genomes/ASSEMBLY_REPORTS/'
                'assembly_summary_refseq.txt',
                database_path)

        # Path to the RefSeq MASH sketches file
        output_file = os.path.join(database_path, 'RefSeqSketchesDefaults.msh')

        # Download the MASH sketches file if it doesn't exist
        if not os.path.isfile(output_file):
            self.database_download(
                output_file,
                'https://obj.umiacs.umd.edu/marbl_publications/mash/'
                'refseq.genomes.k21s1000.msh',
                database_path,
                complete=True)

    def univec(self, databasepath: str) -> None:
        """
        Download the UniVec core database.

        Args:
            databasepath (str): Path to the directory to store the UniVec
            database.
        """
        logging.info('Downloading univec database')

        # Create the UniVec database folder
        databasepath = self.create_database_folder(databasepath, 'univec')

        # Path to the UniVec core file
        output_file = os.path.join(databasepath, 'UniVec_core.tfa')

        # URL to download the UniVec core file
        target_url = 'ftp://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec_Core'

        # Download the UniVec core file
        self.database_download(output_file, target_url, databasepath)

        # Rename the downloaded file to have a .fasta extension
        if os.path.isfile(output_file):
            renamed = os.path.splitext(output_file)[0] + '.fasta'
            shutil.copy(output_file, renamed)

    @staticmethod
    def url_request(target_url: str, output_file: str, context: Any) -> None:
        """
        Use urllib to download the requested file from the target URL.

        Args:
            target_url (str): URL of the file to download.
            output_file (str): Path to save the downloaded file.
            context (Any): SSL context for secure connections.
        """
        # Open the URL and download the file in chunks
        request = urllib.request.urlopen(target_url, context=context)
        with open(output_file, 'wb') as targets:
            while True:
                data = request.read(4096)
                if not data:
                    break
                targets.write(data)

    def custom_databases(
        self,
        database_path: str,
        database_name: str,
        download_id: str,
        f_type: str = 'files',
        post_id: Optional[str] = None,
        compression: str = 'tar',
        nested: bool = False,
        complete: bool = False
    ) -> None:
        """
        Download and extract a .tar.gz file from figshare.

        Args:
            database_path (str): Path to the directory to store the database.
            database_name (str): Name of the database.
            download_id (str): Figshare download ID.
            f_type (str): Type of figshare file (default is 'files').
            post_id (Optional[str]): Post ID for figshare (default is None).
            compression (str): Compression type (default is 'tar').
            nested (bool): Whether the database is nested (default is False).
            complete (bool): Whether to create a 'complete' file after
                extraction (default is False).
        """
        logging.info('Downloading %s databases', database_name)

        # Adjust the database path if nested
        if nested:
            database_path = os.path.join(database_path, database_name)

        # Path to the 'complete' file indicating download completion
        complete_file = os.path.join(database_path, 'complete')

        # Create the database directory if it doesn't exist
        os.makedirs(database_path, exist_ok=True)

        # Path to the downloaded tar file
        tar_file = os.path.join(database_path, download_id)

        # Construct the target URL for downloading the file
        target_url = f'https://ndownloader.figshare.com/{f_type}/{download_id}'
        if post_id:
            target_url += f'/{post_id}'
        logging.debug(target_url)

        # Download the file if the 'complete' file doesn't exist
        if not os.path.isfile(complete_file):
            self.url_request(target_url, tar_file, self.context)

        # Decompress the downloaded file
        self.decompress(
            database_path=database_path,
            database_name=database_name,
            compression=compression,
            compressed_file=tar_file
        )

        # Create the 'complete' file if specified
        if complete:
            with open(complete_file, 'w', encoding='utf-8') as complete:
                complete.write('')

    @staticmethod
    def decompress(
        database_path: str,
        database_name: str,
        compression: str,
        compressed_file: str
    ) -> None:
        """
        Decompress the provided file using the appropriate library.

        Args:
            database_path (str): Path to the directory to store the
                decompressed files.
            database_name (str): Name of the database.
            compression (str): Compression type ('tar', 'gz', or 'zip').
            compressed_file (str): Path to the compressed file.
        """
        # Guard statement to check if the compressed file exists
        if not os.path.isfile(compressed_file):
            logging.warning(
                'Compressed file %s does not exist', compressed_file
            )
            raise FileNotFoundError

        # Decompress based on the specified compression type
        if compression == 'tar':
            logging.info('Extracting %s from archives', database_name)
            with tarfile.open(compressed_file, 'r') as tar:
                tar.extractall(path=database_path)
        elif compression == 'gz':
            logging.info('Extracting %s from archives', database_name)
            with gzip.open(compressed_file, 'rb') as gz:
                file_name = os.path.basename(
                    os.path.splitext(compressed_file)[0]
                )

                # Create the output directory if required
                specific_database_path = os.path.join(
                    database_path,
                    database_name
                )

                os.makedirs(specific_database_path, exist_ok=True)
                output_file = os.path.join(
                    specific_database_path, file_name
                )

                # Write the decompressed output to the file
                with open(output_file, 'wb') as output:
                    shutil.copyfileobj(gz, output)

        elif compression == 'zip':
            logging.info('Extracting %s from archives', database_name)
            with zipfile.ZipFile(compressed_file, 'r') as zip_file:
                zip_file.extractall(path=database_path)
        else:
            logging.error('Unsupported compression type: %s', compression)
            return

        # Remove the compressed file after extraction
        os.remove(compressed_file)

    def cge_db_downloader(
        self,
        database_path: str,
        analysistype: str,
        dbname: str,
        extension_in: str = 'fsa',
        extension_out: str = 'tfa'
    ) -> None:
        """
        Clones CGE databases into appropriate folder.

        Args:
            database_path (str): Path to the directory to store the database.
            analysistype (str): Type of analysis (e.g., serosippr).
            dbname (str): Name of the database to clone.
            extension_in (str): Input file extension (default is 'fsa').
            extension_out (str): Output file extension (default is 'tfa').
        """
        logging.info('Downloading %s database', analysistype)

        # Adjust the database path based on the analysis type
        if analysistype == 'serosippr':
            database_path = os.path.join(
                database_path, analysistype, 'Escherichia')
        else:
            database_path = os.path.join(database_path, analysistype)

        # Construct the git clone command
        target_call = (
            f'git clone https://bitbucket.org/genomicepidemiology/'
            f'{dbname}.git {database_path}'
        )

        # Clone the database repository
        self.database_clone(target_call, database_path)

        # Rename files if the input and output extensions differ
        if extension_in != extension_out:
            fasta_files = glob(
                os.path.join(database_path, f'*.{extension_in}')
            )
            for fasta in fasta_files:
                filename = os.path.splitext(fasta)[0]
                os.rename(fasta, f'{filename}.{extension_out}')

        # Combine target files into a single file if it doesn't already exist
        combined_file = os.path.join(database_path, 'combinedtargets.fasta')
        if not os.path.isfile(combined_file):
            database_files = glob(
                os.path.join(
                    database_path,
                    f'*.{extension_out}'))
            combinetargets(database_files, database_path)

    @staticmethod
    def create_database_folder(databasepath: str, database: str) -> str:
        """
        Create an appropriately named folder in which the database is to be
        stored.

        Args:
            databasepath (str): Path to the parent directory for the database.
            database (str): Name of the database.

        Returns:
            str: Path to the created database folder.
        """
        logging.info('Setting up %s database', database)

        # Construct the full path to the database folder
        database_folder_path = os.path.join(databasepath, database)

        # Guard statement to check if the path is valid
        if not database_folder_path:
            logging.error(
                'Invalid database folder path: %s', database_folder_path
            )
            return ''

        # Create the database folder if it doesn't exist
        os.makedirs(database_folder_path, exist_ok=True)

        return database_folder_path

    def database_download(
        self,
        output_file: str,
        target_url: str,
        database_path: str,
        complete: bool = False
    ) -> None:
        """
        Check to see if the download has previously been completed.

        Args:
            output_file (str): Path to save the downloaded file.
            target_url (str): URL of the file to download.
            database_path (str): Path to the directory to store the database.
            complete (bool): Whether to create a 'complete' file after download
        """
        # Path to the 'complete' file indicating download completion
        complete_file = os.path.join(database_path, 'complete')

        # Check if the 'complete' file exists
        if not os.path.isfile(complete_file):
            # Download the file from the target URL
            self.url_request(target_url, output_file, self.context)

            # Create the 'complete' file if specified
            if complete:
                with open(complete_file, 'w', encoding='utf-8') as complete_f:
                    complete_f.write('')

    @staticmethod
    def database_clone(
        target_call: str,
        database_path: str,
        complete: bool = False
    ) -> None:
        """
        Checks to see if the database has already been downloaded.

        Args:
            target_call (str): Command to clone the database repository.
            database_path (str): Path to the directory to store the database.
            complete (bool): Whether to create a 'complete' file after cloning
        """
        # Path to the 'complete' file indicating clone completion
        complete_file = os.path.join(database_path, 'complete')

        # Check if the 'complete' file exists
        if not os.path.isfile(complete_file):
            # Clone the database repository
            out, err = run_subprocess(target_call)

            # Create the 'complete' file if specified
            if complete:
                with open(complete_file, 'w', encoding='utf-8') as complete_f:
                    complete_f.write(out)
                    complete_f.write(err)

    def download_date(self):
        """
        Write the current date to file.
        """

        # Set the date file
        date_file = os.path.join(self.database_path, 'download_date')

        with open(date_file, 'w', encoding='utf-8') as download:
            download.write(f'{datetime.today():%Y-%m-%d}')

    def _mlst_wrapper(self, databasepath: str) -> None:
        """
        Wrapper for MLST download to handle Enterobase option.

        Args:
            databasepath (str): Path to the directory to store MLST databases.
        """
        # Check if Enterobase option is enabled
        if self.enterobase:
            # Download MLST databases using Enterobase
            self.enterobase_mlst(databasepath)
            # Download MLST databases from PubMLST
            self.mlst(databasepath)
        else:
            # Download MLST databases from PubMLST
            self.mlst(databasepath)

    def _cge_db_downloader_wrapper(
            self,
            databasepath: str,
            analysistype: str) -> None:
        """
        Wrapper for CGE database downloader to handle different analysis types.

        Args:
            databasepath (str): Path to the directory to store CGE databases.
            analysistype (str): Type of analysis (e.g., serosippr).
        """
        # Construct the database name based on the analysis type
        dbname = f'{analysistype}_db'

        # Download the CGE database for the specified analysis type
        self.cge_db_downloader(databasepath, analysistype, dbname)


# If the script is called from the command line, then call the argument parser
if __name__ == '__main__':
    parser = ArgumentParser(
        description='Downloads and sets up required databases'
    )
    parser.add_argument(
        '-d', '--database_path',
        required=True,
        help='Absolute path to location to store database files. '
        'Include any version numbers if required.'
    )
    parser.add_argument(
        '-c,', '--credentials',
        help='Name and path of folder containing required rMLST credentials.'
    )
    parser.add_argument(
        '-o', '--overwrite',
        action='store_true',
        help='Optionally allow for the overwriting of database files in '
        'the database path. Defaults to False, so if the output folder '
        'exists, that part of the download will be skipped.'
    )
    parser.add_argument(
        '-s', '--sipprverse_full',
        action='store_true',
        help='Optionally only download the databases used in the sipprverse. '
        'These include: genesippr, GDCS, sixteenS, ConFindr, MASH, MLST, '
        'rMLST, ResFindr, VirulenceFinder, and SerotypeFinder'
    )
    parser.add_argument(
        '-m', '--sipprverse_method',
        action='store_true',
        help='Optionally only download the databases used by the sipprverse '
        'method: genesippr, sixteenS, GDCS, MASH, and ConFindr'
    )
    parser.add_argument(
        '-e', '--enterobase',
        action='store_true',
        help='Use Enterobase to download MLST definitions for Escherichia, '
        'Salmonella, and Yersinia, as well as cgMLST schemes for '
        'Escherichia and Yersinia. Disabled by default'
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Option to include debug level logging messages.'
    )
    parser.add_argument(
        '-r', '--rmlst',
        action='store_true',
        help='Optionally only download the rMLST database'
    )
    parser.add_argument(
        '-res', '--resfinder',
        action='store_true',
        help='Only download the Resfinder database'
    )
    arguments = parser.parse_args()

    # Create a logger
    logging_level = 'INFO' if not arguments.verbose else 'DEBUG'
    logger = logging.getLogger()
    logger.setLevel(getattr(logging, logging_level))

    # Create a formatter without milliseconds
    formatter = logging.Formatter(
        '%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

    # Create an object
    pipeline = DatabaseSetup(
        database_path=arguments.database_path,
        credentials=arguments.credentials,
        overwrite=arguments.overwrite,
        enterobase=arguments.enterobase
    )

    # Run the appropriate analyses
    if arguments.resfinder:
        pipeline.cge_db_downloader(
            database_path=pipeline.database_path,
            analysistype='resfinder',
            dbname='resfinder_db'
        )
        raise SystemExit

    # Raise an error if the credentials file was not provided
    if not arguments.credentials:
        logging.error(
            'Please provide the name and path of the folder containing your '
            'rMLST credentials with the -c argument.'
        )
        raise SystemExit

    # Run the appropriate method
    if arguments.sipprverse_full:
        pipeline.sipprverse_full()
    elif arguments.sipprverse_method:
        pipeline.sipprverse_method()
    elif arguments.rmlst:
        pipeline.rmlst_method()
    else:
        pipeline.cowbat()
