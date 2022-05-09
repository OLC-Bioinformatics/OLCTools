from olctools.accessoryFunctions.accessoryFunctions import make_path, SetupLogging
import olctools.databasesetup.settings
try:
    from olctools.databasesetup.settings import ECOLI
except (NameError, ImportError):
    ECOLI = str()
try:
    from olctools.databasesetup.settings import SENTERICA
except (NameError, ImportError):
    SENTERICA = str()
try:
    from olctools.databasesetup.settings import YERSINIA
except (NameError, ImportError):
    YERSINIA = str()
from urllib3.exceptions import HTTPError
from argparse import ArgumentParser
from subprocess import call
import logging
import urllib3
import shutil
import json
import ssl
import os


class DownloadScheme(object):

    def main(self):
        self.download_profile()
        if self.name_file:
            self.read_names()
        self.download_alleles()

    def create_request(self, request_str):
        http = urllib3.PoolManager(cert_reqs='CERT_NONE')
        headers = urllib3.util.make_headers(basic_auth='{token}:'.format(token=self.api_key))
        request = http.request(method='GET', url=request_str, headers=headers, preload_content=False)
        return request

    def read_names(self):
        """
        Read in all the names of the genes of interest
        """
        logging.info('Reading names file')
        with open(self.name_file, 'r') as names:
            self.names = [name.rstrip() for name in names.readlines()]

    def download_profile(self):
        logging.info('Downloading {genus} {scheme} profile'.format(genus=self.genus,
                                                                   scheme=self.scheme))
        if self.organism == 'senterica':
            profile_output = os.path.join(self.output_path, '{scheme}-profiles.list.gz'.format(scheme=self.scheme))
            profile_list = os.path.join(self.output_path, '{scheme}-profiles.list'.format(scheme=self.scheme))
        else:
            profile_output = os.path.join(self.output_path, '{scheme}-profiles.gz'.format(scheme=self.scheme))
            profile_list = os.path.join(self.output_path, '{scheme}-profiles'.format(scheme=self.scheme))
        profile_text = os.path.join(self.output_path, 'profile.txt')
        if not os.path.isfile(profile_output) and not os.path.isfile(profile_text):
            address = '{address}{org}/schemes?scheme_name={sn}&limit={limit}&only_fields=download_sts_link' \
                .format(address=self.server_address,
                        org=self.organism,
                        sn=self.scheme,
                        limit=400000)
            try:
                response = self.create_request(address)
                try:
                    data = json.loads(response.data.decode('utf-8'))
                except json.decoder.JSONDecodeError:
                    data = dict()
                    print(response.data.decode('utf-8'))
                    quit()
                logging.debug(json.dumps(data, sort_keys=True, indent=4, separators=(',', ': ')))
                response.release_conn()
                for scheme_record in data['Schemes']:
                    profile_link = scheme_record.get('download_sts_link', None)
                    if profile_link:
                        logging.info('Downloading {genus} profiles from {profile_link} to: {profile}'
                                     .format(genus=self.genus,
                                             profile_link=profile_link,
                                             profile=profile_output))
                        profile_response = self.create_request(profile_link)
                        with open(profile_output, 'wb') as output_profile:
                            while True:
                                profile_data = profile_response.read()
                                if not profile_data:
                                    break
                                output_profile.write(profile_data)
            except (KeyError, HTTPError) as Response_error:
                error_string = str()
                for key, value in vars(Response_error).items():
                    error_string += '{key}: {value}\n'.format(key=key,
                                                              value=value)
                print('HTTPError: {error}'.format(error=error_string))
                quit()
        if os.path.isfile(profile_output) and not os.path.isfile(profile_list) and not os.path.isfile(profile_text):
            logging.info('Decompressing {allele}'.format(allele=profile_output))
            pigz_cmd = 'pigz -d -f {profile_output}'.format(profile_output=profile_output)
            call(pigz_cmd, shell=True)
        try:
            shutil.move(profile_list, profile_text)
        except (FileNotFoundError, FileExistsError):
            pass

    def download_alleles(self):
        logging.info('Downloading {genus} {scheme} alleles'.format(genus=self.genus,
                                                                   scheme=self.scheme))
        address = '{address}{org}/{sn}/loci?&limit={limit}&scheme={sn}'.format(address=self.server_address,
                                                                               org=self.organism,
                                                                               sn=self.scheme,
                                                                               limit=400000)
        try:
            response = self.create_request(address)
            try:
                data = json.loads(response.data.decode('utf-8'))
            except json.decoder.JSONDecodeError:
                data = dict()
                print(response.data.decode('utf-8'))
                quit()
            logging.debug(json.dumps(data, sort_keys=True, indent=4, separators=(',', ': ')))
            response.release_conn()
            for locus_record in data['loci']:
                locus_link = locus_record['download_alleles_link']
                locus_file_name = locus_link.split('/')[-1]
                print(locus_file_name)
                quit()
                locus_output = os.path.join(self.output_path, locus_file_name)
                locus_file = os.path.splitext(locus_output)[0]
                locus_tfa = locus_file.replace('.fasta', '.tfa')
                if locus_link:
                    if not os.path.isfile(locus_output) and not os.path.isfile(locus_file) \
                            and not os.path.isfile(locus_tfa):
                        logging.info('Downloading {scheme} allele {allele}'.format(scheme=self.scheme,
                                                                                   allele=locus_file))
                        locus_reponse = self.create_request(locus_link)
                        with open(locus_output, 'wb') as output_locus:
                            while True:
                                locus_data = locus_reponse.read()
                                if not locus_data:
                                    break
                                output_locus.write(locus_data)
                if os.path.isfile(locus_output) and not os.path.isfile(locus_file):
                    logging.info('Decompressing {allele}'.format(allele=locus_file))
                    pigz_cmd = 'pigz -d {gz_file}'.format(gz_file=locus_output)
                    call(pigz_cmd, shell=True)
                if os.path.isfile(locus_file) and not os.path.isfile(locus_tfa):
                    shutil.move(locus_file, locus_tfa)
        except (KeyError, HTTPError) as Response_error:
            error_string = str()
            for key, value in vars(Response_error).items():
                error_string += '{key}: {value}\n'.format(key=key,
                                                          value=value)
            print('HTTPError: {error}'.format(error=error_string))
            quit()

    def __init__(self, databasepath, organism, scheme, genes=None):
        self.server_address = 'http://enterobase.warwick.ac.uk/api/v2.0/'
        self.organism = organism
        self.scheme = scheme
        if databasepath.startswith('~'):
            self.databasepath = os.path.expanduser(os.path.abspath(os.path.join(databasepath)))
        else:
            self.databasepath = os.path.abspath(os.path.join(databasepath))
        genus_dict = {
            'ecoli': 'Escherichia',
            'senterica': 'Salmonella',
            'yersinia': 'Yersinia'
        }
        self.genus = genus_dict[self.organism]
        self.output_path = os.path.join(self.databasepath, self.scheme.split('_')[0], self.genus)
        make_path(self.output_path)
        if genes:
            if genes.startswith('~'):
                self.name_file = os.path.abspath(os.path.expanduser(os.path.join(genes)))
            else:
                self.name_file = os.path.abspath(os.path.join(genes))
            assert os.path.isfile(self.name_file), f'Cannot find the supplied file with gene names: {self.name_file}'
        self.names = list()
        if self.organism == 'ecoli':
            self.api_key = ECOLI
            if not self.api_key:
                # Use the user input to set the verifier code
                self.api_key = input('Enter API token from http://enterobase.warwick.ac.uk/species/index/ecoli ')
                try:
                    with open(olctools.databasesetup.settings.__file__, 'a+') as env:
                        env.write("ECOLI = '{api}'\n".format(api=self.api_key))
                except:
                    raise
        elif self.organism == 'senterica':
            self.api_key = SENTERICA
            if not self.api_key:
                # Use the user input to set the verifier code
                self.api_key = input('Enter API token from http://enterobase.warwick.ac.uk/species/index/senterica ')
                try:
                    with open(olctools.databasesetup.settings.__file__, 'a+') as env:
                        env.write("SENTERICA = '{api}'\n".format(api=self.api_key))
                except:
                    raise
        else:
            self.api_key = YERSINIA
            if not self.api_key:
                # Use the user input to set the verifier code
                self.api_key = input('Enter API token from http://enterobase.warwick.ac.uk/species/index/yersinia ')
                try:
                    with open(olctools.databasesetup.settings.__file__, 'a+') as env:
                        env.write("YERSINIA = '{api}'\n".format(api=self.api_key))
                except:
                    raise


def cli():
    # Parser for arguments
    parser = ArgumentParser(description='Download typing schemes and alleles from Enterobase')
    parser.add_argument('-d', '--databasepath',
                        required=True,
                        help='The path to the folder in which the typing scheme is to be installed. The program will '
                             'create sub-folders as necessary. So, if you specify '
                             '/mnt/nas2/databases/assemblydatabases/0.5.0.0, that will be used as the root for the '
                             'SCHEME/ORGANISM subfolder, e.g. cgMLST/Escherichia')
    parser.add_argument('-o', '--organism',
                        required=True,
                        choices=['ecoli', 'senterica', 'yersinia'])
    parser.add_argument('-s', '--scheme',
                        required=True,
                        choices=['MLST_Achtman', 'cgMLST', 'wgMLST'])
    parser.add_argument('-g', '--genes',
                        help='Name and path to a file containing the gene names (one per line) to be extracted '
                             'from the profile')
    parser.add_argument('-v', '--verbose',
                        action='store_true',
                        help='Print debug level messages')
    # Get the arguments into an object
    arguments = parser.parse_args()
    # Setup logging
    SetupLogging(debug=arguments.verbose)
    download = DownloadScheme(databasepath=arguments.databasepath,
                              organism=arguments.organism,
                              scheme=arguments.scheme,
                              genes=arguments.genes)
    download.download_profile()
    download.download_alleles()


# If the script is called from the command line, then call the argument parser
if __name__ == '__main__':
    cli()
