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
import urllib3
import json
import os


class DownloadScheme(object):

    def create_request(self, request_str):
        http = urllib3.PoolManager()
        headers = urllib3.util.make_headers(basic_auth='{token}:'.format(token=self.api_key))
        request = http.request(method='GET', url=request_str, headers=headers, preload_content=False)
        return request

    def download_assemblies(self):
        try:
            response = self.create_request(self.server_address)
            try:
                data = json.loads(response.data.decode('utf-8'))
            except json.decoder.JSONDecodeError:
                data = dict()
                print(response.data.decode('utf-8'))
                quit()
            print(json.dumps(data, sort_keys=True, indent=4, separators=(',', ': ')))
            response.release_conn()
            for record in data['straindata']:
                record_values = data['straindata'][record]
                assembly_response = self.create_request(record_values['download_fasta_link'])
                assembly_file = os.path.join(self.outputpath, '{sn}.fasta'.format(sn=record_values['assembly_barcode']))
                if not os.path.isfile(assembly_file):
                    with open(assembly_file, 'wb') as out_assembly:
                        out_assembly.write(assembly_response.read())
        except HTTPError as Response_error:
            error_string = str()
            for key, value in vars(Response_error).items():
                error_string += '{key}: {value}\n'.format(key=key,
                                                          value=value)
            print('HTTPError: {error}'.format(error=error_string))
            quit()

    def __init__(self, serovar, organism, outputpath):
        self.serovar = serovar
        self.organism = organism
        self.server_address = 'http://enterobase.warwick.ac.uk/api/v2.0/{organism}/straindata?serotype={sero}' \
                              '&assembly_status=Assembled&limit={limit}&only_fields=strain_name,download_fasta_link'\
            .format(organism=self.organism,
                    sero=self.serovar,
                    limit=1000)
        if outputpath.startswith('~'):
            self.outputpath = os.path.expanduser(os.path.abspath(os.path.join(outputpath)))
        else:
            self.outputpath = os.path.abspath(os.path.join(outputpath))
        make_path(self.outputpath)
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
    parser.add_argument('-o', '--outputpath',
                        # required=True,
                        default='/mnt/nas2/processed_sequence_data/enterobase_assemblies/Litchfield',
                        help='The path to the folder in which the typing scheme is to be installed. The program will '
                             'create sub-folders as necessary. So, if you specify '
                             '/mnt/nas2/databases/assemblydatabases/0.5.0.0, that will be used as the root for the '
                             'SCHEME/ORGANISM subfolder, e.g. cgMLST/Escherichia')
    parser.add_argument('-g', '--genus',
                        default='senterica',
                        choices=['ecoli', 'senterica', 'yersinia'])
    parser.add_argument('-s', '--serovar',
                        # required=True,
                        default='Litchfield',
                        choices=['Litchfield'])
    parser.add_argument('-v', '--verbose',
                        action='store_true',
                        help='Print debug level messages')
    # Get the arguments into an object
    arguments = parser.parse_args()
    # Setup logging
    SetupLogging(debug=arguments.verbose)
    download = DownloadScheme(outputpath=arguments.outputpath,
                              organism=arguments.genus,
                              serovar=arguments.serovar)
    download.download_assemblies()


# If the script is called from the command line, then call the argument parser
if __name__ == '__main__':
    cli()
