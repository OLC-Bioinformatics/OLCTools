#!/usr/bin/env python
from olctools.accessoryFunctions.accessoryFunctions import make_path, SetupLogging
from Bio import SeqIO
from argparse import ArgumentParser
import logging
import os


class Split(object):

    def read_split_file(self):
        for record in SeqIO.parse(self.allele_file, "fasta"):
            with open(os.path.join(self.output_path, f'{record.id.replace(":", "_")}.fasta'), 'w') as split:
                SeqIO.write(sequences=record,
                            handle=split,
                            format='fasta')

    def __init__(self, allele_file):
        if allele_file.startswith('~'):
            self.allele_file = os.path.abspath(os.path.expanduser(os.path.join(allele_file)))
        else:
            self.allele_file = os.path.abspath(os.path.join(allele_file))
        assert os.path.isfile(self.allele_file), f'Cannot locate the supplied allele file: {self.allele_file}'
        self.output_path = os.path.join(os.path.dirname(self.allele_file), 'outputs')
        make_path(self.output_path)


def cli():
    # Parser for arguments
    parser = ArgumentParser(description='Split each gene in a combined FASTA file into individual files')
    parser.add_argument('-a', '--allele_file',
                        required=True,
                        help='Specify name and path of allele file')
    # Get the arguments into an object
    args = parser.parse_args()
    SetupLogging(debug=True)
    reformat = Split(allele_file=args.allele_file)
    reformat.read_split_file()
    logging.info('Gene splitting complete!')


if __name__ == '__main__':
    cli()
