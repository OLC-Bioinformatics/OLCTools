#!/usr/bin/env python
from olctools.accessoryFunctions.accessoryFunctions import make_path, SetupLogging
from Bio import SeqIO
from argparse import ArgumentParser
import logging
import os


class Split(object):

    def read_split_file(self):
        # Create a FASTA file for each record in the file
        for record in SeqIO.parse(self.allele_file, 'fasta'):
            # Create the file. Clean up the name a bit
            record.id = record.id.replace(':', '-').replace('_', '-')
            with open(os.path.join(self.output_path, f'{record.id}{self.extension}'), 'w') as split:
                # Update the record.id with a fake allele number to resemble an allele file
                if self.number:
                    record.id = record.id + '_0'
                    record.description = str()
                # Write the record in FASTA format to the allele file
                SeqIO.write(sequences=record,
                            handle=split,
                            format='fasta')

    def __init__(self, allele_file, extension='.fasta', number=True):
        if allele_file.startswith('~'):
            self.allele_file = os.path.abspath(os.path.expanduser(os.path.join(allele_file)))
        else:
            self.allele_file = os.path.abspath(os.path.join(allele_file))
        assert os.path.isfile(self.allele_file), f'Cannot locate the supplied allele file: {self.allele_file}'
        self.output_path = os.path.join(os.path.dirname(self.allele_file), 'alleles')
        make_path(self.output_path)
        self.extension = extension
        self.number = number


def cli():
    # Parser for arguments
    parser = ArgumentParser(description='Split each gene in a combined FASTA file into individual files')
    parser.add_argument('-a', '--allele_file',
                        required=True,
                        help='Specify name and path of allele file')
    parser.add_argument('-e', '--extension',
                        choices=['.fasta', '.fa', '.tfa'],
                        default='.fasta',
                        help='Set the file extension to use for the split files')
    parser.add_argument('-n', '--number',
                        action='store_true',
                        help='Do not add an _0 to the end of the gene name in the FASTA header to resemble an '
                             'allele file')
    # Get the arguments into an object
    args = parser.parse_args()
    SetupLogging(debug=True)
    reformat = Split(allele_file=args.allele_file,
                     extension=args.extension,
                     number=args.number)
    reformat.read_split_file()
    logging.info('Gene splitting complete!')


if __name__ == '__main__':
    cli()
