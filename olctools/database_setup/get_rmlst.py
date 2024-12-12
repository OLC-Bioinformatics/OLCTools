#!/usr/bin/env python3

"""
Wrapper script to download the most up-to-date rMLST profile and alleles
and combine them into a single file
"""

# Standard imports
from argparse import ArgumentParser
from glob import glob
import logging
import os
from typing import List

# Third party imports
from Bio import SeqIO

# Local imports
from olctools.database_setup import rest_auth_class

__author__ = 'adamkoziol'


class Get(object):
    """
    Class to download the most up-to-date rMLST profile and alleles
    """

    def get_rmlst_helper(self):
        """
        Makes a system call to rest_auth.py, a Python script modified from
        https://github.com/kjolley/BIGSdb/tree/develop/scripts/test
        And downloads the most up-to-date rMLST profile and alleles
        """
        # Set the path/name of the folder to contain the new alleles and
        # profile
        new_folder = os.path.join(self.path, self.analysistype)

        # Create the path
        os.makedirs(new_folder, exist_ok=True)

        # Create arguments to feed into the rest_auth_class script
        args = ArgumentParser
        args.secret_file = os.path.join(self.credentials, 'secret.txt')
        args.file_path = os.path.join(self.credentials)
        args.output_path = new_folder
        args.logging = self.logging
        rmlst = rest_auth_class.REST(args)

        # Download the profile and alleles
        rmlst.main()

        # Get the new alleles into a list, and create the combinedAlleles file
        alleles = glob(os.path.join(new_folder, '*.tfa'))
        self.combine_alleles(
            allele_path=new_folder,
            alleles=alleles
        )

    @staticmethod
    def combine_alleles(
        *,  # Force the use of named arguments
        allele_path: str,
        alleles=List
    ) -> None:
        """
        Combine the individual rMLST allele files into a single file
        :param allele_path: Path to the folder containing the rMLST alleles
        :param alleles: List of the rMLST alleles
        :return: None
        """
        # Set the path to the combinedAlleles file
        combined_alleles = os.path.join(allele_path, 'rMLST_combined.fasta')

        # Open the combinedAlleles file
        with open(combined_alleles, 'w', encoding='utf-8') as combined_file:
            # Open each allele file
            for allele_file in sorted(alleles):
                # Open each allele file
                with open(allele_file, 'r', encoding='utf-8') as fasta:
                    # Extract the sequence record from each entry in the
                    # multifasta
                    for record in SeqIO.parse(fasta, "fasta"):
                        # Replace and dashes in the record.id with underscores
                        record.id = record.id.replace('-', '_')
                        # Remove and dashes or 'N's from the sequence data as
                        # makeblastdb can't handle sequences with gaps
                        try:
                            record.seq._data = record.seq._data.replace(
                                '-', ''
                            ).replace('N', '')
                        except TypeError:
                            record.seq._data = record.seq._data.replace(
                                b'-', b''
                            ).replace(b'N', b'')
                        # Clear the name and description attributes of the
                        # record
                        record.name = ''
                        record.description = ''
                        # Write each record to the combined file
                        SeqIO.write(record, combined_file, 'fasta')

    def __init__(self, args):
        self.path = os.path.join(args.path)
        os.makedirs(self.path, exist_ok=True)
        self.logging = args.logging
        self.analysistype = 'rMLST'
        self.credentials = args.credentials
        self.get_rmlst_helper()


if __name__ == '__main__':
    # Argument parser for user-inputted values, and a nifty help menu
    # Parser for arguments
    parser = ArgumentParser(description='')
    parser.add_argument(
        '--path',
        required=True,
        help='Specify input directory'
    )
    parser.add_argument(
        '--credentials',
        required=True,
        help='Specify directory storing rMLST credentials'
    )
    # Get the arguments into an object
    arguments = parser.parse_args()
    arguments.logging = logging
    # Run the script
    Get(arguments)
