#!/usr/bin/env python3
from olctools.accessoryFunctions.accessoryFunctions import SetupLogging
from argparse import ArgumentParser
import logging
import pandas
import os
__author__ = 'adamkoziol'


class Check(object):

    def main(self):
        self.load_biosample_table()
        self.load_curation_table()
        self.load_sra_table()
        self.compare_tables()

    def load_biosample_table(self):
        """
        Use pandas to load the Biosample Excel file. Extract the values for the accession, Title, and OrganismName
        columns
        """
        logging.info('Loading Biosample table')
        # Use pandas to create a nested dictionary for the biosample table
        biosample_dictionary = self.load_table(self.biosample_table)
        # Iterate through the dictionary, and extract values for the 'accession', 'Title', and 'OrganismName' columns
        for primary_key, header in biosample_dictionary.items():
            # Add the accession, genus and title values to the dictionary
            self.biosample_title_dict[primary_key] = {
                'title': header['Title'],
                'accession': header['accession'],
                'organism': header['OrganismName']
            }

    def load_curation_table(self):
        """
        Parse the Access database outputs with pandas, and create a dictionary containing the values for the required
        keys
        """
        logging.info('Loading Access database curation table')
        # Use pandas to create a nested dictionary for the curation table
        curation_dictionary = self.load_table(self.curation_table)
        # Iterate through the dictionary, and add the required key: value pairs to the dictionary
        for primary_key, header in curation_dictionary.items():
            self.curation_title_dict[primary_key] = {
                'title': header['CFIAID'],
                'seqid': header['SEQID'],
                'genus': header['Genus'],
                'species': header['Species'],
                'mash': str(header['MASH_ReferenceGenome']),
                'curator_flag': header['CuratorFlag']
            }

    def load_sra_table(self):
        """
        Parse the SRA metadata table with pandas, and create a dictionary containing the values for the required keys
        """
        logging.info('Loading SRA metadata table')
        # Use pandas to create a nested dictionary for the SRA metadata table
        sra_dictionary = self.load_table(self.sra_table,
                                         sheet_name='SRA_data')
        # Iterate through the dictionary, and add the required key: value pairs to the dictionary
        for primary_key, header in sra_dictionary.items():
            # The 'title' column contains the CFIA ID joined to the organism e.g. CFIAFB20180146_L. monocytogenes
            # split on the underscore
            title, organism = header['title'].split('_')
            self.sra_title_dict[primary_key] = {
                'title': title,
                'seqid': header['library_ID'],
                'accession': header['biosample_accession'],
                'organism': organism,
                'forward': header['filename'],
                'reverse': header['filename2']
            }

    def compare_tables(self):
        """

        """
        logging.info('Finding errors and discrepancies between tables')
        curation_fails = set()
        curation_passes = set()
        for sra_key, sra_dict in self.sra_title_dict.items():
            # Extract the CFIA ID, SEQ ID and accession used in the SRA table
            cfia_id = sra_dict['title']
            seqid = sra_dict['seqid']
            accession = sra_dict['accession']
            # Check to see if the genus is the right one, and that it has been entered correctly
            try:
                assert self.organism_conversion_dict[self.genus] in sra_dict['organism']
            except AssertionError:
                logging.warning('SRA incorrect organism:', cfia_id, self.genus, sra_dict['organism'])
            # Ensure that the file names are formatted properly
            try:
                assert sra_dict['seqid'] in sra_dict['forward'] and sra_dict['seqid'] in sra_dict['reverse']
            except AssertionError:
                logging.warning('Filename mismatch:', cfia_id, self.genus, sra_dict['seqid'], sra_dict['forward'],
                                sra_dict['reverse'])
            # Compare the SRA table entries against the BioSample table
            bio_sample_present = False
            for biosample_key, biosample_dict in self.biosample_title_dict.items():
                # Find the matching CFIA ID
                if biosample_dict['title'] == cfia_id:
                    bio_sample_present = True
                    # Make sure the genus is correct
                    try:
                        assert self.genus in biosample_dict['organism']
                    except AssertionError:
                        logging.warning('Genus mismatch biosample table:',  cfia_id, self.genus,
                                        biosample_dict['organism'])
                    # Ensure that the accessions match
                    try:
                        assert accession == biosample_dict['accession']
                    except AssertionError:
                        logging.warning('Accession mismatch biosample table:', cfia_id, sra_dict['accession'],
                                        biosample_dict['accession'])
            # Indicate that the BioSample table is missing an entry
            if not bio_sample_present:
                logging.warning('No entry in the Biosample table:', cfia_id, sra_dict['seqid'], sra_dict['accession'])
            # Compare the SRA table and the Access database curation table
            curation_present = False
            for curation_key, curation_dict in self.curation_title_dict.items():
                curation_pass = True
                # Find the matching CFIA ID in the curation table. Note that the table can have multiple samples with
                # the same ID
                if curation_dict['title'] == cfia_id:
                    curation_present = True
                    # Check to see if there is an entry with a matching SEQ ID
                    try:
                        assert seqid == curation_dict['seqid']
                    except AssertionError:
                        curation_pass = False
                    # Ensure that the genus is correct
                    try:
                        assert self.genus in curation_dict['genus']
                    except AssertionError:
                        curation_pass = False
                    # Ensure that the MASH prediction is correct
                    try:
                        assert self.genus in curation_dict['mash']
                    except AssertionError:
                        curation_pass = False
                    # Confirm that the sample passes curation
                    try:
                        assert 'REFERENCE' in curation_dict['curator_flag'] or 'PASS' in curation_dict['curator_flag']
                    except AssertionError:
                        curation_pass = False
                    # Add the entry to the appropriate list
                    if not curation_pass:
                        curation_fails.add(cfia_id)
                    else:
                        curation_passes.add(cfia_id)
            # Indicate that there was no entry in the curation table
            if not curation_present:
                logging.warning('No entry in the curation table:  ', cfia_id, 'SEQ ID:', sra_dict['seqid'], 'Accession:',
                                sra_dict['accession'])
        # Iterate through all the entries that
        for cfia_id in sorted(curation_fails - curation_passes):
            logging.warning('Curation fail:', cfia_id)
            # Iterate through the SRA table to return the correct entry
            for sra_key, sra_dict in self.sra_title_dict.items():
                if sra_dict['title'] == cfia_id:
                    logging.warning('\tSRA entry:     ', cfia_id, 'SEQ ID:', sra_dict['seqid'], 'Accession:',
                                    sra_dict['accession'])
            # Iterate through the Access curation table, and return the entry/entries corresponding to the CFIA ID
            for curation_key, curation_dict in self.curation_title_dict.items():
                if curation_dict['title'] == cfia_id:
                    logging.warning('\tCuration entry:', cfia_id, 'SEQ ID:', curation_dict['seqid'], 'Genus:',
                                    curation_dict['genus'], 'MASH genus:', curation_dict['mash'],
                                    'Curator flag:', curation_dict['curator_flag'])

    @staticmethod
    def load_table(table_name, sheet_name=None):
        """
        Create nested dictionary from Excel file using pandas
        :param: table_name: Name and path of Excel file to process
        :param: sheet_name: Optional name of sheet in the Excel file to use. If not supplied, pandas will read the first
        sheet
        :return nested dictionary of primary key: header: value
        """
        # A dictionary to store the parsed excel file in a more readable format
        nesteddictionary = dict()
        if sheet_name:
            dictionary = pandas.read_excel(table_name, sheet_name=sheet_name).to_dict()
        else:
            dictionary = pandas.read_excel(table_name).to_dict()
        for header in dictionary:
            # primary_key is the primary key, and value is the value of the cell for that key + header combination
            for primary_key, value in dictionary[header].items():
                # Update the dictionary with the new data
                try:
                    nesteddictionary[primary_key].update({header: value})
                # Create the nested dictionary if it hasn't been created yet
                except KeyError:
                    nesteddictionary[primary_key] = dict()
                    nesteddictionary[primary_key].update({header: value})
        return nesteddictionary

    def __init__(self, path, biosampletable, curationtable, sratable, genus):
        if path.startswith('~'):
            self.path = os.path.abspath(os.path.expanduser(os.path.join(path)))
        else:
            self.path = os.path.abspath(os.path.join(path))
        self.biosample_table = os.path.join(self.path, biosampletable)
        self.curation_table = os.path.join(self.path, curationtable)
        self.sra_table = os.path.join(self.path, sratable)
        self.genus = genus
        self.biosample_title_dict = dict()
        self.curation_title_dict = dict()
        self.sra_title_dict = dict()
        self.organism_conversion_dict = {
            'Camplylobacter': 'C. jejuni',
            'Escherichia': 'E. coli',
            'Listeria': 'L. monocytogenes',
            'Salmonella': 'S. enterica',
        }
        SetupLogging(filehandle=os.path.join(self.path, 'log.txt'),
                     log_level=logging.WARN)


def cli():
    parser = ArgumentParser(description='Confirm that SRA metadata table has correct biosample accessions, SEQIDs, '
                                        'only contains the genus of interest, and passes curation')
    parser.add_argument('-p', '--path',
                        required=True,
                        help='Path to folder containing necessary tables')
    parser.add_argument('-b', '--biosampletable',
                        required=True,
                        help='Name of Biosample table from NCBI (must be in the supplied path)')
    parser.add_argument('-c', '--curationtable',
                        required=True,
                        help='Name of Curation table from Access database(must be in the supplied path)')
    parser.add_argument('-s', '--sratable',
                        required=True,
                        help='Name of SRA metadata table (must be in the supplied path)')
    parser.add_argument('-g', '--genus',
                        choices=['Campylobacter', 'Escherichia', 'Listeria', 'Salmonella'],
                        required=True,
                        help='Genus of samples being uploaded')
    arguments = parser.parse_args()
    check = Check(path=arguments.path,
                  biosampletable=arguments.biosampletable,
                  curationtable=arguments.curationtable,
                  sratable=arguments.sratable,
                  genus=arguments.genus)
    check.main()


if __name__ == '__main__':
    cli()
