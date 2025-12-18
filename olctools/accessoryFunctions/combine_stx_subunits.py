#!/usr/bin/env python
from olctools.accessoryFunctions.accessoryFunctions import make_path, SetupLogging
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from argparse import ArgumentParser
import logging
import os


class Reformat(object):

    def main(self):
        self.read_file()
        self.find_accessions()
        self.match_accessions()

    def read_file(self):
        self.record_dict = SeqIO.to_dict(SeqIO.parse(self.allele_file, "fasta"))

    def find_accessions(self):
        for record in self.record_dict:
            # Create an empty string to store the subunit information. Sometimes, it is not included in the gene header
            subunit = str()
            # Split the header information e.g. stx1A:1:AB050958:d yields gene: stx1A, allele: 1, accession: AB050958,
            # subunit: d
            try:
                gene, allele, accession, subunit = record.split(':')
            # If there is no subunit information e.g. stx1A:13:JQ327854 yields gene: stx1A, allele 13, accession
            # JQ327854
            except ValueError:
                gene, allele, accession = record.split(':')
            # Populate accession_dict with accession: gene: string of gene sequence as required
            if accession not in self.accession_dict:
                self.accession_dict[accession] = dict()
            if gene not in self.accession_dict[accession]:
                self.accession_dict[accession].update({
                    gene: str(self.record_dict[record].seq),
                }
                )
            # Store the link between accession and subunit, and accession and allele
            self.subunit_dict[accession] = subunit
            self.allele_dict[accession] = allele

    def match_accessions(self):
        """
        Find allele subunits with matching accessions
        """
        # Create a list to store SeqRecords if a combined file is being written at the end
        allele_list = list()
        for accession, allele_dict in self.accession_dict.items():
            # Initialise variables to store gene and sequence information
            gene = str()
            seq_str = str()
            # If the length of the dictionary is two, then there are two matching subunits
            if len(allele_dict) == 2:
                # Extract the gene (plus subunit) and its corresponding sequence from the sorted dictionary - the
                # sorting ensures that subunits are in order e.g. stx2A then stx2B
                for gene_subunit, seq in sorted(allele_dict.items()):
                    gene = gene_subunit[:-1]
                    seq_str += seq
            # Only accession with both subunits need to be written to file
            if gene:
                # Create a string to store all the necessary information for each gene
                gene_id = '{gene}_{acc}_{sub}'.format(gene=gene,
                                                      acc=accession,
                                                      sub=self.subunit_dict[accession]
                                                      )
                # Create a SeqRecord of the sequence
                seq_record = SeqRecord(seq=Seq(seq_str),
                                       id=gene_id,
                                       name='',
                                       description='')
                # Write the SeqRecord if individual output files are requested
                if not self.combined:
                    with open(os.path.join(self.output_path, f'{gene_id}.fasta'), 'w') as mashed:
                        SeqIO.write(seq_record, mashed, 'fasta')
                # Otherwise append it to the list
                else:
                    allele_list.append(seq_record)
        if self.combined:
            with open(os.path.join(self.output_path, 'combined_alleles.fasta'), 'w') as combined:
                SeqIO.write(sequences=allele_list,
                            handle=combined,
                            format='fasta')

    def __init__(self, allele_file, combined):
        if allele_file.startswith('~'):
            self.allele_file = os.path.abspath(os.path.expanduser(os.path.join(allele_file)))
        else:
            self.allele_file = os.path.abspath(os.path.join(allele_file))
        assert os.path.isfile(self.allele_file), f'Cannot locate the supplied allele file: {self.allele_file}'
        self.combined = combined
        self.output_path = os.path.join(os.path.dirname(self.allele_file), 'outputs')
        make_path(self.output_path)
        self.record_dict = dict()
        self.accession_dict = dict()
        self.subunit_dict = dict()
        self.allele_dict = dict()


def cli():
    # Parser for arguments
    parser = ArgumentParser(description='Find stx subunits (A and B) with matching accessions. '
                                        'Combine the subunits into one ORF')
    parser.add_argument('-a', '--allele_file',
                        required=True,
                        help='Specify name and path of allele file')
    parser.add_argument('-c', '--combined',
                        action='store_true',
                        help='Create a single, combined output file rather than individual files. '
                             'Default is individual files')
    # Get the arguments into an object
    args = parser.parse_args()
    SetupLogging(debug=True)
    reformat = Reformat(allele_file=args.allele_file,
                        combined=args.combined)
    reformat.main()
    logging.info('Allele combining complete!')


if __name__ == '__main__':
    cli()
