#!/usr/bin/env python3
from Bio import SeqIO
from glob import glob
import os

__author__ = 'adamkoziol'


class ResistanceNotes(object):

    @staticmethod
    def classes(targetpath):
        """
        Uses .tfa files included in the ResFinder database to determine the resistance class of gene matches
        :param targetpath: Path to database files
        :return: Dictionary of resistance class: gene set
        """
        # Initialise dictionary to store results
        resistance_dict = dict()
        # Find all the .tfa files in the folder
        resistance_files = sorted(glob(os.path.join(targetpath, '*.tfa')))
        # Iterate through each file
        for fasta in resistance_files:
            # Extract the resistance class from the file name and path
            resistance_class = os.path.splitext(os.path.basename(fasta))[0]
            # Initialise the resistance class as a set in the dictionary
            resistance_dict[resistance_class] = set()
            # Open the file
            with open(fasta) as resistance:
                # Iterate through the FASTA records
                for record in SeqIO.parse(resistance, 'fasta'):
                    # Replace dashes with underscores
                    record.id = record.id.replace('-', '_')
                    # Add the gene name to the set
                    resistance_dict[resistance_class].add(record.id)
        return resistance_dict

    @staticmethod
    def gene_name(name):
        """
        Split the FASTA header string into its components, including gene name, allele, and accession
        :param name: FASTA header
        :return: gname, genename, accession, allele: name of gene. Often the same as genename, but for certain entries
        it is longer, full gene name, accession, and allele extracted from the FASTA header
        """
        if 'Van' in name or 'mcr' in name or 'aph' in name or 'ddlA' in name or 'ant' in name or 'aadE_Cc' in name:
            try:
                if name == "ant(3'')_Ih_aac(6')_IId_1_AF453998":
                    # >aac(3)_Ib_aac(6')_Ib_1_AF355189 yields gname, genename: aac(3)-Ib-aac(6')-Ib, allele:1,
                    # accession: AF355189
                    gene1, version1, gene2, version2, allele, accession = name.split('_')
                    gname = '{g1}-{v1}-{g2}-{v2}'.format(g1=gene1,
                                                         v1=version1,
                                                         g2=gene2,
                                                         v2=version2)
                    genename = gname
                elif name == 'ant(3'')_Ia_1_X02340':
                    # >ant(3'')_Ia_1_X02340
                    gene, version, allele, accession = name.split('_')
                    gname = '{g}-{v}'.format(g=gene,
                                             v=version)
                    genename = gname
                elif 'mcr_3' in name or 'mcr_2' in name or 'mcr_1.10' in name:
                    # >mcr_3.3_1_NG055492 yields genename, gname: mcr-3, allele: 1, accession: NG055492
                    gene, combinedversion, allele, accession = name.split('_')
                    version = combinedversion.split('.')[0]
                    gname = '{gene}-{version}'.format(gene=gene,
                                                      version=version)
                    genename = gname
                else:
                    # Allow for an additional part to the gene name aph(3'')_Ib_5_AF321551 yields gname: aph(3''),
                    # genename: aph(3'')-Ib, allele: 5, accession AF321551
                    try:
                        pregene, postgene, allele, accession = name.split('_')
                        gname = '{pre}-{post}'.format(pre=pregene,
                                                      post=postgene)
                        genename = gname
                    except ValueError:
                        # Allow for underscores in the accession: aac(2')_Ie_1_NC_011896 yields gname: aac(2'),
                        # genename:  aac('2)-1e, allele: 1, accession NC_011896
                        pregene, postgene, allele, preaccession, postaccession = name.split('_')
                        genename = '{pre}-{post}'.format(pre=pregene,
                                                         post=postgene)
                        accession = '{pre}_{post}'.format(pre=preaccession,
                                                          post=postaccession)
                        gname = pregene
            except ValueError:
                # VanC_2_DQ022190
                genename, allele, accession = name.split('_')
                gname = genename
        else:
            if 'bla' in name or 'aac' in name or 'ARR' in name or 'POM' in name:
                if 'OKP' in name or 'CTX' in name or 'OXY' in name:
                    # >blaOKP_B_11_1_AM051161 yields gname: blaOKP-B-11, genename: blaOXP, allele: 1,
                    # accession: AM051161
                    gene, version1, version2, allele, accession = name.split('_')
                    gname = '{g}-{v1}-{v2}'.format(g=gene,
                                                   v1=version1,
                                                   v2=version2)
                    genename = gname
                elif 'CMY' in name:
                    # >blaCMY_12_1_Y16785 yields gname, genename: blaCMY, allele: 12
                    try:
                        gname, allele, version, accession = name.split('_')
                    except ValueError:
                        # blaCMY_59_1_NG_048854
                        gname, allele, version, pre_accession, post_accession = name.split('_')
                        accession = '{pre}_{post}'.format(pre=pre_accession,
                                                          post=post_accession)
                    genename = gname
                elif name == "aac(3)_Ib_aac(6')_Ib_1_AF355189":
                    # >aac(3)_Ib_aac(6')_Ib_1_AF355189 yields gname, genename: aac(3)-Ib-aac(6')-Ib, allele:1,
                    # accession: AF355189
                    gene1, version1, gene2, version2, allele, accession = name.split('_')
                    gname = '{g1}-{v1}-{g2}-{v2}'.format(g1=gene1,
                                                         v1=version1,
                                                         g2=gene2,
                                                         v2=version2)
                    genename = gname
                elif 'alias' in name:
                    # >blaSHV_5a_alias_blaSHV_9_1_S82452
                    gene1, version1, alias, gene2, version2, allele, accession = name.split('_')
                    gname = '{g1}-{v1}'.format(g1=gene1,
                                               v1=version1)
                    genename = gname
                else:
                    # Split the name on '_'s: ARR-2_1_HQ141279; gname, genename: ARR-2, allele: 1, accession: HQ141279
                    try:
                        genename, allele, accession = name.split('_')
                        gname = genename
                    except ValueError:
                        try:
                            # >blaACC_1_2_AM939420 yields gname: blaACC-1, genename: blaACC, allele: 2,
                            # accession: AM939420
                            genename, version, allele, accession = name.split('_')
                            gname = '{g}-{v}'.format(g=genename,
                                                     v=version)
                        except ValueError:
                            # >aac(2')_Ie_1_NC_011896 yields gname, genename: aac(2')-Ie, allele: 1,
                            # accession: NC_011896
                            genename, version, allele, preaccession, postaccession = name.split('_')
                            gname = '{g}-{v}'.format(g=genename,
                                                     v=version)
                            genename = gname
                            accession = '{preaccess}_{postaccess}'.format(preaccess=preaccession,
                                                                          postaccess=postaccession)
            else:
                # Split the name on '_'s: ARR-2_1_HQ141279; gname, genename: ARR-2, allele: 1, accession: HQ141279
                try:
                    genename, allele, accession = name.split('_')
                    gname = genename
                # Some names have a slightly different naming scheme:
                except ValueError:
                    # tet(44)_1_NZ_ABDU01000081 yields gname, genename: tet(44), allele: 1,
                    # accession: NZ_ABDU01000081
                    genename, allele, preaccession, postaccession = name.split('_')
                    accession = '{preaccess}_{postaccess}'.format(preaccess=preaccession,
                                                                  postaccess=postaccession)
                    gname = genename
        return gname, genename, accession, allele

    @staticmethod
    def resistance(genename, resistance_dict):
        """
        Determine the resistance class of the gene by searching the sets of genes included in every resistance FASTA
        file
        :param genename: Header string returned from analyses
        :param resistance_dict: Dictionary of resistance class: header
        :return: resistance class of the gene
        """
        # Initialise a list to store the resistance class(es) for the gene
        resistance_list = list()
        # Iterate through the dictionary of the resistance class: set of gene names
        for resistance_class, gene_set in resistance_dict.items():
            # If the gene is presence in the set
            if genename in gene_set:
                # Set the resistance class appropriately
                resistance_list.append(resistance_class)
        # Create a comma-separated string of the sorted genes
        resistance = ','.join(sorted(resistance_list))
        # Return the calculated resistance class
        return resistance
