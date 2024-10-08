#!/usr/bin/env python3
from olctools.accessoryFunctions.accessoryFunctions import make_path
from argparse import ArgumentParser
from urllib.parse import urlparse
import xml.dom.minidom as xml
import urllib.request as url
import ssl
import os
"""
Software Copyright License Agreement (BSD License)
The copyrights to the software code for SRST2 and accompanying code are licensed
under the following terms:
Copyright (c) 2013, Michael Inouye, Bernie Pope, Harriet Dashnow, Kathryn Holt
All rights reserved.
"""

'''
Download MLST datasets from this site: http://pubmlst.org/data/ by
parsing an xml file (http://pubmlst.org/data/dbases.xml).
Data is downloaded for a species determined by the user:
- profiles (maps STs to allele numbers)
- numbered sequences for each locus in the scheme
In addition, the alleles are concatenated together for use with SRST2.
A log file is also generated in the working directory, detailing the
time, date and location of all files downloaded, as well as the <retrieved> 
tag which tells us when the XML entry was last updated. 
If the species name input by the user matches multiple <species> in the
xml file, the script simply reports the possible matches so the user can
try again.
'''

ssl._create_default_https_context = ssl._create_unverified_context

def parse_args():
    parser = ArgumentParser(description='Download MLST datasets by species'
                                        'from pubmlst.org.')

    parser.add_argument('--repository_url',
                        metavar='URL',
                        default='http://pubmlst.org/data/dbases.xml',
                        help='URL for MLST repository XML index')

    parser.add_argument('--genus',
                        metavar='NAME',
                        required=True,
                        help='The name of the species that you want to download (e.g. "Escherichia")')

    parser.add_argument('--force_scheme_name',
                        action="store_true",
                        default=False,
                        help='Flag to force downloading of specific scheme name (e.g. "Clostridium difficile")')

    parser.add_argument('--path',
                        metavar='PATH',
                        default=os.getcwd(),
                        help='Path in which to store the downloaded alleles and profiles')

    return parser.parse_args()


# test if a node is an Element and that it has a specific tag name
def testelementtag(node, name):
    return node.nodeType == node.ELEMENT_NODE and node.localName == name


# Get the text from an element node with a text node child
def gettext(element):
    result = ''
    for node in element.childNodes:
        if node.nodeType == node.TEXT_NODE:
            result += node.data
    return normalisetext(result)


# remove unwanted whitespace including linebreaks etc.
def normalisetext(string):
    return ' '.join(string.split())


# A collection of interesting information about a taxa
class SpeciesInfo(object):
    def __init__(self):
        self.name = None  # String name of species
        self.database_url = None  # URL as string
        self.retrieved = None  # date as string
        self.profiles_url = None  # URL as string
        self.profiles_count = None  # positive integer
        self.loci = []  # list of loci


class LocusInfo(object):
    def __init__(self):
        self.url = None
        self.name = None


# retrieve the interesting information for a given sample element
def getspeciesinfo(species_node, species, exact):
    this_name = gettext(species_node)
    store = False
    if exact:
        if this_name == species:
            store = True
    else:
        if this_name.startswith(species):
            store = True
    if store:
        info = SpeciesInfo()
        info.name = this_name
        for mlst_node in species_node.getElementsByTagName('mlst'):
            for database_node in mlst_node.getElementsByTagName('database'):
                for database_child_node in database_node.childNodes:
                    if testelementtag(database_child_node, 'url'):
                        info.database_url = gettext(database_child_node)
                    elif testelementtag(database_child_node, 'retrieved') or testelementtag(database_child_node, b'retrieved'):
                        info.retrieved = gettext(database_child_node)
                    elif testelementtag(database_child_node, 'profiles'):
                        for profile_count in database_child_node.getElementsByTagName('count'):
                            info.profiles_count = gettext(profile_count)
                        for profile_url in database_child_node.getElementsByTagName('url'):
                            info.profiles_url = gettext(profile_url)
                    elif testelementtag(database_child_node, 'loci'):
                        for locus_node in database_child_node.getElementsByTagName('locus'):
                            locus_info = LocusInfo()
                            locus_info.name = gettext(locus_node)
                            for locus_url in locus_node.getElementsByTagName('url'):
                                locus_info.url = gettext(locus_url)
                            info.loci.append(locus_info)
        return info
    else:
        return None


def main(args):
    # Create the path to store the schemes (if necessary)
    make_path(args.path)
    # Allow for Shigella to use the Escherichia MLST profile/alleles
    args.genus = args.genus if args.genus != 'Shigella' else 'Escherichia'
    # As there are multiple profiles for certain organisms, this dictionary has the schemes I use as values
    organismdictionary = {'Escherichia': 'Escherichia coli#1',
                          'Vibrio': 'Vibrio parahaemolyticus',
                          'Campylobacter': 'Campylobacter jejuni',
                          'Listeria': 'Listeria monocytogenes',
                          'Bacillus': 'Bacillus cereus',
                          'Staphylococcus': "Staphylococcus aureus",
                          'Salmonella': 'Salmonella enterica',
                          'Yersinia': 'Yersinia ruckeri'}
    # Set the appropriate profile based on the dictionary key:value pairs
    try:
        args.genus = organismdictionary[args.genus]
    except (KeyError, AttributeError):
        pass
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE
    with url.urlopen(args.repository_url, context=ctx) as docfile:
        doc = xml.parse(docfile)
        root = doc.childNodes[0]
        found_species = []
        for species_node in root.getElementsByTagName('species'):
            info = getspeciesinfo(species_node, args.genus, args.force_scheme_name)
            if info is not None:
                found_species.append(info)
        if len(found_species) == 0:
            print("No species matched your query.")
            return
        if len(found_species) > 1:
            print("The following {} species match your query, please be more specific:".format(len(found_species)))
            for info in found_species:
                print(info.name)
            return
        # exit(2)
    # output information for the single matching species
    assert len(found_species) == 1
    species_info = found_species[0]
    species_name_underscores = species_info.name.replace(' ', '_')
    species_name_underscores = species_name_underscores.replace('/', '_')
    species_all_fasta_filename = species_name_underscores + '.fasta'
    species_all_fasta_file = open('{}/{}'.format(args.path, species_all_fasta_filename), 'w')
    log_filename = "mlst_data_download_{}.log".format(species_name_underscores)
    log_file = open('{}/{}'.format(args.path, log_filename), "w")
    #log_file.write(species_info.retrieved + '\n')
    profile_path = urlparse(species_info.profiles_url).path
    profile_filename = profile_path.split('/')[-1]
    log_file.write("definitions: {}\n".format(profile_filename))
    log_file.write("{} profiles\n".format(species_info.profiles_count))
    log_file.write("sourced from: {}\n\n".format(species_info.profiles_url))
    #
    # with url.urlopen(species_info.profiles_url) as profile_doc:
    #     with open(os.path.join(args.path, profile_filename), 'w') as profile_file:
    localfile, headers = url.urlretrieve(species_info.profiles_url)
    with open(localfile, 'r') as profile_doc:
        with open(os.path.join(args.path, profile_filename), 'w') as profile_file:
            profile_file.write(profile_doc.read())
    for locus in species_info.loci:
        locus_path = urlparse(locus.url).path
        locus_filename = locus_path.split('/')[-1]
        log_file.write("locus {}\n".format(locus.name))
        log_file.write(locus_filename + '\n')
        log_file.write("Sourced from {}\n\n".format(locus.url))
        #
        local_locus_doc, headers = url.urlretrieve(locus.url)
        with open(local_locus_doc, 'r') as locus_doc:
            with open(os.path.join(args.path, locus_filename), 'w') as locus_file:
                # locus_doc = url.urlopen(locus.url)
                # locus_file = open('{}/{}'.format(args.path, locus_filename), 'w')
                locus_fasta_content = locus_doc.read()
                locus_file.write(locus_fasta_content)
                species_all_fasta_file.write(locus_fasta_content)
                # locus_file.close()
                # locus_doc.close()
    log_file.write("all loci: {}\n".format(species_all_fasta_filename))
    log_file.close()
    species_all_fasta_file.close()


if __name__ == '__main__':
    arguments = parse_args()
    main(arguments)
