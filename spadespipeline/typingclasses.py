#!/usr/bin/env python3
from accessoryFunctions.accessoryFunctions import combinetargets, filer, GenObject, MetadataObject, printtime, \
    make_path, run_subprocess, write_to_logfile
from accessoryFunctions.metadataprinter import MetadataPrinter
from spadespipeline.GeneSeekr import GeneSeekr
from sipprCommon.objectprep import Objectprep
from sipprCommon.sippingmethods import Sippr
from genesippr.genesippr import GeneSippr
from serosippr.serosippr import SeroSippr
from reporter.reports import Reports
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from csv import DictReader
from glob import glob
import xlsxwriter
import threading
import pandas
import shutil
import os
import re

__author__ = 'adamkoziol'


class GDCS(Sippr):

    def reporter(self):
        # Create the report object
        report = Reports(self)
        report.gdcsreporter()

    def __init__(self, inputobject):
        self.reports = str()
        self.samples = inputobject.runmetadata
        self.starttime = inputobject.starttime
        self.completemetadata = inputobject.runmetadata
        self.path = inputobject.path
        self.analysescomplete = True
        self.reportpath = inputobject.reportpath
        self.runmetadata = inputobject.runmetadata
        self.reporter()


class Plasmids(GeneSippr):

    def reporter(self):
        """
        Creates a report of the results
        """
        # Create the path in which the reports are stored
        make_path(self.reportpath)
        data = 'Strain,Gene,PercentIdentity,Length,FoldCoverage\n'
        with open(os.path.join(self.reportpath, self.analysistype + '.csv'), 'w') as report:
            for sample in self.runmetadata.samples:
                data += sample.name + ','
                try:
                    if sample[self.analysistype].results:
                        multiple = False
                        for name, identity in sample[self.analysistype].results.items():
                            if not multiple:
                                data += '{},{},{},{}\n'.format(name, identity,
                                                               len(sample[self.analysistype].sequences[name]),
                                                               sample[self.analysistype].avgdepth[name])
                            else:
                                data += ',{},{},{},{}\n'.format(name, identity,
                                                                len(sample[self.analysistype].sequences[name]),
                                                                sample[self.analysistype].avgdepth[name])
                            multiple = True
                    else:
                        data += '\n'
                except KeyError:
                    data += '\n'
            report.write(data)


class PlasmidExtractor(object):

    def main(self):
        """
        Run the methods in the correct order
        """
        self.run_plasmid_extractor()
        self.parse_report()

    def run_plasmid_extractor(self):
        """
        Create and run the plasmid extractor system call
        """
        printtime('Extracting plasmids', self.start)
        # Define the system call
        extract_command = 'PlasmidExtractor.py -i {inf} -o {outf} -p {plasdb} -d {db} -t {cpus} -nc'\
            .format(inf=self.path,
                    outf=self.plasmid_output,
                    plasdb=os.path.join(self.plasmid_db, 'plasmid_db.fasta'),
                    db=self.plasmid_db,
                    cpus=self.cpus)
        # Only attempt to extract plasmids if the report doesn't already exist
        if not os.path.isfile(self.plasmid_report):
            # Run the system calls
            out, err = run_subprocess(extract_command)
            # Acquire thread lock, and write the logs to file
            self.threadlock.acquire()
            write_to_logfile(extract_command, extract_command, self.logfile)
            write_to_logfile(out, err, self.logfile)
            self.threadlock.release()

    def parse_report(self):
        """
        Parse the plasmid extractor report, and populate metadata objects
        """
        printtime('Parsing Plasmid Extractor outputs', self.start)
        # A dictionary to store the parsed excel file in a more readable format
        nesteddictionary = dict()
        # Use pandas to read in the CSV file, and convert the pandas data frame to a dictionary (.to_dict())
        dictionary = pandas.read_csv(self.plasmid_report).to_dict()
        # Iterate through the dictionary - each header from the CSV file
        for header in dictionary:
            # Sample is the primary key, and value is the value of the cell for that primary key + header combination
            for sample, value in dictionary[header].items():
                # Update the dictionary with the new data
                try:
                    nesteddictionary[sample].update({header: value})
                # Create the nested dictionary if it hasn't been created yet
                except KeyError:
                    nesteddictionary[sample] = dict()
                    nesteddictionary[sample].update({header: value})
        # Get the results into the metadata object
        for sample in self.metadata:
            # Initialise the plasmid extractor genobject
            setattr(sample, self.analysistype, GenObject())
            # Initialise the list of all plasmids
            sample[self.analysistype].plasmids = list()
            # Iterate through the dictionary of results
            for line in nesteddictionary:
                # Extract the sample name from the dictionary in a manner consistent with the rest of the COWBAT
                # pipeline e.g. 2014-SEQ-0276_S2_L001 becomes 2014-SEQ-0276
                sample_name = nesteddictionary[line]['Sample']
                # Use the filer method to extract the name
                name = list(filer([sample_name]))[0]
                # Ensure that the names match
                if name == sample.name:
                    # Append the plasmid name extracted from the dictionary to the list of plasmids
                    sample[self.analysistype].plasmids.append(nesteddictionary[line]['Plasmid'])
        # Copy the report to the folder containing all reports for the pipeline
        try:
            shutil.copyfile(self.plasmid_report, os.path.join(self.reportpath, 'plasmidReport.csv'))
        except IOError:
            pass

    def __init__(self, inputobject):
        self.path = inputobject.path
        self.reportpath = inputobject.reportpath
        self.reffilepath = inputobject.reffilepath
        self.analysistype = 'plasmidextractor'
        self.plasmid_output = os.path.join(self.path, self.analysistype)
        self.plasmid_db = os.path.join(self.reffilepath, self.analysistype)
        self.plasmid_report = os.path.join(self.plasmid_output, 'plasmidReport.csv')
        self.start = inputobject.starttime
        self.cpus = inputobject.cpus
        self.logfile = inputobject.logfile
        self.threadlock = threading.Lock()
        self.metadata = inputobject.runmetadata.samples


class ResistanceNotes(object):

    @staticmethod
    def notes(targetpath):
        """
        Populates resistance dictionaries for different styles of gene:resistance entries
        :param targetpath: Directory in which the notes.txt file is located
        :return: the three resistance dictionaries
        """
        # Create a set of all the gene names without alleles or accessions e.g. sul1_18_AY260546 becomes sul1
        genedict = dict()
        altgenedict = dict()
        revaltgenedict = dict()
        # Load the notes file to a dictionary
        notefile = os.path.join(targetpath, 'notes.txt')
        with open(notefile, 'r') as notes:
            for line in notes:
                # Ignore comment lines - they will break the parsing
                if line.startswith('#'):
                    continue
                # Split the line on colons e.g. QnrB53:Quinolone resistance: has three variables after the split:
                # gene(QnrB53), resistance(Quinolone resistance), and _(\n) (unless there is an alternate name)
                gene, resistance, alternate = line.split(':')
                # Set up the resistance dictionary
                genedict[gene] = resistance
                # strA:Aminoglycoside resistance:Alternate name; aph(3'')-Ib - yields gene:resistance of aph(3'')-Ib
                # Aminoglycoside resistance
                if 'Alternate name' in line:
                    try:
                        altgene = alternate.split(';')[1].rstrip().lstrip()
                    except IndexError:
                        # blaGES-8:Beta-lactam resistance:Alternate name IBC-2
                        altgene = alternate.split()[-1].rstrip()
                    # Populate the dictionaries
                    genedict[altgene] = resistance
                    altgenedict[gene] = altgene
                    revaltgenedict[altgene] = gene
        return genedict, altgenedict, revaltgenedict

    @staticmethod
    def gene_name(name):
        """
        Split the FASTA header string into its components, including gene name, allele, and accession
        :param name: FASTA header
        :return:
        """
        # Allow for an additional part to the gene name aph(3'')_Ib_5_AF321551 yields gname: aph(3''), genename:
        # aph(3'')-Ib, allele: 5, accession AF321551
        if '_I' in name or 'Van' in name:
            pregene, postgene, allele, accession = name.split('_')
            genename = '{pre}-{post}'.format(pre=pregene,
                                             post=postgene)
            gname = pregene
        else:
            # Split the name on '_'s: ARR-2_1_HQ141279; gname, genename: ARR-2, allele: 1, accession: HQ141279
            try:
                genename, allele, accession = name.split('_')
                gname = genename
            # Some names have a slightly different naming scheme:
            except ValueError:
                try:
                    if 'bla' in name or 'aac' in name:
                        # >blaACC_1_2_AM939420 yields gname, genename: blaACC-1, allele: 2, accession: AM939420
                        genename, version, allele, accession = name.split('_')
                        gname = '{g}-{v}'.format(g=genename,
                                                 v=version)
                    else:
                        # tet(44)_1_NZ_ABDU01000081 yields gname, genename: tet(44), allele: 1,
                        # accession: NZ_ABDU01000081
                        genename, allele, preaccession, postaccession = name.split('_')
                        accession = '{preaccess}_{postaccess}'.format(preaccess=preaccession,
                                                                      postaccess=postaccession)
                        gname = genename
                except ValueError:
                    # Beta-lactamases have their own naming scheme
                    if name.split('_')[1].isdigit():
                        # blaOXY_1_1_1_Z30177 yields gname: blaOXY-1-1, genename: blaOXY, allele: 1, accession: Z30177
                        genename, version, allele, duplicate, accession = name.split('_')
                    else:
                        try:
                            # blaOKP_B_15_1_AM850917 yields gname: blaOKP-B-15, genename: blaOKP, allele: 1,
                            # accession: AM850917
                            genename, version, allele, unknown, accession = name.split('_')
                        except ValueError:
                            try:
                                # blaSHV_5a_alias_blaSHV_9_1_S82452 yields genename: blaSHV, version: 5a, alias: alias
                                # unknown: blaSHV, allele: 9, unknown: 1, accession: S82452
                                genename, version, alias, unknown, allele, unknown, accession = name.split('_')
                            except ValueError:
                                genename, version, alias, allele, unknown, accession = name.split('_')
                    gname = '{g}-{ver}-{a}'.format(g=genename,
                                                   ver=version,
                                                   a=allele)
        return gname, genename, accession, allele

    @staticmethod
    def resistance(gname, genename, genedict, altgenedict, revaltgenedict):
        """
        Extracts the resistance phenotype from the dictionaries using the gene name
        :param gname: Name of gene. Often the same as genename, but for certain entries it is longer
        e.g. blaOKP-B-15 instead of blaOKP
        :param genename: Name of gene e.g. blaOKP
        :param genedict: Dictionary of gene:resistance
        :param altgenedict: Dictionary of gene alternate name:resistance
        :param revaltgenedict: Dictionary of gene alternate name: gene
        :return: finalgene: gene name to be used in the report, the resistance phenotype
        """
        # If the gene name is present in the altgenedict dictionary, adjust the string to output
        # to include the alternate name in parentheses e.g. strA (aph(3'')-Ib
        try:
            finalgene = '{namegene} ({genealt})'.format(namegene=genename,
                                                        genealt=altgenedict[genename])
        except KeyError:
            # Similar to above except with revaltdict
            try:
                finalgene = '{namegene} ({genealt})'.format(namegene=revaltgenedict[genename],
                                                            genealt=genename)
            except KeyError:
                finalgene = genename
        # Extract the resistance from the genedict dictionary
        try:
            res = genedict[genename]
        except KeyError:
            try:
                res = genedict[gname]
            except KeyError:
                res = genedict[gname]
        return finalgene, res


class Serotype(SeroSippr):

    def runner(self):
        """
        Run the necessary methods in the correct order
        """
        printtime('Starting {} analysis pipeline'.format(self.analysistype), self.starttime)
        # Run the analyses
        ShortKSippingMethods(self, self.cutoff)
        printer = MetadataPrinter(self)
        printer.printmetadata()
        self.serotype_escherichia()
        self.serotype_salmonella()
        # Create the reports
        self.reporter()
        # Print the metadata
        printer = MetadataPrinter(self)
        printer.printmetadata()


class ShortKSippingMethods(Sippr):

    def main(self):
        """
        Run the methods in the correct order for pipelines
        """
        # Find the target files
        self.targets()
        # Use bbduk to bait the FASTQ reads matching the target sequences
        self.bait(maskmiddle='t', k=17)
        # If desired, use bbduk to bait the target sequences with the previously baited FASTQ files
        if self.revbait:
            self.reversebait(maskmiddle='t', k=17)
        # Run the bowtie2 read mapping module
        self.mapping()
        # Use samtools to index the sorted bam file
        self.indexing()
        # Parse the results
        self.parsing()
        # Clear out the large attributes that will difficult to handle objects
        self.clear()
        # Filter out any sequences with cigar features such as internal soft-clipping from the results
        self.clipper()


class ResSippr(GeneSippr):

    def runner(self):
        """
        Run the necessary methods in the correct order
        """
        printtime('Starting {} analysis pipeline'.format(self.analysistype), self.starttime)
        if not self.pipeline:
            general = None
            for sample in self.runmetadata.samples:
                general = getattr(sample, 'general')
            if general is None:
                # Create the objects to be used in the analyses
                objects = Objectprep(self)
                objects.objectprep()
                self.runmetadata = objects.samples
        # Run the analyses
        ShortKSippingMethods(self, self.cutoff)
        # Create the reports
        self.reporter()
        # Print the metadata
        printer = MetadataPrinter(self)
        printer.printmetadata()


class Resistance(ResSippr):

    def reporter(self):
        """
        Creates a report of the results
        """
        genedict, altgenedict, revaltgenedict = ResistanceNotes.notes(self.targetpath)
        # Find unique gene names with the highest percent identity
        for sample in self.runmetadata.samples:
            try:
                if sample[self.analysistype].results:
                    # Initialise a dictionary to store the unique genes, and their percent identities
                    sample[self.analysistype].uniquegenes = dict()
                    for name, identity in sample[self.analysistype].results.items():
                        # Split the name of the gene from the string e.g. ARR-2_1_HQ141279 yields ARR-2
                        genename = name.split('_')[0]
                        # Set the best observed percent identity for each unique gene
                        try:
                            # Pull the previous best identity from the dictionary
                            bestidentity = sample[self.analysistype].uniquegenes[genename]
                            # If the current identity is better than the old identity, save it
                            if float(identity) > float(bestidentity):
                                sample[self.analysistype].uniquegenes[genename] = float(identity)
                        # Initialise the dictionary if necessary
                        except KeyError:
                            sample[self.analysistype].uniquegenes[genename] = float(identity)
            except KeyError:
                pass
        # Create the path in which the reports are stored
        make_path(self.reportpath)
        # Initialise strings to store the results
        data = 'Strain,Resistance,Gene,Allele,Accession,PercentIdentity,Length,FoldCoverage\n'
        with open(os.path.join(self.reportpath, self.analysistype + '.csv'), 'w') as report:
            for sample in self.runmetadata.samples:
                data += sample.name + ','
                if sample[self.analysistype].results:
                    # Create an attribute to store the string for the eventual pipeline report
                    sample[self.analysistype].pipelineresults = list()
                    # If there are multiple results for a sample, don't write the name in each line of the report
                    multiple = False
                    for name, identity in sorted(sample[self.analysistype].results.items()):
                        # Extract the necessary variables from the gene name string
                        gname, genename, accession, allele = ResistanceNotes.gene_name(name)
                        # Retrieve the best identity for each gene
                        try:
                            percentid = sample[self.analysistype].uniquegenes[gname]
                        # Beta-lactamases will not have the allele and version from the gene name defined above
                        except KeyError:
                            percentid = sample[self.analysistype].uniquegenes[gname.split('-')[0]]
                        # If the percent identity of the current gene matches the best percent identity, add it to
                        # the report - there can be multiple occurrences of genes e.g.
                        # sul1,1,AY224185,100.00,840 and sul1,2,CP002151,100.00,927 are both included because they
                        # have the same 100% percent identity
                        if float(identity) == percentid:
                            try:
                                # Determine the name of the gene to use in the report, as well as its associated
                                # resistance phenotype
                                finalgene, res = ResistanceNotes.resistance(gname, genename, genedict, altgenedict,
                                                                            revaltgenedict)
                                # Treat the initial vs subsequent results for each sample slightly differently - instead
                                # of including the sample name, use an empty cell instead
                                if multiple:
                                    data += ','
                                # Populate the results
                                data += '{},{},{},{},{},{},{}\n'.format(
                                    res,
                                    finalgene,
                                    allele,
                                    accession,
                                    identity,
                                    len(sample[self.analysistype].sequences[name]),
                                    sample[self.analysistype].avgdepth[name])
                                sample[self.analysistype].pipelineresults.append(
                                    '{rgene} ({pid}%) {rclass}'.format(rgene=finalgene,
                                                                       pid=identity,
                                                                       rclass=res)
                                )
                                multiple = True
                            except KeyError:
                                pass
                else:
                    data += '\n'
            # Write the strings to the file
            report.write(data)


class ResFinder(GeneSeekr):

    @staticmethod
    def sequencenames(contigsfile):
        """
        Takes a multifasta file and returns a list of sequence names
        :param contigsfile: multifasta of all sequences
        :return: list of all sequence names
        """
        sequences = list()
        for record in SeqIO.parse(open(contigsfile, "rU", encoding="iso-8859-15"), "fasta"):
            sequences.append(record.id)
        return sequences

    def strainer(self):
        for sample in self.metadata:
            if sample.general.bestassemblyfile != 'NA':
                setattr(sample, self.analysistype, GenObject())
                targets = glob(os.path.join(self.targetpath, '*.tfa'))
                targetcheck = glob(os.path.join(self.targetpath, '*.tfa'))
                if targetcheck:
                    try:
                        combinedtargets = glob(os.path.join(self.targetpath, '*.fasta'))[0]
                    except IndexError:
                        combinetargets(targets, self.targetpath)
                        combinedtargets = glob(os.path.join(self.targetpath, '*.fasta'))[0]
                    sample[self.analysistype].targets = targets
                    sample[self.analysistype].combinedtargets = combinedtargets
                    sample[self.analysistype].targetpath = self.targetpath
                    sample[self.analysistype].targetnames = self.sequencenames(combinedtargets)
                    sample[self.analysistype].reportdir = os.path.join(sample.general.outputdirectory,
                                                                       self.analysistype)
                    make_path(sample[self.analysistype].reportdir)
                else:
                    # Set the metadata file appropriately
                    sample[self.analysistype].targets = 'NA'
                    sample[self.analysistype].combinedtargets = 'NA'
                    sample[self.analysistype].targetpath = 'NA'
                    sample[self.analysistype].targetnames = 'NA'
                    sample[self.analysistype].reportdir = 'NA'
                    sample[self.analysistype].blastresults = 'NA'
            else:
                # Set the metadata file appropriately
                setattr(sample, self.analysistype, GenObject())
                sample[self.analysistype].targets = 'NA'
                sample[self.analysistype].combinedtargets = 'NA'
                sample[self.analysistype].targetpath = 'NA'
                sample[self.analysistype].targetnames = 'NA'
                sample[self.analysistype].reportdir = 'NA'
                sample[self.analysistype].blastresults = 'NA'

    def resfinderreporter(self):
        """
        Custom reports for ResFinder analyses. These reports link the gene(s) found to their resistance phenotypes
        """
        # Initialise resistance dictionaries from the notes.txt file
        genedict, altgenedict, revaltgenedict = ResistanceNotes.notes(self.targetpath)
        # Create a workbook to store the report. Using xlsxwriter rather than a simple csv format, as I want to be
        # able to have appropriately sized, multi-line cells
        workbook = xlsxwriter.Workbook(os.path.join(self.reportpath, '{}.xlsx'.format(self.analysistype)))
        # New worksheet to store the data
        worksheet = workbook.add_worksheet()
        # Add a bold format for header cells. Using a monotype font size 10
        bold = workbook.add_format({'bold': True, 'font_name': 'Courier New', 'font_size': 8})
        # Format for data cells. Monotype, size 10, top vertically justified
        courier = workbook.add_format({'font_name': 'Courier New', 'font_size': 8})
        courier.set_align('top')
        # Initialise the position within the worksheet to be (0,0)
        row = 0
        col = 0
        # A dictionary to store the column widths for every header
        columnwidth = dict()
        extended = False
        headers = ['Strain', 'Gene', 'Resistance', 'PercentIdentity', 'PercentCovered', 'Contig', 'Location',
                   'nt_sequence']
        for sample in self.metadata:
            sample[self.analysistype].sampledata = list()
            # Process the sample only if the script could find targets
            if sample[self.analysistype].blastresults != 'NA':
                for result in sample[self.analysistype].blastresults:
                    # Set the name to avoid writing out the dictionary[key] multiple times
                    name = result['subject_id']
                    # Use the ResistanceNotes gene name extraction method to get the necessary variables
                    gname, genename, accession, allele = ResistanceNotes.gene_name(name)
                    # Initialise a list to store all the data for each strain
                    data = list()
                    # Determine the name of the gene to use in the report and the resistance using the resistance
                    # method
                    finalgene, resistance = ResistanceNotes.resistance(gname, genename, genedict, altgenedict,
                                                                       revaltgenedict)
                    # Append the necessary values to the data list
                    data.append(finalgene)
                    data.append(resistance)
                    percentid = result['percentidentity']
                    data.append(percentid)
                    data.append(result['alignment_fraction'])
                    data.append(result['query_id'])
                    data.append('...'.join([str(result['low']), str(result['high'])]))
                    try:
                        # Only if the alignment option is selected, for inexact results, add alignments
                        if self.align and percentid != 100.00:

                            # Align the protein (and nucleotide) sequences to the reference
                            self.alignprotein(sample, name)
                            if not extended:
                                # Add the appropriate headers
                                headers.extend(['aa_Identity',
                                                'aa_Alignment',
                                                'aa_SNP_location',
                                                'nt_Alignment',
                                                'nt_SNP_location'
                                                ])
                                extended = True
                            # Create a FASTA-formatted sequence output of the query sequence
                            record = SeqRecord(sample[self.analysistype].dnaseq[name],
                                               id='{}_{}'.format(sample.name, name),
                                               description='')

                            # Add the alignment, and the location of mismatches for both nucleotide and amino
                            # acid sequences
                            data.extend([record.format('fasta'),
                                         sample[self.analysistype].aaidentity[name],
                                         sample[self.analysistype].aaalign[name],
                                         sample[self.analysistype].aaindex[name],
                                         sample[self.analysistype].ntalign[name],
                                         sample[self.analysistype].ntindex[name]
                                         ])
                        else:
                            record = SeqRecord(Seq(result['subject_sequence'], IUPAC.unambiguous_dna),
                                               id='{}_{}'.format(sample.name, name),
                                               description='')
                            data.append(record.format('fasta'))
                            if self.align:
                                # Add '-'s for the empty results, as there are no alignments for exact matches
                                data.extend(['-', '-', '-', '-', '-'])
                    # If there are no blast results for the target, add a '-'
                    except (KeyError, TypeError):
                        data.append('-')
                    sample[self.analysistype].sampledata.append(data)

        if 'nt_sequence' not in headers:
            headers.append('nt_sequence')
        # Write the header to the spreadsheet
        for header in headers:
            worksheet.write(row, col, header, bold)
            # Set the column width based on the longest header
            try:
                columnwidth[col] = len(header) if len(header) > columnwidth[col] else columnwidth[
                    col]
            except KeyError:
                columnwidth[col] = len(header)
            worksheet.set_column(col, col, columnwidth[col])
            col += 1
        # Increment the row and reset the column to zero in preparation of writing results
        row += 1
        col = 0
        # Write out the data to the spreadsheet
        for sample in self.metadata:
            worksheet.write(row, col, sample.name, courier)
            columnwidth[col] = len(sample.name)
            worksheet.set_column(col, col, columnwidth[col])
            col += 1
            multiple = False
            if not sample[self.analysistype].sampledata:
                # Increment the row and reset the column to zero in preparation of writing results
                row += 1
                col = 0
                # Set the width of the row to be the number of lines (number of newline characters) * 12
                worksheet.set_row(row)
                worksheet.set_column(col, col, columnwidth[col])
            for data in sample[self.analysistype].sampledata:
                if multiple:
                    col += 1
                # List of the number of lines for each result
                totallines = list()
                for results in data:
                    #
                    worksheet.write(row, col, results, courier)
                    try:
                        # Counting the length of multi-line strings yields columns that are far too wide, only count
                        # the length of the string up to the first line break
                        alignmentcorrect = len(str(results).split('\n')[1])
                        # Count the number of lines for the data
                        lines = results.count('\n') if results.count('\n') >= 1 else 1
                        # Add the number of lines to the list
                        totallines.append(lines)
                    except IndexError:
                        try:
                            # Counting the length of multi-line strings yields columns that are far too wide, only count
                            # the length of the string up to the first line break
                            alignmentcorrect = len(str(results).split('\n')[0])
                            # Count the number of lines for the data
                            lines = results.count('\n') if results.count('\n') >= 1 else 1
                            # Add the number of lines to the list
                            totallines.append(lines)
                        # If there are no newline characters, set the width to the length of the string
                        except AttributeError:
                            alignmentcorrect = len(str(results))
                            lines = 1
                            # Add the number of lines to the list
                            totallines.append(lines)
                    # Increase the width of the current column, if necessary
                    try:
                        columnwidth[col] = alignmentcorrect if alignmentcorrect > columnwidth[col] else \
                            columnwidth[col]
                    except KeyError:
                        columnwidth[col] = alignmentcorrect
                    worksheet.set_column(col, col, columnwidth[col])
                    col += 1
                    multiple = True
                # Set the width of the row to be the number of lines (number of newline characters) * 12
                worksheet.set_row(row, max(totallines) * 11)
                # Increase the row counter for the next strain's data
                row += 1
                col = 0
        # Close the workbook
        workbook.close()

    def object_clean(self):
        """

        """
        for sample in self.metadata:
            try:
                delattr(sample[self.analysistype], 'aaidentity')
                delattr(sample[self.analysistype], 'aaalign')
                delattr(sample[self.analysistype], 'aaindex')
                delattr(sample[self.analysistype], 'ntalign')
                delattr(sample[self.analysistype], 'ntindex')
                delattr(sample[self.analysistype], 'dnaseq')
                delattr(sample[self.analysistype], 'blastresults')
            except KeyError:
                pass

    def __init__(self, inputobject):
        # qseqid sacc stitle positive mismatch gaps evalue bitscore slen length
        self.resfinderfields = ['query_id', 'subject_id', 'positives', 'mismatches', 'gaps', 'evalue', 'bit_score',
                                'subject_length', 'alignment_length', 'query_start', 'query_end']
        self.analysistype = 'resfinder_assembled'
        self.metadata = inputobject.runmetadata.samples
        self.cutoff = 70
        self.start = inputobject.starttime
        self.reportdir = inputobject.reportpath
        self.pipeline = True
        self.referencefilepath = inputobject.reffilepath
        self.targetpath = os.path.join(self.referencefilepath, 'resfinder')
        self.threads = inputobject.cpus
        self.align = True
        self.logfile = inputobject.logfile
        self.unique = True
        self.strainer()
        self.runmetadata = MetadataObject()
        self.runmetadata.samples = self.metadata
        GeneSeekr.__init__(self, self)
        self.resfinderreporter()
        self.object_clean()


class Prophages(GeneSeekr):

    def reporter(self):
        with open(os.path.join(self.reportpath, 'prophages.csv'), 'w') as report:
            data = 'Strain,Gene,Host,PercentIdentity,PercentCovered,Contig,Location\n'
            # Set the required variables to load prophage data from a summary file
            targetpath = os.path.join(self.referencefilepath, self.analysistype)
            overview = glob(os.path.join(targetpath, '*.txt'))[0]
            fieldnames = ['id_prophage', 'file_name', 'host', 'host_id', 'number_of_prophages_in_host',
                          'start_position_of_prophage', 'end_position_of_prophage', 'length_of_prophage']
            for sample in self.metadata:
                # Create a set to ensure that genes are only entered into the report once
                genes = set()
                if sample.general.bestassemblyfile != 'NA':
                    # Open the prophage file as a dict - I do this here, as if I open it earlier, it looks like the
                    # file remains partially-read through for the next iteration. Something like prophagedata.seek(0)
                    # would probably work, but Dictreader objects don't have a .seek attribute
                    prophagedata = DictReader(open(overview), fieldnames=fieldnames, dialect='excel-tab')
                    try:
                        if sample[self.analysistype].blastresults:
                            data += '{},'.format(sample.name)
                            # Allow for formatting multiple hits for the same sample
                            multiple = False
                            for result in sample[self.analysistype].blastresults:
                                gene = result['subject_id']
                                if gene not in genes:
                                    if multiple:
                                        data += ','
                                    # Iterate through the phage data in the dictionary
                                    for phage in prophagedata:
                                        if phage['id_prophage'] == gene:
                                            # Add the data to the row
                                            data += '{},{},{},{},{},{}..{}\n' \
                                                .format(gene,
                                                        phage['host'],
                                                        result['percentidentity'],
                                                        result['alignment_fraction'] if float(
                                                            result['alignment_fraction']) <= 100 else '100.0',
                                                        result['query_id'],
                                                        result['low'],
                                                        result['high'])
                                    genes.add(gene)
                                    # Set multiple to true for any additional hits for this sample
                                    multiple = True
                        else:
                            data += '{}\n'.format(sample.name)
                    except KeyError:
                        data += '{}\n'.format(sample.name)
                else:
                    data += '{}\n'.format(sample.name)
            report.write(data)


class Univec(GeneSeekr):

    def reporter(self):
        with open(os.path.join(self.reportpath, 'univec.csv'), 'w') as report:
            data = 'Strain,Gene,Description,PercentIdentity,PercentCovered,Contig,Location\n'
            for sample in self.metadata:
                if sample.general.bestassemblyfile != 'NA':
                    # Create a set to ensure that genes are only entered into the report once
                    genes = set()
                    try:
                        if sample[self.analysistype].blastresults:
                            data += '{},'.format(sample.name)
                            # If multiple hits are returned for a sample, don't re-add the sample name on the next row
                            multiple = False
                            for result in sample[self.analysistype].blastresults:
                                gene = result['subject_id']
                                # Parse the reference file in order to extract the description of the BLAST hits
                                for entry in SeqIO.parse(sample[self.analysistype].combinedtargets, 'fasta'):
                                    # Find the corresponding entry for the gene
                                    if entry.id == gene:
                                        # Cut out the description from the entry.description using regex
                                        # e.g. for 'gnl|uv|X66730.1:1-2687-49 B.bronchiseptica plasmid pBBR1 genes for
                                        # mobilization and replication' only save the string after '2687-49'
                                        description = re.findall('\d+-\d+\s(.+)', entry.description)[0]
                                        # Don't add the same gene more than once to the report
                                        if gene not in genes:
                                            if multiple:
                                                data += ','
                                            data += '{},{},{},{},{},{}..{}\n' \
                                                .format(gene.split('|')[-1],
                                                        description,
                                                        result['percentidentity'],
                                                        result['alignment_fraction'] if float(
                                                            result['alignment_fraction'])
                                                        <= 100 else '100.0',
                                                        result['query_id'],
                                                        result['low'],
                                                        result['high'])
                                            # Allow for the proper formatting
                                            multiple = True
                                            genes.add(gene)
                        else:
                            data += '{}\n'.format(sample.name)
                    except KeyError:
                        data += '{}\n'.format(sample.name)
                else:
                    data += '{}\n'.format(sample.name)
            report.write(data)


class Virulence(GeneSippr):

    def reporter(self):
        """
        Creates a report of the results
        """
        # Create a set of all the gene names without alleles or accessions e.g. sul1_18_AY260546 becomes sul1
        genedict = dict()
        # Load the notes file to a dictionary
        notefile = os.path.join(self.targetpath, 'notes.txt')
        with open(notefile, 'r') as notes:
            for line in notes:
                # Ignore comment lines - they will break the parsing
                if line.startswith('#'):
                    continue
                # Split the line on colons e.g. stx1Aa:  Shiga toxin 1, subunit A, variant a: has three variables after
                # the split: gene(stx1Aa), description(Shiga toxin 1, subunit A, variant a), and _(\n)
                try:
                    gene, description, _ = line.split(':')
                # There are exceptions to the parsing. Some lines only have one :, while others have three. Allow for
                # these possibilities.
                except ValueError:
                    try:
                        gene, description = line.split(':')
                    except ValueError:
                        gene, description, _, _ = line.split(':')
                # Set up the description dictionary
                genedict[gene] = description.replace(', ', '_').strip()
        # Find unique gene names with the highest percent identity
        for sample in self.runmetadata.samples:
            try:
                if sample[self.analysistype].results:
                    # Initialise a dictionary to store the unique genes, and their percent identities
                    sample[self.analysistype].uniquegenes = dict()
                    for name, identity in sample[self.analysistype].results.items():
                        # Split the name of the gene from the string e.g. stx1:11:Z36899:11 yields stx1
                        genename = name.split('_')[0]
                        # Only allow matches of 100% identity for stx genes
                        if 'stx' in genename and float(identity) < 100.0:
                            pass
                        else:
                            # Set the best observed percent identity for each unique gene
                            try:
                                # Pull the previous best identity from the dictionary
                                bestidentity = sample[self.analysistype].uniquegenes[genename]
                                # If the current identity is better than the old identity, save it
                                if float(identity) > float(bestidentity):
                                    sample[self.analysistype].uniquegenes[genename] = float(identity)
                            # Initialise the dictionary if necessary
                            except KeyError:
                                sample[self.analysistype].uniquegenes[genename] = float(identity)
            except KeyError:
                pass
        # Create the path in which the reports are stored
        make_path(self.reportpath)
        # Initialise strings to store the results
        data = 'Strain,Gene,Subtype/Allele,Description,Accession,PercentIdentity,FoldCoverage\n'
        with open(os.path.join(self.reportpath, self.analysistype + '.csv'), 'w') as report:
            for sample in self.runmetadata.samples:
                data += sample.name + ','
                try:
                    if sample[self.analysistype].results:
                        # If there are many results for a sample, don't write the sample name in each line of the report
                        multiple = False
                        for name, identity in sorted(sample[self.analysistype].results.items()):
                            try:
                                # Split the name on colons: stx2A:63:AF500190:d; gene: stx2A, allele: 63, accession:
                                # AF500190, subtype: d
                                genename, allele, accession, subtype = name.split('_')
                            # Treat samples without a subtype e.g. icaC:intercellular adhesion protein C: differently.
                            # Extract the allele as the 'subtype', and the gene name, and accession as above
                            except ValueError:
                                genename, subtype, accession = name.split('_')
                            # Retrieve the best identity for each gene
                            percentid = sample[self.analysistype].uniquegenes[genename]
                            # If the percent identity of the current gene matches the best percent identity, add it to
                            # the report - there can be multiple occurrences of genes e.g.
                            # sul1,1,AY224185,100.00,840 and sul1,2,CP002151,100.00,927 are both included because they
                            # have the same 100% percent identity
                            if float(identity) == percentid:
                                # Treat the initial vs subsequent results for each sample slightly differently - instead
                                # of including the sample name, use an empty cell instead
                                if multiple:
                                    data += ','
                                try:
                                    description = genedict[genename]
                                except KeyError:
                                    description = 'na'
                                # Populate the results
                                data += '{},{},{},{},{},{}\n'.format(
                                    genename,
                                    subtype,
                                    description,
                                    accession,
                                    identity,
                                    sample[self.analysistype].avgdepth[name])
                                multiple = True
                    else:
                        data += '\n'
                except KeyError:
                    data += '\n'
            # Write the strings to the file
            report.write(data)
