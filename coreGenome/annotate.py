#!/usr/bin/env python
from accessoryFunctions.accessoryFunctions import *
import spadespipeline.metadataprinter as metadataprinter
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from threading import Thread
import threading
from subprocess import call
from glob import glob
__author__ = 'adamkoziol'


class Annotate(object):

    def annotatethreads(self):
        """
        Perform multi-threaded prokka annotations of each strain
        """
        import spadespipeline.createobject as createobject
        # Move the files to subfolders and create objects
        self.runmetadata = createobject.ObjectCreation(self)
        # Fix headers
        self.headerthreads()
        printtime('Performing prokka analyses', self.start)
        # Create and start threads
        for i in range(self.cpus):
            # Send the threads to the appropriate destination function
            threads = Thread(target=self.annotate, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        for sample in self.runmetadata.samples:
            # Create the prokka attribute in the metadata object
            setattr(sample, 'prokka', GenObject())
            # docker run -v /path/to/sequences:/path/to/sequences coreGenome
            # prokka 2014-SEQ-0275.fasta --force --genus Escherichia --species coli --usegenus --addgenes
            # --prefix 2014-SEQ-0275 --locustag EC0275 --outputdir /path/to/sequences/2014-SEQ-0275/prokka
            sample.prokka.outputdir = os.path.join(sample.general.outputdirectory, 'prokka')
            # TODO Incorporate MASH/rMLST/user inputted genus, species results in the system call
            # Create the system call
            sample.prokka.command = 'docker run -v {}:{} {} ' \
                                    'prokka {} ' \
                                    '--force ' \
                                    '--genus {} ' \
                                    '--species {} ' \
                                    '--usegenus ' \
                                    '--addgenes ' \
                                    '--prefix {} ' \
                                    '--locustag {} ' \
                                    '--outdir {}' \
                .format(self.sequencepath, self.sequencepath, self.dockerimage, sample.general.fixedheaders,
                        self.genus, self.species, sample.name, sample.name, sample.prokka.outputdir)
            # sample.name.split('-')[-1]
            self.queue.put(sample)
        self.queue.join()
        # Create the core genome
        self.codingthreads()

    def annotate(self):
        while True:
            sample = self.queue.get()
            threadlock = threading.Lock()
            if not os.path.isfile('{}/{}.gff'.format(sample.prokka.outputdir, sample.name)):
                # call(sample.prokka.command, shell=True, stdout=self.devnull, stderr=self.devnull)
                out, err = run_subprocess(sample.prokka.command)
                threadlock.acquire()
                write_to_logfile(sample.prokka.command, sample.prokka.command, self.logfile)
                write_to_logfile(out, err, self.logfile)
                threadlock.release()
            # List of the file extensions created with a prokka analysis
            files = ['err', 'faa', 'ffn', 'fna', 'fsa', 'gbk', 'gff', 'log', 'sqn', 'tbl', 'txt']
            # List of the files created for the sample by prokka
            prokkafiles = glob('{}/*'.format(sample.prokka.outputdir))
            # Find out which files have been created in the analysis
            for extension in files:
                # If the file was created, set the file path/name as the data for the attribute
                if extension in [prokka.split('.')[1] for prokka in prokkafiles]:
                    for output in prokkafiles:
                        setattr(sample.prokka, output.split('.')[1], output)
                # Otherwise, populate the attribute with 'NA'
                else:
                    setattr(sample.prokka, extension, 'NA')
            self.queue.task_done()

    def headerthreads(self):
        """
        The contig ID must be twenty characters or fewer. The names of the headers created following SPAdes assembly
        are usually far too long. This renames them as the sample name
        """
        # Create and start threads
        for i in range(self.cpus):
            # Send the threads to the appropriate destination function
            threads = Thread(target=self.headers, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        # from Bio import SeqIO
        for sample in self.runmetadata.samples:
            # Create an attribute to store the path/file name of the fasta file with fixed headers
            sample.general.fixedheaders = sample.general.bestassemblyfile.replace('.fasta', '.ffn')
            self.headerqueue.put(sample)
        self.headerqueue.join()

    def headers(self):
        while True:
            sample = self.headerqueue.get()
            # A list of contigs with modified record.id values
            fixedheaders = list()
            # Only do this if the file with fixed headers hasn't previously been created
            if not os.path.isfile(sample.general.fixedheaders):
                # Refseq genomes don't necessarily have underscores (or contig numbers) in the headers
                count = 0
                formatcount = '{:04d}'.format(count)
                for record in SeqIO.parse(open(sample.general.bestassemblyfile, "rU"), "fasta"):
                    # Split off anything following the contig number
                    # >2013-SEQ-0129_1_length_303005_cov_13.1015_ID_17624 becomes
                    # >2013-SEQ-0129_1
                    record.id = record.id.split('_length')[0]
                    # Prokka has a requirement that the header is unique and less than or equal to 20 characters
                    if len(record.id) > 20:
                        # Extract the contig number from the string - assumption is that this number is the final
                        # entry in the string, and that there are underscore separating the different components
                        contignumber = record.id.split('_')[-1] if '_' in record.id else formatcount
                        # Subtract the length of the contig number (and an additional one for the underscore) from
                        # 20 for the string slice, and add the contig number at the end
                        record.id = record.id[:(20 - len(contignumber) - 1)] + '_{}'.format(formatcount)
                    # Clear the name and description attributes of the record
                    record.name = ''
                    record.description = ''
                    # Add this record to our list
                    fixedheaders.append(record)
                # Open the filtered assembly file
                with open(sample.general.fixedheaders, 'w') as formatted:
                    # Write the records in the list to the file
                    SeqIO.write(fixedheaders, formatted, 'fasta')
            self.headerqueue.task_done()

    def codingthreads(self):
        """
        Find CDS features in .gff files to filter out non-coding sequences from the analysis
        """
        printtime('Extracting CDS features', self.start)
        # Create and start threads
        for i in range(self.cpus):
            # Send the threads to the appropriate destination function
            threads = Thread(target=self.codingsequences, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        for sample in self.runmetadata.samples:
            self.codingqueue.put(sample)
        self.codingqueue.join()
        # Create CDS files and determine gene presence/absence
        self.corethreads()

    def codingsequences(self):
        while True:
            sample = self.codingqueue.get()
            with open(sample.prokka.gff, 'r') as gff:
                for feature in gff:
                    # Only interested in the sequence name if it is a CDS
                    if 'CDS' in feature:
                        # Extract the sequence name from the string. Example below
                        # 2013-SEQ-0123-2014_1	Prodigal:2.6	CDS	443	1741	.	+	0
                        # ID=0279_00002;Parent=0279_00002_gene;gene=kgtP_1;
                        # inference=ab initio prediction:Prodigal:2.6,similar to AA sequence:UniProtKB:P0AEX3;
                        # locus_tag=0279_00002;product=Alpha-ketoglutarate permease
                        name = feature.split('ID=')[1].split(';')[0]
                        # Add number and names of genes to dictionaries
                        try:
                            gene = feature.split('gene=')[1].split(';')[0]
                            # Remove duplicate genes e.g. aceE_1 and aceE_2, as the numbers seem to be added arbitrarily
                            if '_' not in gene:
                                try:
                                    self.genenames[name] = gene
                                    self.genes[gene] += 1
                                except KeyError:
                                    self.genes[gene] = 1
                                try:
                                    self.cdsset[sample.name].add(name)
                                except KeyError:
                                    self.cdsset[sample.name] = set()
                                    self.cdsset[sample.name].add(name)
                        except IndexError:
                            pass
            self.codingqueue.task_done()

    def corethreads(self):
        """
        Create a .cds file consisting of fasta records of CDS features for each strain
        """
        printtime('Creating CDS files and finding core genes', self.start)
        # Create and start threads
        for i in range(self.cpus):
            # Send the threads to the appropriate destination function
            threads = Thread(target=self.coregroups, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        for sample in self.runmetadata.samples:
            # Define the name of the file to store the CDS nucleotide sequences
            sample.prokka.cds = '{}/{}.cds'.format(sample.prokka.outputdir, sample.name)
            self.corequeue.put(sample)
        self.corequeue.join()
        # Write the core .fasta files for each gene
        self.corewriter()

    def coregroups(self):
        while True:
            sample = self.corequeue.get()
            if not os.path.isfile(sample.prokka.cds):
                with open(sample.prokka.cds, 'w') as cds:
                    # Use BioPython to iterate through the records in the .ffn file
                    for record in SeqIO.parse(open(sample.prokka.ffn, 'r'), 'fasta'):
                        # Extract records present in the CDS set
                        if record.id in self.cdsset[sample.name]:
                            # Write each record to the file containing the nucleotide CDSs
                            SeqIO.write(record, cds, 'fasta')
                            # Parse each record
                            self.cdsparse(record)
            else:
                for record in SeqIO.parse(open(sample.prokka.cds, 'r'), 'fasta'):
                    #
                    self.cdsparse(record)
            self.corequeue.task_done()

    def cdsparse(self, record):
        """
        Finds core genes, and records gene names and sequences in dictionaries
        :param record: SeqIO record
        """
        try:
            # Find genes that are present in all strains of interest - the number of times the gene is found is
            # equal to the number of strains. Earlier parsing ensures that the same gene is not present in a strain
            # more than once
            if self.genes[self.genenames[record.id]] == len(self.runmetadata.samples):
                # Add the gene names and sequences to the appropriate dictionaries
                try:
                    self.genesequence[self.genenames[record.id]].add(str(record.seq))
                # Initialise the dictionary as required, then populate as above
                except KeyError:
                    self.genesequence[self.genenames[record.id]] = set()
                    self.genesequence[self.genenames[record.id]].add(str(record.seq))
                try:
                    self.coresequence[str(record.seq)].add(record.id)
                except KeyError:
                    self.coresequence[str(record.seq)] = set()
                    self.coresequence[str(record.seq)].add(record.id)
        except KeyError:
            pass

    def corewriter(self):
        """
        Creates .fasta files containing all alleles for each gene
        """
        printtime('Creating core allele files', self.start)
        for gene in sorted(self.genesequence):
            self.geneset.add(gene)
            # Set the name of the allele file
            genefile = '{}/{}.fasta'.format(self.coregenelocation, gene)
            # If the file doesn't exist, create it
            if not os.path.isfile(genefile):
                with open(genefile, 'w') as core:
                    for count, sequence in enumerate(self.genesequence[gene]):
                        # The definition line is the gene name, and the allele number (count (+ 1 to compensate for
                        # base zero))
                        definitionline = '{}-{}'.format(gene, count + 1)
                        # Create a sequence record using BioPython
                        fasta = SeqRecord(Seq(sequence),
                                          # Without this, the header will be improperly formatted
                                          description='',
                                          # Use >:definitionline as the header
                                          id=definitionline)
                        # Use the SeqIO module to properly format the new sequence record
                        SeqIO.write(fasta, core, 'fasta')
                        for strain in self.coresequence[sequence]:
                            # Record the strain name, the gene name, and the allele number.
                            # [:-6] removes the contig number: 2014-SEQ-0276_00001 becomes 2014-SEQ-0276
                            try:
                                self.corealleles[strain[:-6]].update({gene: count + 1})
                            except KeyError:
                                self.corealleles[strain[:-6]] = {gene: count + 1}
            else:
                # If the file exists, don't recreate it; only iterate through the dictionary of gene sequences
                for count, sequence in enumerate(self.genesequence[gene]):
                    for strain in self.coresequence[sequence]:
                        # Populate the dictionary as above
                        try:
                            self.corealleles[strain[:-6]].update({gene: count + 1})
                        except KeyError:
                            self.corealleles[strain[:-6]] = {gene: count + 1}
        # Create a combined file of all the core genes to be used in typing strain(s) of interest
        if not os.path.isfile('{}/core_combined.fasta'.format(self.coregenelocation)):
            fastafiles = glob('{}/*.fasta'.format(self.coregenelocation))
            # Run the method for each allele
            self.combinealleles(fastafiles)
        # Run the profiler
        self.profiler()

    def combinealleles(self, alleles):
        """
        Creates a large multi-fasta file from all core genes in the analysis
        :param alleles: .fasta file for each core gene
        """
        printtime('Creating combined core allele file', self.start)
        if not os.path.isfile('{}/core_combined.tfa'.format(self.coregenelocation)):
            with open('{}/core_combined.tfa'.format(self.coregenelocation), 'w') as combinedfile:
                # Open each allele file
                for allele in sorted(alleles):
                    for record in SeqIO.parse(open(allele, "rU"), "fasta"):
                        # Extract the sequence record from each entry in the multifasta
                        # Remove and dashes or 'N's from the sequence data - makeblastdb can't handle sequences
                        # with gaps
                        # noinspection PyProtectedMember
                        record.seq._data = record.seq._data.replace('-', '').replace('N', '')
                        # Clear the name and description attributes of the record
                        record.name = ''
                        record.description = ''
                        # Write each record to the combined file
                        SeqIO.write(record, combinedfile, 'fasta')

    def profiler(self):
        """
        Calculates the core profile for each strain
        """
        printtime('Calculating core profiles', self.start)
        # Only create the profile if it doesn't exist already
        # if not os.path.isfile('{}/profile.txt'.format(self.profilelocation)):
        for strain in self.corealleles:
            # Add the gene name and allele number pair for each core gene in each strain
            self.coreset.add(tuple(sorted(self.corealleles[strain].items())))
        # Set the header to be similar to an MLST profile - ST,gene1,gene2,etc
        header = 'ST,{}\n'.format(','.join(sorted(self.geneset)))
        data = ''
        for count, core in sorted(enumerate(self.coreset)):
            # Increment count now to account for 0-based numbering
            count += 1
            # Add the sequence type number to the profile
            data += '{}'.format(count)
            # Store the sequence type for each strain
            for strain in self.corealleles:
                if tuple(sorted(self.corealleles[strain].items())) == core:
                    self.profiles[strain] = count
            # Add the allele number for each gene
            for gene in sorted(core):
                data += ',{}'.format(gene[1])
            data += '\n'
        # Write the profile
        with open('{}/profile.txt'.format(self.profilelocation), 'w') as profile:
            profile.write(header)
            profile.write(data)
        # Create a list of which strains correspond to the sequence types
        self.linker()

    def linker(self):
        """
        Link the sequence types to the strains. Create a .csv file of the linkages
        """
        import operator
        strainprofile = '{}/strainprofiles.txt'.format(self.profilelocation)
        if not os.path.isfile(strainprofile):
            header = 'Strain,SequenceType\n'
            data = ''
            # Sort the profiles based on sequence type
            sortedprofiles = sorted(self.profiles.items(), key=operator.itemgetter(1))
            # Associate the sequence type with each strain
            for strain, seqtype in sortedprofiles:
                for sample in self.runmetadata.samples:
                    if sample.name == strain:
                        sample.general.coretype = seqtype
                        data += '{},{}\n'.format(strain, seqtype)
            # Write the results to file
            with open(strainprofile, 'w') as profile:
                profile.write(header)
                profile.write(data)

    def __init__(self, inputobject):
        from queue import Queue
        self.path = inputobject.path
        self.sequencepath = inputobject.databasesequencepath
        self.start = inputobject.start
        self.cpus = inputobject.cpus
        self.genus = inputobject.genus
        self.species = inputobject.species
        self.runmetadata = MetadataObject()
        self.dockerimage = inputobject.dockerimage
        # Set and create necessary folders
        self.coregenelocation = os.path.join(self.path, 'coregenes', self.genus)
        self.profilelocation = os.path.join(self.path, 'profile', self.genus)
        make_path(self.profilelocation)
        make_path(self.coregenelocation)
        # Create class variables
        self.genes = dict()
        self.genenames = dict()
        self.genesequence = dict()
        self.cdsset = dict()
        self.coresequence = dict()
        self.geneset = set()
        self.corealleles = dict()
        self.coreset = set()
        self.profiles = dict()
        self.queue = Queue()
        self.corequeue = Queue()
        self.codingqueue = Queue()
        self.headerqueue = Queue()
        # self.devnull = open(os.devnull, 'wb')
        self.logfile = inputobject.logfile
        # Run the analyses
        self.annotatethreads()
        # Print the metadata to file
        metadataprinter.MetadataPrinter(self)
