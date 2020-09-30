#!/usr/bin/env python3
from olctools.accessoryFunctions.accessoryFunctions import make_path, run_subprocess, SetupLogging, write_to_logfile
from decimal import Decimal, ROUND_DOWN, ROUND_UP
from argparse import ArgumentParser
from threading import Thread
from queue import Queue
import multiprocessing
from glob import glob
from Bio import SeqIO
import logging
import random
import gzip
import csv
import os


class CreateRun(object):

    def main(self):
        self.file_list()
        self.populate_dictionary()
        self.art_illumina()
        # self.compress()
        if self.mixfile:
            self.read_mixfile()
            self.subsample_reads()
            self.mix_reads()

    def file_list(self):
        """
        Create a list of all the assemblies from which reads are to be simulated
        """
        logging.info('Locating assemblies in {path}'.format(path=self.assembly_path))
        # Use glob to find all the assemblies in the supplied folder
        self.assemblies = sorted(glob(os.path.join(self.assembly_path, '*.fasta')))
        # Ensure that there are assemblies in the folder before proceeding
        assert self.assemblies, 'Could not locate .fasta files in the supplied assembly path: {path}' \
            .format(path=self.assembly_path)

    def populate_dictionary(self):
        """
        Create a nested dictionary containing all the relevant input and output details for each assembly
        """
        logging.info('Calculating necessary values from assemblies')
        for i in range(self.threads):
            # Send the threads to
            threads = Thread(target=self.threaded_dictionary_populate, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        for assembly in sorted(self.assemblies):
            # Need to create outputs for both the empirical profile as well as the over-clustered profile
            for quality in self.qualities:
                make_path(self.profile_dictionary[quality]['output'])
                # Populate the dictionaries
                basename = os.path.splitext(os.path.basename(assembly))[0]
                self.basenames.append(basename)
                self.base_assemblies[basename] = assembly
                if assembly not in self.read_dictionary:
                    self.read_dictionary[assembly] = dict()
                    self.genome_size[basename] = dict()
                if quality not in self.read_dictionary[assembly]:
                    self.read_dictionary[assembly][quality] = dict()
                self.dict_queue.put((assembly, quality, basename))
        self.dict_queue.join()

    def threaded_dictionary_populate(self):
        while True:
            assembly, quality, basename = self.dict_queue.get()
            self.genome_size[basename] = self.genome_length(assembly=assembly)
            # Set the prefix of the reads to be :assemblyname_R e.g. 2015-SEQ-1364_R
            self.read_dictionary[assembly][quality]['prefix'] = \
                os.path.join(self.profile_dictionary[quality]['output'], '{basename}_R'
                             .format(basename=basename))
            # Add the base name of the assembly to the dictionary
            self.read_dictionary[assembly][quality]['basename'] = basename
            # Output files with have the direction as either '1' or '2'
            for direction in ['1', '2']:
                # Create the file names of the outputs from art_illumina (reads that are inputs for gzip), and
                # the outputs from gzip e.g. 2015-SEQ-1364_R1
                self.read_dictionary[assembly][quality][direction] = {
                    'input_reads': '{prefix}{direction}.fq'.format(
                        prefix=self.read_dictionary[assembly][quality]['prefix'],
                        direction=direction),
                    'output_reads': '{prefix}{direction}.fastq.gz'.format(
                        prefix=self.read_dictionary[assembly][quality]['prefix'],
                        direction=direction)
                }
            self.dict_queue.task_done()

    def art_illumina(self):
        """
        Use art_illumina to simulate the FASTQ reads with the supplied parameters
        """
        logging.info('Simulating FASTQ files with art_illumina')
        for i in range(self.threads):
            # Send the threads to
            threads = Thread(target=self.threaded_simulate, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        for assembly in sorted(self.read_dictionary):
            for quality, info_dict in self.read_dictionary[assembly].items():
                self.simulate_queue.put((quality, info_dict, assembly))
        self.simulate_queue.join()

    def threaded_simulate(self):
        while True:
            quality, info_dict, assembly = self.simulate_queue.get()
            # Create the system call to art_illumina and reformat.sh to simulate, and subsequently compress FASTQ files
            # -1: forward read profile file
            # -2: reverse read profile file
            # -i: input FASTA reference file
            # -l: read length of reads
            # -f: coverage of reads
            # -m: fragment length
            # -s: standard deviation of fragment length
            # -p: paired-end reads
            # -na: do not output ALN alignment file
            # -rs: seed for random number generator
            # -o: prefix of output filename
            # zl=9: Compression level set to 9 (maximum)
            # t=1: Cap the number of compression threads to one
            art_cmd = 'art_illumina -1 {forward_profile} -2 {reverse_profile} -i {assembly} -l {read_length} ' \
                      '-f {fcov} -m {frag_length} -s {frag_stdev} -p -na -rs {seed} -o {reads} ' \
                      '&& reformat.sh zl=9 in={in_for} in2={in_rev} out={out_for} out2={out_rev} t=1' \
                .format(forward_profile=self.profile_dictionary[quality]['forward'],
                        reverse_profile=self.profile_dictionary[quality]['reverse'],
                        read_length=self.read_length,
                        fcov=self.coverage,
                        frag_length=self.frag_length,
                        frag_stdev=self.frag_stdev,
                        assembly=assembly,
                        seed=self.random_seed,
                        reads=info_dict['prefix'],
                        in_for=info_dict['1']['input_reads'],
                        in_rev=info_dict['2']['input_reads'],
                        out_for=info_dict['1']['output_reads'],
                        out_rev=info_dict['2']['output_reads']
                        )
            # Only run the command if neither the art_illumina, or the compressed reads are present in the
            # output folder
            if not os.path.isfile(info_dict['2']['input_reads']) and not \
                    os.path.isfile(info_dict['2']['output_reads']):
                self.run_subprocess_write_log(cmd=art_cmd)
                for direction in ['1', '2']:
                    # Delete the uncompressed reads
                    if os.path.isfile(info_dict[direction]['output_reads']):
                        try:
                            os.remove(path=info_dict[direction]['input_reads'])
                        except FileNotFoundError:
                            pass
            self.simulate_queue.task_done()

    def read_mixfile(self):
        """
        Parse the supplied .csv file containing the FASTQ read mixing details
        """
        logging.info('Parsing supplied mix file')
        # Load the file into a dictionary using the csv library
        csv_dict = csv.DictReader(open(self.mixfile), delimiter=',')
        # Create sets and dictionaries to store any issues encountered with the mix file
        illegal_headers = set()
        missing_headers = set()
        missing_base_ids = set()
        missing_cont_ids = set()
        empirical_ids = set()
        duplicate_empirical = set()
        soc_ids = set()
        voc_ids = set()
        duplicate_soc = set()
        duplicate_voc = set()
        illegal_rows = dict()
        for i, entry in enumerate(csv_dict):
            # The row number will be stored for certain errors. It is the current enumerator + 2 (0-based vs 1-based,
            # and allowing for header)
            row_number = i + 2
            # Check to see if there are any illegal headers
            for header in entry:
                if header not in self.headers:
                    illegal_headers.add(header)
            # Check to see if any required headers are missing
            for header in self.headers:
                if header not in entry:
                    missing_headers.add(header)
            for quality in self.qualities:
                try:
                    # Check for base IDs missing corresponding assemblies
                    if entry['BaseID'] and entry['BaseID'] not in self.basenames:
                        missing_base_ids.add(entry['BaseID'])
                    # If a contaminant is specified, check to see if the assembly is missing
                    if entry['ContaminantID']:
                        if entry['ContaminantID'] not in self.basenames:
                            missing_cont_ids.add(entry['ContaminantID'])
                    if quality == 'empirical':
                        # Store the empirical sample ID
                        sample_id = entry['EID']
                        # Make sure that all the empirical IDs are unique
                        if entry['EID'] not in empirical_ids:
                            empirical_ids.add(entry['EID'])

                        else:
                            duplicate_empirical.add(entry['EID'])
                    elif quality == 'slightly-over-clustered':
                        # Store the slightly over-clustered sample ID
                        sample_id = entry['SOCID']
                        # Make sure that all the slightly over-clustered IDs are unique
                        if entry['SOCID'] not in soc_ids:
                            soc_ids.add(entry['SOCID'])
                        else:
                            duplicate_soc.add(entry['SOCID'])
                    else:
                        # Store the very over-clustered sample ID
                        sample_id = entry['VOCID']
                        # Make sure that all the over-clustered IDs are unique
                        if entry['OCID'] not in voc_ids:
                            voc_ids.add(entry['VOCID'])
                        else:
                            duplicate_voc.add(entry['VOCID'])
                    # Ensure that a base ID has been supplied
                    if not entry['BaseID']:
                        self.dict_append(dictionary=illegal_rows,
                                         row=row_number,
                                         string='Base ID missing')
                    # Validate the base fraction entries
                    if entry['BaseFraction']:
                        # Make sure that the supplied base fraction is a float, and is <= 1.0
                        try:
                            if float(entry['BaseFraction']) > 1.0:
                                self.dict_append(dictionary=illegal_rows,
                                                 row=row_number,
                                                 string='Base fraction ({bf}) is > 1.0'
                                                 .format(bf=entry['BaseFraction']))
                        except ValueError:
                            self.dict_append(dictionary=illegal_rows,
                                             row=row_number,
                                             string='The supplied base fraction, {bf}, is not a valid float'
                                             .format(bf=entry['BaseFraction']))
                    else:
                        # If the base fraction is not supplied, assume it is set to 1.0
                        if not entry['ContaminantFraction']:
                            entry['BaseFraction'] = '1.0'
                        # Ensure that there is no contamination fraction if there is a no base fraction
                        else:
                            self.dict_append(dictionary=illegal_rows,
                                             row=row_number,
                                             string='Cannot include a contaminant fraction without a base fraction')
                    # Validate the optional contaminant fraction entries
                    if entry['ContaminantFraction']:
                        # Ensure that a corresponding contaminant ID had been provided
                        if not entry['ContaminantID'] and float(entry['ContaminantFraction']) != 0:
                            self.dict_append(dictionary=illegal_rows,
                                             row=row_number,
                                             string='Contaminant fraction supplied, but no contaminant ID')
                        # Check to see that the contaminant fraction is a valid float, and that base fraction +
                        # contaminant fraction == 1.0
                        if entry['BaseFraction']:
                            try:
                                base_fraction = float(entry['BaseFraction'])
                            except ValueError:
                                base_fraction = 0
                            try:
                                if base_fraction + float(entry['ContaminantFraction']) != 1.0:
                                    self.dict_append(dictionary=illegal_rows,
                                                     row=row_number,
                                                     string='Base fraction ({bf}) + contaminant fraction ({cf}) = '
                                                            '{total}, not 1.0'
                                                     .format(bf=entry['BaseFraction'],
                                                             cf=entry['ContaminantFraction'],
                                                             total=base_fraction + float(
                                                                 entry['ContaminantFraction'])))
                            except ValueError:
                                self.dict_append(dictionary=illegal_rows,
                                                 row=row_number,
                                                 string='The supplied contaminant fraction, {bf}, is not a valid float'
                                                 .format(bf=entry['ContaminantFraction']))
                    else:
                        if entry['BaseFraction'] == '0':
                            entry['ContaminantFraction'] = '1.0'
                        else:
                            entry['ContaminantFraction'] = '0'
                    # Populate the dictionary
                    if quality not in self.mix_dict:
                        self.mix_dict[quality] = dict()
                    if sample_id not in self.mix_dict[quality]:
                        self.mix_dict[quality][sample_id] = dict()
                    self.mix_dict[quality][sample_id] = {
                        'base_ID': entry['BaseID'],
                        'base_fraction': float(entry['BaseFraction']),
                        'contaminant_id': entry['ContaminantID'],
                        'contaminant_fraction': float(entry['ContaminantFraction'])
                    }
                except KeyError:
                    pass
        # Create a list to store all the required fixes to the mix file
        quit_message = list()
        # Report all the errors
        if illegal_headers or missing_headers or illegal_rows or missing_base_ids or missing_cont_ids or \
                duplicate_empirical or duplicate_soc or duplicate_voc:
            logging.error('Errors in supplied mix file: {mf}'.format(mf=self.mixfile))
        if missing_base_ids or missing_cont_ids:
            logging.error('Possibly missing assemblies in folder {ap}'.format(ap=self.assembly_path))
            quit_message.append('Please add the required assemblies to {ap} or fix {mf}'
                                .format(ap=self.assembly_path,
                                        mf=self.mixfile))
        if missing_base_ids:
            logging.error('The following Base IDs are present in the mix file, but not in the assembly folder: {baseid}'
                          .format(baseid=', '.join(sorted(missing_base_ids))))
        if missing_cont_ids:
            logging.error('The following Contamination IDs are present in the mix file, but not in the '
                          'assembly folder: {contid}'
                          .format(contid=', '.join(sorted(missing_cont_ids))))
        if duplicate_empirical:
            logging.error('The following Empirical Sample IDs are present more than once in your mix file: {dupid}'
                          .format(dupid=', '.join(sorted(duplicate_empirical))))
            quit_message.append('Please remove the duplicate Empirical Sample IDs ({dupid}) from {mf}'
                                .format(dupid=', '.join(sorted(duplicate_empirical)),
                                        mf=self.mixfile))
        if duplicate_soc:
            logging.error('The following Slightly-Over-Clustered Sample IDs are present more than once in your mix '
                          'file: {dupid}'
                          .format(dupid=', '.join(sorted(duplicate_soc))))
            quit_message.append('Please remove the duplicate Slightly-Over-clustered Sample IDs ({dupid}) from {mf}'
                                .format(dupid=', '.join(sorted(duplicate_soc)),
                                        mf=self.mixfile))
        if duplicate_voc:
            logging.error('The following Very-Over-Clustered Sample IDs are present more than once in your mix file: '
                          '{dupid}'
                          .format(dupid=', '.join(sorted(duplicate_voc))))
            quit_message.append('Please remove the duplicate Very-Over-clustered Sample IDs ({dupid}) from {mf}'
                                .format(dupid=', '.join(sorted(duplicate_voc)),
                                        mf=self.mixfile))
        if missing_headers:
            logging.error('The following headers are missing from the mix file: {mh}'
                          .format(mh=', '.join(sorted(missing_headers))))
            quit_message.append('Please add the following headers to the mix file: {mh}'
                                .format(mh=', '.join(sorted(missing_headers))))
        if illegal_headers:
            logging.error('The following headers in the mix file are not valid: {ih}'
                          .format(ih=', '.join(sorted(illegal_headers))))
            quit_message.append('Please remove the following headers from the mix file: {ih}'
                                .format(ih=', '.join(sorted(illegal_headers))))
        if illegal_rows:
            for row in illegal_rows:
                logging.error('Row {rn} errors: {errs}'
                              .format(rn=row,
                                      errs=', '.join([entry for entry in illegal_rows[row]])))
            quit_message.append('Please see {lf} for list of issues with mix file entries'.format(lf=self.log))
        # Quit the script, and display the list of error messages
        if quit_message:
            quit('\n'.join(quit_message))

    def subsample_reads(self):
        """
        Subsample FASTQ files using reformat.sh
        """
        logging.info('Subsampling FASTQ reads with reformat.sh')
        for i in range(self.reformat_threads):
            # Send the threads to
            threads = Thread(target=self.threaded_subsample, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        for quality in self.qualities:
            count = 0
            for sample in sorted(self.mix_dict[quality]):
                count += 1
                self.subsample_queue.put((quality, sample, count))
        self.subsample_queue.join()

    def threaded_subsample(self):
        while True:
            quality, sample, count = self.subsample_queue.get()
            # Extract necessary variables from dictionaries
            sample_dict = self.mix_dict[quality][sample]
            base_id = sample_dict['base_ID']
            base_fraction = sample_dict['base_fraction']
            contaminant_id = sample_dict['contaminant_id']
            contaminant_fraction = sample_dict['contaminant_fraction']
            base_assembly = self.base_assemblies[base_id]
            # Set and create the path to store the mixed reads
            mix_path = os.path.join(self.output_path, quality, 'mix')
            make_path(mix_path)
            # Create the name and path of the final mixed FASTQ files
            sample_f = os.path.join(mix_path, '{sn}_S{count}_L001_R1_001.fastq.gz'
                                    .format(sn=sample,
                                            count=count))
            sample_r = os.path.join(mix_path, '{sn}_S{count}_L001_R2_001.fastq.gz'
                                    .format(sn=sample,
                                            count=count))
            # Initialise strings to store the name and path information of the subsampled FASTQ files
            base_f = str()
            base_r = str()
            cont_f = str()
            cont_r = str()
            # If the base fraction is 0, don't bother subsampling
            if base_fraction:
                # The number of total bases to subsample is defined as the total number of bases in the
                # reference genome * the base fraction of the final sample * the desired coverage level
                num_bases = int(self.genome_size[base_id] * base_fraction * self.coverage)
                # Update the name and path of the temporary files
                base_f = os.path.join(self.output_path, quality, '{sn}_base_R1.fastq.gz'.format(sn=sample))
                base_r = os.path.join(self.output_path, quality, '{sn}_base_R2.fastq.gz'.format(sn=sample))
                # Create the subsampling call to reformat.sh
                # sampleseed=seed for random number generator
                # samplebasestarget=number of output bases desired
                # t=1: Use one thread
                sample_cmd_base = 'reformat.sh in={base_in_for} in2={base_in_rev} out={base_out_for} ' \
                                  'out2={base_out_rev} sampleseed={seed} samplebasestarget={num_bases} t=1' \
                    .format(base_in_for=self.read_dictionary[base_assembly][quality]['1']['output_reads'],
                            base_in_rev=self.read_dictionary[base_assembly][quality]['2']['output_reads'],
                            base_out_for=base_f,
                            base_out_rev=base_r,
                            seed=self.random_seed,
                            num_bases=num_bases)
                # Only run the system call if the temporary subsampled files, as well as the final mixed files
                # do not exist
                if not os.path.isfile(base_f) and not os.path.isfile(base_r) and not os.path.isfile(sample_f) \
                        and not os.path.isfile(sample_r):
                    self.run_subprocess_write_log(cmd=sample_cmd_base)
            # Subsample the contaminant fraction as requested
            if contaminant_fraction:
                # Define the assembly file here - as this may be blank, there will be a KeyError if
                # contaminant_fraction does not exist
                cont_assembly = self.base_assemblies[contaminant_id]
                # Determine the number of bases to sample - same as above, but for the contaminant genome
                num_bases = int(self.genome_size[base_id] * contaminant_fraction * self.coverage)
                cont_f = os.path.join(self.output_path, quality, '{sn}_cont_R1.fastq.gz'.format(sn=sample))
                cont_r = os.path.join(self.output_path, quality, '{sn}_cont_R2.fastq.gz'.format(sn=sample))
                sample_cmd_cont = 'reformat.sh in={cont_in_for} in2={cont_in_rev} out={cont_out_for} ' \
                                  'out2={cont_out_rev} sampleseed={seed} samplebasestarget={num_bases}' \
                    .format(cont_in_for=self.read_dictionary[cont_assembly][quality]['1']['output_reads'],
                            cont_in_rev=self.read_dictionary[cont_assembly][quality]['2']['output_reads'],
                            cont_out_for=cont_f,
                            cont_out_rev=cont_r,
                            seed=self.random_seed,
                            num_bases=num_bases)
                # Create the subsampled files if they, and the final mixes do not exist
                if not os.path.isfile(cont_f) and not os.path.isfile(cont_r) and not os.path.isfile(sample_f) \
                        and not os.path.isfile(sample_r):
                    self.run_subprocess_write_log(cmd=sample_cmd_cont)
            # Populate the sample dictionary with the file path information
            self.sample_dict[sample] = {
                'sample_f': sample_f,
                'sample_r': sample_r,
                'base_f': base_f,
                'base_r': base_r,
                'cont_f': cont_f,
                'cont_r': cont_r
            }
            self.subsample_queue.task_done()

    def mix_reads(self):
        """
        Create mixes as defined in the mix file in a multithreaded manner
        """
        for i in range(self.reformat_threads):
            # Send the threads to
            threads = Thread(target=self.threaded_mix, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        logging.info('Mixing FASTQ reads')
        # Initialise variables used in renaming FASTQ reads
        instrument_id = 'GCRSR'
        control_num = 0
        lane = '1'
        # There are 38 tiles in an Illumina MiSeq run.
        tiles = [
            '1101', '1102', '1103', '1104', '1105', '1106', '1107', '1108', '1109', '1110', '1111', '1112', '1113',
            '1114', '1115', '1116', '1117', '1118', '1119',
            '2101', '2102', '2103', '2104', '2105', '2106', '2107', '2108', '2109', '2110', '2111', '2112', '2113',
            '2114', '2115', '2116', '2117', '2118', '2119'
        ]
        # Empirical and over-clustered runs will be considered consecutive runs. Set the run number and flowcell ID to
        # reflect this
        for quality in self.qualities:
            if quality == 'empirical':
                run_num = 1
                flowcell_id = 'GCRSR001'
            elif quality == 'very-over-clustered':
                run_num = 2
                flowcell_id = 'GCRSR002'
            else:
                run_num = 3
                flowcell_id = 'GCRSR003'
            # Initialise a count for setting the file names
            count = 1
            # Perform mixing for the requested qualities
            for sample in sorted(self.mix_dict[quality]):
                self.mix_queue.put((sample, quality, instrument_id, control_num, lane, tiles, run_num, flowcell_id,
                                    count))
                # Increment the count in preparation for the next sample
                count += 1
        self.mix_queue.join()

    def threaded_mix(self):
        while True:
            sample, quality, instrument_id, control_num, lane, tiles, run_num, flowcell_id, count = self.mix_queue.get()
            # Ensure that the output files do not exist before taking time to read in FASTQ files
            if not os.path.isfile(self.sample_dict[sample]['sample_f']) \
                    and not os.path.isfile(self.sample_dict[sample]['sample_r']):
                # https://wiki-bsse.ethz.ch/display/DBSSEQGF/Header+description+for+FASTQs+generated+with+
                # Illumina+bclToFastq+converter
                # Initialise lists to store all the FASTQ read records
                forward_read_list = list()
                reverse_read_list = list()
                # Load the reads from the subsampled FASTQ files associated with the sample
                if self.mix_dict[quality][sample]['base_fraction']:
                    # Open the compressed file with gzip
                    with gzip.open(self.sample_dict[sample]['base_f'], 'rt') as forward_base:
                        # Append to the list of all the FASTQ reads using SeqIO
                        forward_read_list += list(SeqIO.parse(forward_base, 'fastq'))
                    # Same for the reverse base reads
                    with gzip.open(self.sample_dict[sample]['base_r'], 'rt') as reverse_base:
                        reverse_read_list += list(SeqIO.parse(reverse_base, 'fastq'))
                # Contaminant files are processed in a similar manner
                if self.mix_dict[quality][sample]['contaminant_fraction']:
                    with gzip.open(self.sample_dict[sample]['cont_f'], 'rt') as forward_cont:
                        forward_read_list += list(SeqIO.parse(forward_cont, 'fastq'))
                    with gzip.open(self.sample_dict[sample]['cont_r'], 'rt') as reverse_cont:
                        reverse_read_list += list(SeqIO.parse(reverse_cont, 'fastq'))
                # Use random to mix the order of the read lists. Include the random seed variable to ensure
                # repeatability
                random.Random(self.random_seed).shuffle(forward_read_list)
                random.Random(self.random_seed).shuffle(reverse_read_list)
                # Calculate the number of reads per tile
                step_size = len(forward_read_list) / len(tiles)
                # If the total number of reads isn't divisible by the number of tiles, round the step size down
                # to the closest integer
                step_size = int(Decimal(step_size).to_integral_value(rounding=ROUND_DOWN))
                # Determine the extra reads that must be added to the final tile due to the rounding down
                remainder = len(forward_read_list) - (step_size * len(tiles))
                # Initialise a dictionary to store the y position for every read associated with that tile
                tile_dict = dict()
                # x and y range min and max calculated from assembly pipeline validation panel - extracted min and
                # max values from FASTQ files x: 1720 29720, y: 1017 29259'
                for tile in tiles:
                    # Add any extra reads to the final tile
                    if tile == '2119':
                        # Create a sorted list of randomly sampled values from the calculated y range. Pull the
                        # number of values corresponding to the number of reads associated with the tile
                        tile_dict[tile] = sorted(random.sample(range(1017, 29260), step_size + remainder))
                    else:
                        tile_dict[tile] = sorted(random.sample(range(1017, 29260), step_size))
                # If the number of reads is divisible by the number of tiles, the step_list will be one element
                # longer than it would be otherwise (e.g. 4 goes into 12 3 times, remainder 0, while 4 goes into
                # 11 twice, remainder 3). Add 1 to the total number of reads to ensure that the list is always
                # the same length
                if remainder != 0:
                    step_list = [i for i in range(0, len(forward_read_list), step_size)]
                else:
                    step_list = [i for i in range(0, len(forward_read_list) + 1, step_size)]
                # Remove the first element of the list it will be 0, and the code below expects it to be removed
                step_list.pop(0)
                # Set the final value as the total number of reads in the FASTQ file
                step_list[-1] = len(forward_read_list)
                # Create a dictionary with the current read position: iterator for use in determining which tile
                # index should be used in analyses
                step_dict = {step: i for i, step in enumerate(step_list)}
                # Use gzip to open the forward and reverse outputs
                with gzip.open(self.sample_dict[sample]['sample_f'], 'wt') as forward_out:
                    with gzip.open(self.sample_dict[sample]['sample_r'], 'wt') as reverse_out:
                        # Initialise a count to track the number of the read within the current tile
                        tile_read_count = 0
                        for read_num, record in enumerate(forward_read_list):
                            # Determine the current step based on the read number
                            step = self.calculate_step(read_number=read_num,
                                                       step_list=step_list)
                            # Use random to select a list of one x value from the calculated range
                            x = random.sample(range(1720, 29720), 1)[0]
                            # Use the current step to extract the index of the tile to use from the tile dictionary
                            tile = tiles[step_dict[step]]
                            # Use the tile and the tile-specific read number to extract the y position
                            y = tile_dict[tile][tile_read_count]
                            # Create the name following the Illumina FASTQ format
                            name = '{iid}:{rn}:{fcid}:{lane}:{tile}:{x_pos}:{y_pos}' \
                                .format(iid=instrument_id,
                                        rn=run_num,
                                        fcid=flowcell_id,
                                        lane=lane,
                                        tile=tile,
                                        x_pos=x,
                                        y_pos=y)
                            # Replace the read name and ID with this name
                            record.name = name
                            record.id = name
                            # The record description also contains the read orientation (1 forward, 2 for reverse),
                            # as well as the control number and index number
                            record.description = name + ' 1:N:{cn}:{index}' \
                                .format(cn=control_num,
                                        index=count)
                            # Rename the corresponding reverse read using the same logic
                            reverse_read_list[read_num].name = name
                            reverse_read_list[read_num].id = name
                            reverse_read_list[read_num].description = name + ' 2:N:{cn}:{index}' \
                                .format(cn=control_num,
                                        index=count)
                            # Write the renamed read to the appropriate output file
                            SeqIO.write(record, forward_out, 'fastq')
                            SeqIO.write(reverse_read_list[read_num], reverse_out, 'fastq')
                            # Determine whether the whether the tile read count needs to be reset for the next tile
                            tile_read_count = self.calculate_tile(read_number=read_num,
                                                                  step_list=step_list,
                                                                  tile_read_count=tile_read_count)
            # Clean up all the temporary files
            if os.path.isfile(self.sample_dict[sample]['sample_f']):
                try:
                    os.remove(self.sample_dict[sample]['base_f'])
                except FileNotFoundError:
                    pass
                try:
                    os.remove(self.sample_dict[sample]['cont_f'])
                except FileNotFoundError:
                    pass
            if os.path.isfile(self.sample_dict[sample]['sample_r']):
                try:
                    os.remove(self.sample_dict[sample]['base_r'])
                except FileNotFoundError:
                    pass
                try:
                    os.remove(self.sample_dict[sample]['cont_r'])
                except FileNotFoundError:
                    pass
            self.mix_queue.task_done()

    @staticmethod
    def calculate_step(read_number, step_list):
        """
        Determine to which step a read number corresponds
        """
        # Initialise an integer to store the index
        step = int()
        # Iterate through all the steps
        for step_value in step_list:
            # If the read number is less than the current step, assign the step value to be the current step
            if read_number <= step_value:
                step = step_value
                # As we are looking in a sorted list, if we continue, the step will always correspond to the last
                # item in the list
                break
        # Return the calculated step value
        return step

    @staticmethod
    def calculate_tile(read_number, step_list, tile_read_count):
        """
        Determine to which read index a read number corresponds
        """
        # Iterate through all the steps
        for step_value in step_list:
            # If the current read number is less than the value of the step (-1 due to 0-based indexing), allow the
            # tile_read_count to be incremented, as there are additional reads associated with the tile at this step
            if read_number < step_value - 1:
                tile_read_count += 1
                break
            # Otherwise the tile associated with the step has no additional reads, and the counter is reset
            elif read_number == step_value - 1:
                tile_read_count = 0
                break
        # Return the updated variable
        return tile_read_count

    def run_subprocess_write_log(self, cmd):
        """
        Run the subprocess, and write outputs to log file
        :param cmd: type STR: string of system call to run
        """
        # Run the subprocess
        out, err = run_subprocess(command=cmd)
        # Write the command, as well as stdout, and stderr to the log file
        write_to_logfile(out='{cmd}\n{out}'.format(cmd=cmd,
                                                   out=out),
                         err=err,
                         logfile=self.log)

    @staticmethod
    def genome_length(assembly):
        """
        Use SeqIO to determine the total genome length of the supplied assembly
        :param assembly: type STR: path to assembly file
        :return genome_length: INT of total calculated genome length
        """
        # Create a variable to store the total extracted genome length
        genome_length = 0
        # Iterate through all the contigs in the assembly, and add the length of each one to the total length
        for record in SeqIO.parse(assembly, 'fasta'):
            genome_length += len(record.seq)
        # Return the total genome length
        return genome_length

    @staticmethod
    def dict_append(dictionary, row, string):
        """
        Initialise (if necessary) the supplied dictionary with the row number: growing list of errors encountered
        :param dictionary: type DICT, dictionary containing illegal rows, and associated errors
        :param row: type INT, current row number of .csv file being parsed
        :param string: type STR, current issue encountered with row
        return updated dictionary
        """
        try:
            dictionary[row].append(string)
        except KeyError:
            dictionary[row] = [string]
        return dictionary

    def __init__(self, assembly_path, read_length, coverage, frag_length, stdev, random_seed, output, mixfile, profile,
                 threads, debug):
        if assembly_path.startswith('~'):
            self.assembly_path = os.path.abspath(os.path.expanduser(os.path.join(assembly_path)))
        else:
            self.assembly_path = os.path.abspath(os.path.join(assembly_path))
        assert os.path.isdir(self.assembly_path), 'Supplied assembly path location is not a valid directory {path}' \
            .format(path=self.assembly_path)
        self.read_length = read_length
        self.coverage = coverage
        self.frag_length = frag_length
        self.frag_stdev = stdev
        self.random_seed = random_seed
        if not output:
            output = 'outputs'
        self.output_path = os.path.join(self.assembly_path, output)
        self.log = os.path.join(self.assembly_path, 'log')
        try:
            os.remove(self.log)
        except FileNotFoundError:
            pass
        SetupLogging(debug=debug,
                     filehandle=self.log,
                     logfile_level=logging.ERROR)
        self.assemblies = list()
        self.read_dictionary = dict()
        self.sample_dict = dict()
        self.genome_size = dict()
        if profile == '0':
            self.qualities = ['empirical']
        elif profile == '1':
            self.qualities = ['slightly-over-clustered']
        elif profile == '2':
            self.qualities = ['very-over-clustered']
        else:
            self.qualities = ['empirical', 'slightly-over-clustered', 'very-over-clustered']
        self.profile_dictionary = {
            'empirical': {
                'forward': os.path.join(self.assembly_path, 'EmpMiSeq250R1.txt'),
                'reverse': os.path.join(self.assembly_path, 'EmpMiSeq250R2.txt'),
                'output': os.path.join(self.output_path, 'empirical')
            },
            'slightly-over-clustered': {
                'forward': os.path.join(self.assembly_path, '200811_profile_1.txt'),
                'reverse': os.path.join(self.assembly_path, '200811_profile_2.txt'),
                'output': os.path.join(self.output_path, 'slightly-over-clustered')
            },
            'very-over-clustered': {
                'forward': os.path.join(self.assembly_path, '181019_profile_1.txt'),
                'reverse': os.path.join(self.assembly_path, '181019_profile_2.txt'),
                'output': os.path.join(self.output_path, 'very-over-clustered')
            }
        }
        self.basenames = list()
        self.base_assemblies = dict()
        if mixfile:
            self.mixfile = os.path.join(self.assembly_path, mixfile)
            assert os.path.isfile(self.mixfile), 'Cannot find supplied file containing FASTQ mixing details {mf}' \
                .format(mf=mixfile)
        else:
            self.mixfile = str()
        self.mix_dict = dict()
        self.threads = threads
        # Reduce the number of threads for CPU/IO-intensive methods
        self.reformat_threads = int(Decimal(self.threads / 4).to_integral_value(rounding=ROUND_UP))
        self.dict_queue = Queue(maxsize=self.threads)
        self.simulate_queue = Queue(maxsize=self.threads)
        self.subsample_queue = Queue(maxsize=self.reformat_threads)
        self.mix_queue = Queue(maxsize=self.threads)
        self.headers = ['SampleNum', 'EID', 'SOCID', 'VOCID', 'BaseID', 'BaseFraction', 'ContaminantID',
                        'ContaminantFraction']
        logging.info('Welcome to the CFIA FASTQ read simulator')


def cli():
    # Parser for arguments
    parser = ArgumentParser(description='Create simulated FASTQ files from over-clustered runs, and from the empirical '
                                        'ART profile. Optionally mix reads to create contaminated samples.')
    parser.add_argument('-a', '--assemblypath',
                        required=True,
                        help='Path to folder containing assemblies, profiles, and mix file')
    parser.add_argument('-o', '--output',
                        type=str,
                        default=str(),
                        help='Name of folder in which outputs are to be created. Default is assemblypath / "outputs". '
                             'Note only the name can be changed, not the path')
    parser.add_argument('-m', '--mixfile',
                        type=str,
                        default=str(),
                        help='Name of .csv file containing details on the optional mixing of FASTQ outputs. '
                             'Column headers: EID,OCID,BaseID,BaseFraction,ContaminantID,ContaminantFraction, '
                             'where EID=empirical ID (e.g. 2020-GCRSR-0001), OCID=over-clustered ID '
                             '(e.g. 2020-GCRSR-0021), BaseID=SEQID of base genome (e.g. 2015-SEQ-1539), '
                             'BaseFraction=fraction of base genome to use in mix (if left blank, will defaul to 1.0, '
                             'and throw an error if a contaminant is included. Accepts any value between 0 and 1.0, '
                             'e.g. 0.99), ContaminantID=SEQID of contaminating genome (may be left blank), '
                             'ContaminantFraction=fraction of contaminant genome to use in mix (may be left blank. If '
                             'included BaseFraction + ContaminantFraction must equal 1.0)')
    parser.add_argument('-p', '--profiles',
                        choices=['0', '1', '2', '3'],
                        default='3',
                        help='Profiles to use for simulating reads. Acceptable options are 0: Empirical MiSeq, '
                             '1: Slightly over-clustered MiSeq, 2: Very-over clustered MiSeq, '
                             '3: Empirical AND both over-clustered MiSeq. Default is 3')
    parser.add_argument('-r', '--readlength',
                        type=float,
                        default=250,
                        help='Length of reads to simulate. Default is 250 bp')
    parser.add_argument('-c', '--coverage',
                        type=float,
                        default=50,
                        help='Average coverage depth of reads to simulate. Default in 50X')
    parser.add_argument('-f', '--fraglength',
                        type=float,
                        default=300,
                        help='Average fragment length used to create simulated reads. Default is 300 bp')
    parser.add_argument('-s', '--stdevfraglength',
                        type=float,
                        default=50,
                        help='Standard deviation of fragment lengths used to simulate reads. Default is 50 bp')
    parser.add_argument('-n', '--number',
                        type=int,
                        default=960,
                        help='Random number seed used to re-create identical datasets. Default is 960')
    parser.add_argument('-t', '--threads',
                        type=int,
                        default=multiprocessing.cpu_count() - 1,
                        help='Number of threads to use for analyses. Default is number of system CPUs - 1')
    parser.add_argument('-d', '--debug',
                        action='store_true',
                        help='Enable debug mode for the pipeline')
    # Get the arguments into an object
    arguments = parser.parse_args()
    # Some assertions to ensure that the supplied arguments are reasonable
    assert arguments.readlength < arguments.fraglength, 'Read length {rl} must be less than fragment length {fl}. ' \
        .format(rl=arguments.readlength,
                fl=arguments.fraglength)
    assert 36 <= arguments.readlength <= 250, 'Read length {rl} must be between 36 and 250 bp' \
        .format(rl=arguments.readlength)
    create_run = CreateRun(assembly_path=arguments.assemblypath,
                           read_length=str(arguments.readlength),
                           coverage=arguments.coverage,
                           frag_length=str(arguments.fraglength),
                           stdev=str(arguments.stdevfraglength),
                           random_seed=str(arguments.number),
                           output=arguments.output,
                           mixfile=arguments.mixfile,
                           profile=arguments.profiles,
                           threads=arguments.threads,
                           debug=arguments.debug)
    create_run.main()
    logging.info('Simulation and mixing complete!')


if __name__ == '__main__':
    cli()
