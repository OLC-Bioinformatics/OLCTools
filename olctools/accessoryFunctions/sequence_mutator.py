#!/usr/bin/env python
from olctools.accessoryFunctions.accessoryFunctions import make_path, SetupLogging
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from argparse import ArgumentParser, RawTextHelpFormatter
import itertools
import logging
import random
import math
import sys
import os


def read_target_file(target_file):
    """
    Read in the target file using SeqIO
    :param target_file: type str: Name and path of file containing target sequence to mutate
    :return: records: First entry in target file in Bio.SeqRecord.SeqRecord format
    """
    logging.info('Extracting target sequence')
    # Create a list to store the record - a generator will be consumed, and we need to access the record more than once
    records = list()
    # Read in the record in the target file
    for record in SeqIO.parse(handle=target_file, format='fasta'):
        # Append the record to the list
        records.append(record)
    # Return the first (and only) record in the list
    logging.debug(f'Target sequence is {records[0].seq}')
    return records[0]


def mutate_bases(location, record, random_seed, mutation_dict, mutated_seq, mutated_header, adjust=True):
    """
    Mutate the base of the target sequence at the requested location
    :param location: type int: Index of the base in the target sequence to be mutated
    :param record: type: Bio.SeqRecord.SeqRecord object created from the target sequence by SeqIO
    :param random_seed: type int: User-provided seed value to allow consistent outputs
    :param mutation_dict: type dict: Dictionary storing the list of potential replacement bases for each nucleotide
    :param mutated_seq: type str: Target sequence to be modified with the mutated bases
    :param mutated_header: type header: Header of the target sequence to be modified with details of the mutations
    :param adjust: type bool: Boolean whether the location needs to be adjusted to base-1. Default is True
    :return: mutated_seq: Target sequence modified with the mutations
    :return: mutated_header: Target header modified with details of the mutations
    """
    # Adjust the location to be consistent with base-1 as required
    if adjust:
        adjusted_location = int(location) + 1
    else:
        adjusted_location = location
    # If a random seed wasn't provided, the choice will not be repeatable. Otherwise, adjust the random seed,
    # so that there is some randomness (but is still repeatable) based on the position, name, and length of
    # the record
    if random_seed:
        # Set the random seed by adding the user-provided seed, the length of the record, and the length of the
        # record header to the adjusted location to create a target, and position-specific repeatable change
        adjusted_random_seed = random_seed + adjusted_location + len(record.seq) + len(record.id)
        logging.debug(f'Using random seed: {adjusted_random_seed} for location {location} in {record.id}')
        random.seed(adjusted_random_seed)
    # Extract the sequence of the current base
    current_base = record.seq[adjusted_location]
    # Use random.choice (this uses the random.seed defined above) to choose an entry from the list of bases
    mutated_base = random.choice(mutation_dict[current_base])
    logging.debug(f'Mutating {record.id} base {location} from {current_base} to {mutated_base}')
    # Create a string of the mutated sequence by using slices of the sequence up to the position + the mutated
    # base + the remainder of the string
    mutated_seq = mutated_seq[:adjusted_location] + mutated_base + mutated_seq[adjusted_location + 1:]
    # Add the current location, the original base, and the mutated base to the header
    mutated_header += f'_base_{location}_{current_base}_to_{mutated_base}'
    return mutated_seq, mutated_header


def mutate_target_locations(record, locations, mutation_dict, random_seed, adjust=True):
    """
    Mutate the target sequence at the requested location(s)
    :param record: type: Bio.SeqRecord.SeqRecord object created from the target sequence by SeqIO
    :param locations: type list: List of locations to mutate
    :param mutation_dict: type dict: Dictionary storing the list of potential replacement bases for each nucleotide
    :param random_seed: type int: User-provided seed value to allow consistent outputs
    :param adjust: type bool: Boolean whether the location needs to be adjusted to base-1. Default is True
    """
    logging.info('Mutating requested base(s)')
    # Initialise the string of the target sequence, so it can be iteratively mutated at each desired location
    mutated_seq = str(record.seq)
    # Initialise the header of the output sequence, so it can be iteratively updated after each mutation
    mutated_header = record.id
    # Update the header and sequence for each mutated location
    for location in locations:
        mutated_seq, mutated_header = \
            mutate_bases(
                location=location,
                record=record,
                random_seed=random_seed,
                mutation_dict=mutation_dict,
                mutated_seq=mutated_seq,
                mutated_header=mutated_header,
                adjust=adjust)
    # Create a SeqRecord object using a Seq object of the sequence and the updated header
    mutated_record = SeqRecord(seq=Seq(mutated_seq),
                               description=str(),
                               id=mutated_header)
    logging.debug(f'Mutated target sequence is {str(mutated_record.seq)}')
    return mutated_record


def create_mutation_location_list(record, number, random_seed):
    """
    Create a list of mutations to locate in the target sequence
    :param record: type: Bio.SeqRecord.SeqRecord object created from the target sequence by SeqIO
    :param number: type int: User-provided number of mutations to introduce
    :param random_seed: type int: User-provided seed value to allow consistent outputs
    :return: locations: Sorted list of bases to mutate in the target sequence
    """
    # Find all target locations by creating a list of the range of the length of the sequence
    target_locations = list(range(len(record.seq)))
    # If a random seed wasn't provided, the choice will not be repeatable. Otherwise, adjust the random seed,
    # so that there is some randomness (but is still repeatable) the name and length of the record
    if random_seed:
        # Set the random seed by adding the user-provided seed, the length of the record, and the length of the
        # record header to the adjusted location to create a target, and position-specific repeatable change
        adjusted_random_seed = random_seed + len(record.seq) + len(record.id)
        logging.debug(f'Using random seed: {adjusted_random_seed} for {record.id}')
        random.seed(adjusted_random_seed)
    locations = sorted(random.sample(target_locations, number))
    logging.debug(f'The following location(s) will be mutated in {record.id}: {", ".join(str(i) for i in locations)}')
    return locations



def write_outputs(output_file, mutated_record):
    """
    Use SeqIO to write the updated sequence(s) to file
    :param output_file: type str: Name and path of file in which the mutated sequence is to be saved
    :param mutated_record: type: Bio.SeqRecord.SeqRecord object created from the mutated sequence by
    Bio.SeqRecord.SeqRecord
    """
    logging.info(f'Outputs written to: {output_file}')
    # Write the mutated record to file
    with open(output_file, 'w') as output:
        SeqIO.write(
            sequences=mutated_record,
            handle=output,
            format='fasta'
        )


class SequenceMutatorLocation(object):

    def main(self):
        # Parse the target sequence with SeqIO
        record = read_target_file(target_file=self.target_file)
        # Ensure that the provided locations are valid
        self.locations = self.assert_location_in_target(
            record=record,
            locations=self.locations
        )
        # Create the mutated SeqRecord from the target sequence
        mutated_record = mutate_target_locations(
            record=record,
            locations=self.locations,
            mutation_dict=self.mutation_dict,
            random_seed=self.random_seed,
        )
        # Write the mutated record to the output file
        write_outputs(output_file=self.output_file,
                      mutated_record=mutated_record)

    @staticmethod
    def assert_location_in_target(record, locations):
        """
        Ensure that the provided locations are numerals, and fall within the length of the target sequence
        :param record: type: Bio.SeqRecord.SeqRecord object created from the target sequence by SeqIO
        :param locations: type list: User-provided locations to mutate
        :return: cleaned_locations: List of locations that pass the checks
        """
        logging.info('Ensuring that all requested locations are valid')
        # Initialise a list to store any missing positions
        missing = list()
        # Initialise a string to store the name of the target
        seq_name = str()
        # Initialise a list to store the curated locations
        cleaned_locations = list()
        for location in locations:
            # If the user provides a list with a hanging comma e.g. 5,67, the split of the provided string to create
            # the list will have an empty value. Ignore this, and remove it from the curated locations
            if location:
                # Ensure that all the locations provided are numbers
                if not location.isnumeric():
                    # If a location is not a number, inform the user of the error and exit
                    logging.error(f'{location} is not a valid number. Please ensure that you entered all information '
                                  f'correctly')
                    raise SystemExit
                # Store the name of the target sequence
                seq_name = record.name
                # Check to see if the location is present in the (base-1 adjusted length of the target sequence)
                if int(location) not in [i + 1 for i in range(len(record.seq))]:
                    missing.append(location)
                else:
                    cleaned_locations.append(location)
        # If there are locations that do not fall within the length of the target sequence, inform the user, and exit
        if missing:
            logging.error(f'Missing the following position(s) in {seq_name}: {", ".join(missing)}')
            raise SystemExit
        return cleaned_locations

    def __init__(self, target_file, locations, output_path, random_seed):
        # Expand the path of the target file as necessary
        if target_file.startswith('~'):
            self.target_file = os.path.abspath(os.path.expanduser(os.path.join(target_file)))
        else:
            self.target_file = os.path.abspath(os.path.join(target_file))
        try:
            assert os.path.isfile(self.target_file), \
                f'Cannot located target file: {self.target_file}. ' \
                f'Please ensure that you entered the correct name and path'
        except AssertionError:
            raise SystemExit
        if not output_path:
            self.output_path = os.path.join(os.path.dirname(self.target_file), 'outputs')
        else:
            if output_path.startswith('~'):
                self.output_path = os.path.abspath(os.path.expanduser(os.path.join(output_path)))
            else:
                self.output_path = os.path.abspath(os.path.join(output_path))
        make_path(self.output_path)
        self.output_file = os.path.join(self.output_path, 'mutated_target.fasta')
        self.locations = locations.split(',')
        self.random_seed = int(random_seed)
        self.mutation_dict = {
            'A': ['C', 'G', 'T'],
            'C': ['A', 'G', 'T'],
            'G': ['A', 'C', 'T'],
            'T': ['A', 'C', 'G']
        }


class SequenceMutatorRange(object):

    def main(self):
        # Parse the target sequence with SeqIO
        record = read_target_file(target_file=self.target_file)
        # Ensure that the proper method is used depending on whether spacing is required
        if not self.spacing:
            # If the minimum and maximum are the same, only a single output is required
            if self.minimum == self.maximum:
                mutated_record = self.mutate_target(
                    record=record,
                    minimum=self.minimum,
                    mutation_dict=self.mutation_dict,
                    random_seed=self.random_seed,
                )
                # Write the mutated record to the output file
                write_outputs(output_file=self.output_file,
                              mutated_record=mutated_record)
            else:
                self.mutate_target_range(
                    record=record,
                    minimum=self.minimum,
                    maximum=self.maximum,
                    mutation_dict=self.mutation_dict,
                    random_seed=self.random_seed,
                    output_file=self.output_file
                )
        else:
            self.mutate_target_range_spacing(
                record=record,
                minimum=self.minimum,
                maximum=self.maximum,
                spacing=self.spacing,
                mutation_dict=self.mutation_dict,
                random_seed=self.random_seed,
                output_file=self.output_file
            )

    @staticmethod
    def mutate_target(record, minimum, mutation_dict, random_seed):
        """
        Mutate the target sequence with the desired number of mutations
        :param record: type: Bio.SeqRecord.SeqRecord object created from the target sequence by SeqIO
        :param minimum: type int: Desired number of mutations to introduce into the target sequence
        :param mutation_dict: type dict: Dictionary storing the list of potential replacement bases for each nucleotide
        :param random_seed: type int: User-provided seed value to allow consistent outputs
        :return: mutated_record: Bio.SeqRecord.SeqRecord object created from the mutated sequence by
        Bio.SeqRecord.SeqRecord
        """
        logging.info(f'Mutating {minimum} base(s) in target sequence')
        # Initialise the string of the target sequence, so it can be iteratively mutated at each desired location
        mutated_seq = str(record.seq)
        # Initialise the header of the output sequence, so it can be iteratively updated after each mutation
        mutated_header = record.id
        # If a random seed wasn't provided, the choice will not be repeatable. Otherwise, adjust the random seed,
        # so that there is some randomness (but is still repeatable) based on the position, name, and length of
        # the record
        if random_seed:
            # Set the random seed by adding the user-provided seed, the length of the record, and the length of the
            # record header to create a target-specific repeatable choice
            adjusted_random_seed = random_seed + len(record.seq) + len(record.id)
            logging.debug(f'Using random seed: {adjusted_random_seed} ')
            random.seed(adjusted_random_seed)
        # Use random.sample to choose (without replacement) the appropriate index locations from the length of
        # the target sequence
        bases_list = sorted(random.sample(
            population=list([i for i in range(len(mutated_seq))]),
            k=minimum)
        )
        logging.debug(f'Mutating base(s) {", ".join([str(i) for i in bases_list])} in {record.id}')
        # Update the sequence and header for each mutation
        for location in bases_list:
            mutated_seq, mutated_header = mutate_bases(location=location,
                                                       record=record,
                                                       random_seed=random_seed,
                                                       mutation_dict=mutation_dict,
                                                       mutated_seq=mutated_seq,
                                                       mutated_header=mutated_header,
                                                       adjust=False)
        # Create a SeqRecord object using a Seq object of the sequence and the updated header
        mutated_record = SeqRecord(seq=Seq(mutated_seq),
                                   description=str(),
                                   id=mutated_header)
        # Return the mutated_record
        return mutated_record

    @staticmethod
    def mutate_target_range(record, minimum, maximum, mutation_dict, random_seed, output_file):
        """
        Create multiple outputs of the mutation of the target sequence by a range of the number of desired mutations
        :param record: type: Bio.SeqRecord.SeqRecord object created from the target sequence by SeqIO
        :param minimum: type int: Desired minimum number of mutations to introduce into the target sequence
        :param maximum: type int: Desired maximum number of mutations to introduce into the target sequence
        :param mutation_dict: type dict: Dictionary storing the list of potential replacement bases for each nucleotide
        :param random_seed: type int: User-provided seed value to allow consistent outputs
        :param output_file: type str: Name and path of file in which the mutated sequences are to be saved
        """
        logging.info(f'Mutating target sequence with permutations of between {minimum} and {maximum}')
        mutated_records = list()
        for i in range(minimum, maximum + 1):
            mutated_record = SequenceMutatorRange.mutate_target(
                record=record,
                minimum=i,
                mutation_dict=mutation_dict,
                random_seed=random_seed,
            )
            mutated_records.append(mutated_record)
        # Write all the mutated records to the output file
        write_outputs(output_file=output_file,
                      mutated_record=mutated_records)

    @staticmethod
    def mutate_target_range_spacing(record, minimum, maximum, spacing, mutation_dict, random_seed, output_file):
        """

        :param record: type: Bio.SeqRecord.SeqRecord object created from the target sequence by SeqIO
        :param minimum: type int: Desired minimum number of mutations to introduce into the target sequence
        :param maximum: type int: Desired maximum number of mutations to introduce into the target sequence
        :param spacing: type int: Desired spacing of mutated bases
        :param mutation_dict: type dict: Dictionary storing the list of potential replacement bases for each nucleotide
        :param random_seed: type int: User-provided seed value to allow consistent outputs
        :param output_file: type str: Name and path of file in which the mutated sequences are to be saved
        """
        logging.info(f'Mutating target sequence with permutations of between {minimum} and {maximum} with spacing of '
                     f'{spacing} bases')
        # Ensure that the spacing is not too large
        if spacing > len(record.seq):
            logging.error(f'The requested spacing, {spacing}, is greater than the length of the target sequence, '
                          f'{len(record.seq)}. Please use a value of {len(record.seq)} or less')
            raise SystemExit
        # Slice the list of all the indexes in a list of the target sequence to determine which bases are to be mutated.
        # Start the slice at the first interval of the spacing, and use the value for spacing as the step
        # e.g. for a 50 bp target, and 10 bp spacing, 10, 20, 30, 40 would be returned. 0 is skipped from starting at
        # the first interval, and 50 is skipped from the nature of slices not including the stop value
        spaced_bases = [i for i in range(len(record.seq))][spacing::spacing]
        logging.debug(f'The following base(s) will be used in the mutation analyses '
                      f'{", ".join([str(base) for base in spaced_bases])}')
        # Initialise a list to store all the mutated records
        mutated_records = list()
        # Iterate through all the minimum to maximum number of desired mutations (+1 from the range not including the
        # stop value)
        for i in range(minimum, maximum + 1):
            # Use itertools.combinations to find all the possible combinations from the list of bases to mutate, and
            # use the current number of desired mutations
            # e.g. spaced_bases = [10, 20, 30, 40], number of mutations = 2 would yield the following combinations
            # [(10, 20), (10, 30), (10, 40), (20, 30), (20, 40), (30, 40)]
            # https://stackoverflow.com/a/170248
            combinations = list(itertools.combinations(spaced_bases, i))
            # Iterate through all the tuples of locations in the list of combinations e.g. (20, 30)
            for locations in combinations:
                # Use SequenceMutatorLocation.mutate_target_locations to mutate the requested locations
                mutated_record = mutate_target_locations(
                    record=record,
                    locations=list(locations),
                    mutation_dict=mutation_dict,
                    random_seed=random_seed,
                    adjust=False,
                )
                mutated_records.append(mutated_record)
        # Write all the mutated records to the output file
        write_outputs(output_file=output_file,
                      mutated_record=mutated_records)

    def __init__(self, target_file, minimum, maximum, spacing, output_path, random_seed):
        # Expand the path of the target file as necessary
        if target_file.startswith('~'):
            self.target_file = os.path.abspath(os.path.expanduser(os.path.join(target_file)))
        else:
            self.target_file = os.path.abspath(os.path.join(target_file))
        try:
            assert os.path.isfile(self.target_file), \
                f'Cannot located target file: {self.target_file}. ' \
                f'Please ensure that you entered the correct name and path'
        except AssertionError:
            raise SystemExit
        if not output_path:
            self.output_path = os.path.join(os.path.dirname(self.target_file), 'outputs')
        else:
            if output_path.startswith('~'):
                self.output_path = os.path.abspath(os.path.expanduser(os.path.join(output_path)))
            else:
                self.output_path = os.path.abspath(os.path.join(output_path))
        make_path(self.output_path)
        self.output_file = os.path.join(self.output_path, 'mutated_target.fasta')
        if minimum <= maximum:
            self.minimum = minimum
            self.maximum = maximum
        else:
            self.minimum = maximum
            self.minimum = minimum
            logging.warning(f'The supplied minimum number of mutations: {minimum} was greater than the supplied '
                            f'maximum number of mutations: {maximum}. Using {maximum} as the minimum and {minimum} '
                            f'as the maximum.')
        self.spacing = spacing
        self.random_seed = int(random_seed)
        self.mutation_dict = {
            'A': ['C', 'G', 'T'],
            'C': ['A', 'G', 'T'],
            'G': ['A', 'C', 'T'],
            'T': ['A', 'C', 'G']
        }


class SequenceMutatorNumber(object):

    def main(self):
        # Parse the target sequence with SeqIO
        record = read_target_file(target_file=self.target_file)
        target_length = len(str(record.seq))
        if self.number > target_length:
            logging.warning(
                f'You have asked for the sequence to be mutated in {self.number} locations, however, the sequence is '
                f'only {target_length}. The number of mutations will be decreased to {target_length}')
            self.number = target_length
        self.locations = create_mutation_location_list(
            record=record,
            number=self.number,
            random_seed=self.random_seed)
        # Create the mutated SeqRecord from the target sequence
        mutated_record = mutate_target_locations(
            record=record,
            locations=self.locations,
            mutation_dict=self.mutation_dict,
            random_seed=self.random_seed,
            adjust=False,
        )
        # Write the mutated record to the output file
        write_outputs(output_file=self.output_file,
                      mutated_record=mutated_record)

    def __init__(self, target_file, number, output_path, random_seed):
        # Expand the path of the target file as necessary
        if target_file.startswith('~'):
            self.target_file = os.path.abspath(os.path.expanduser(os.path.join(target_file)))
        else:
            self.target_file = os.path.abspath(os.path.join(target_file))
        try:
            assert os.path.isfile(self.target_file), \
                f'Cannot located target file: {self.target_file}. ' \
                f'Please ensure that you entered the correct name and path'
        except AssertionError:
            raise SystemExit
        if not output_path:
            self.output_path = os.path.join(os.path.dirname(self.target_file), 'outputs')
        else:
            if output_path.startswith('~'):
                self.output_path = os.path.abspath(os.path.expanduser(os.path.join(output_path)))
            else:
                self.output_path = os.path.abspath(os.path.join(output_path))
        make_path(self.output_path)
        self.output_file = os.path.join(self.output_path, 'mutated_target.fasta')
        self.number = number
        self.locations = list()
        self.random_seed = int(random_seed)
        self.mutation_dict = {
            'A': ['C', 'G', 'T'],
            'C': ['A', 'G', 'T'],
            'G': ['A', 'C', 'T'],
            'T': ['A', 'C', 'G']
        }


class SequenceMutatorPercentIdentity(object):

    def main(self):
        # Parse the target sequence with SeqIO
        record = read_target_file(target_file=self.target_file)
        self.number = self.determine_number_of_mutations(
            record=record,
            percent_id=self.percent_id)

        self.locations = create_mutation_location_list(
            record=record,
            number=self.number,
            random_seed=self.random_seed)
        # Create the mutated SeqRecord from the target sequence
        mutated_record = mutate_target_locations(
            record=record,
            locations=self.locations,
            mutation_dict=self.mutation_dict,
            random_seed=self.random_seed,
            adjust=False
        )
        # Write the mutated record to the output file
        write_outputs(output_file=self.output_file,
                      mutated_record=mutated_record)

    @staticmethod
    def determine_number_of_mutations(record, percent_id):
        """

        :param record:
        :param percent_id:
        :return:
        """
        # Determine the length of the target sequence
        target_length = len(str(record.seq))
        # Calculate the number of bases matching the supplied percent identity. Subtract from the total target length
        # the length of the sequence multiplied by the supplied percent identity, and divide by 100.
        # Convert it to int to remove decimal places. Due to the rounding, small sequences will yield the same outputs
        # with different percent identities e.g. 50 bp sequence will return 1 base with both 98% and 99% identity
        number = target_length - int(target_length * float(percent_id) / 100)
        logging.debug(f'{number} base(s) will be mutated in {record.id}')
        return number


    def __init__(self, target_file, percent_id, output_path, random_seed):
        # Expand the path of the target file as necessary
        if target_file.startswith('~'):
            self.target_file = os.path.abspath(os.path.expanduser(os.path.join(target_file)))
        else:
            self.target_file = os.path.abspath(os.path.join(target_file))
        try:
            assert os.path.isfile(self.target_file), \
                f'Cannot located target file: {self.target_file}. ' \
                f'Please ensure that you entered the correct name and path'
        except AssertionError:
            raise SystemExit
        if not output_path:
            self.output_path = os.path.join(os.path.dirname(self.target_file), 'outputs')
        else:
            if output_path.startswith('~'):
                self.output_path = os.path.abspath(os.path.expanduser(os.path.join(output_path)))
            else:
                self.output_path = os.path.abspath(os.path.join(output_path))
        make_path(self.output_path)
        self.output_file = os.path.join(self.output_path, 'mutated_target.fasta')
        self.percent_id = percent_id
        try:
            assert 0 < self.percent_id < 100
        except AssertionError:
            logging.error('The percent identity must be between 1% and 99%')
            raise SystemExit
        self.number = int()
        self.locations = list()
        self.random_seed = int(random_seed)
        self.mutation_dict = {
            'A': ['C', 'G', 'T'],
            'C': ['A', 'G', 'T'],
            'G': ['A', 'C', 'T'],
            'T': ['A', 'C', 'G']
        }


def mutation_location(args):
    """
    Run the SequenceMutatorLocation class
    :param args: type ArgumentParser arguments
    """
    logging.info(f'Mutating target sequence in {args.target_file} at the following location(s): {args.locations}')
    mut_location = SequenceMutatorLocation(
        target_file=args.target_file,
        locations=args.locations,
        output_path=args.output_path,
        random_seed=args.random_seed
    )
    mut_location.main()


def mutation_range(args):
    """
    Run the SequenceMutatorRange class
    :param args: type ArgumentParser arguments
    """
    logging_string = f'Mutating target sequence as follows: minimum mutations {args.min}, maximum mutations {args.max}'
    if args.spacing:
        logging_string += f', with spacing of mutations: {args.spacing}'
    logging.info(logging_string)
    mut_range = SequenceMutatorRange(
        target_file=args.target_file,
        minimum=args.min,
        maximum=args.max,
        spacing=args.spacing,
        output_path=args.output_path,
        random_seed=args.random_seed
    )
    mut_range.main()


def mutation_number(args):
    """
    Run the SequenceMutatorNumber class
    :param args: type ArgumentParser arguments
    """
    logging.info(f'Inserting {args.number} mutations into target sequence in {args.target_file}')
    mut_location = SequenceMutatorNumber(
        target_file=args.target_file,
        number=args.number,
        output_path=args.output_path,
        random_seed=args.random_seed
    )
    mut_location.main()


def mutation_percent_id(args):
    """
    Run the SequenceMutatorPercentIdentity class
    :param args: type ArgumentParser arguments
    """
    logging.info(f'Creating output sequence with {args.percent_identity}% identity to target sequence in '
                 f'{args.target_file}')
    mut_location = SequenceMutatorPercentIdentity(
        target_file=args.target_file,
        percent_id=args.percent_identity,
        output_path=args.output_path,
        random_seed=args.random_seed
    )
    mut_location.main()

def cli():
    parser = ArgumentParser(description='Mutate a target sequence. Specify the number of mutations')
    subparsers = parser.add_subparsers(title='Available functionality')
    # Create a parental parser that can be inherited by subparsers
    parent_parser = ArgumentParser(add_help=False)
    parent_parser.add_argument('-t', '--target_file',
                               required=True,
                               metavar='',
                               help='Name and path of file containing a single target sequence to mutate'
                               )
    parent_parser.add_argument('-r', '--random_seed',
                               default=0,
                               type=int,
                               metavar='',
                               help='Random seed to use when the replacement base is being set. This will be adjusted '
                                    'for each sequence, so there is still some randomness'
                               )
    parent_parser.add_argument('-o', '--output_path',
                               metavar='',
                               type=str,
                               help='Name and path of folder in which the outputs are to be written. Default is:\n'
                                    'path of target file/outputs'
                               )
    parent_parser.add_argument('-d', '--debug',
                               action='store_true',
                               help='Enable debug-level messages to be printed to the console'
                               )
    # Defined location(s) subparser
    location_subparser = subparsers.add_parser(
        parents=[parent_parser],
        name='location',
        description='Specify the location(s) of the desired mutation(s)',
        formatter_class=RawTextHelpFormatter,
        help='Specify the location(s) of the desired mutation(s)'
    )
    location_subparser.add_argument(
        '-l', '--locations',
        required=True,
        type=str,
        metavar='',
        help='Mutate specific base(s) of the target sequence. Provide a comma-separated list of position(s)\n'
             'e.g. 1 or 17,22 or 5,23,73,149'
    )
    location_subparser.set_defaults(func=mutation_location)
    # Range subparser
    range_subparser = subparsers.add_parser(
        parents=[parent_parser],
        name='range',
        description='Specify the range of both the number and spacing of the desired mutation(s)',
        formatter_class=RawTextHelpFormatter,
        help='Specify the range of both the number and spacing of the desired mutation(s)'
    )
    range_subparser.add_argument(
        '-min', '--min',
        default=1,
        type=int,
        choices=[1, 2, 3, 4, 5],
        metavar='',
        help='Minimum number of mutations to insert into target sequence. Allowed range is 1 - 5'
    )
    range_subparser.add_argument(
        '-max', '--max',
        default=1,
        type=int,
        choices=[1, 2, 3, 4, 5],
        metavar='',
        help='Maximum number of mutations to insert into target sequence. Allowed range is 1 - 5. If the maximum is '
             'greater than the minimum, multiple outputs will be created'
    )
    range_subparser.add_argument(
        '-s', '--spacing',
        default=0,
        type=int,
        metavar='',
        help='If provided, the script will create permutations of the target sequence with mutations at the desired '
             'spacing. Minimum is 5'
    )
    range_subparser.set_defaults(func=mutation_range)
    # Number subparser
    number_subparser = subparsers.add_parser(
        parents=[parent_parser],
        name='number',
        description='Specify the number of desired mutation(s)',
        formatter_class=RawTextHelpFormatter,
        help='Specify the number of desired mutation(s)'
    )
    number_subparser.add_argument(
        '-n', '--number',
        default=0,
        type=int,
        metavar='',
        help='The number of mutations to introduce into the target sequence. If the target has fewer nucleotides than '
             'the number of mutations, the entire sequence will be mutated.'
    )
    number_subparser.set_defaults(func=mutation_number)
    # Percent ID subparser
    percent_id_subparser = subparsers.add_parser(
        parents=[parent_parser],
        name='percent_id',
        description='Specify the percent identity of the sequence following the insertion of mutations to the '
                    'target sequence',
        formatter_class=RawTextHelpFormatter,
        help='Specify the percent identity of the sequence following the insertion of mutations to the '
                    'target sequence'
    )
    percent_id_subparser.add_argument(
        '-p', '--percent_identity',
        default=99,
        type=int,
        metavar='',
        help='The percent identity of the output sequence to the original sequence following the insertion of mutations'
             ' e.g. if the sequence is 100 bp, and you wish to have a 90% identity, 10 mutations will be inserted. '
             'Default is 99.'
    )
    percent_id_subparser.set_defaults(func=mutation_percent_id)
    # Get the arguments into an object
    arguments = parser.parse_args()
    # Set up the logging
    SetupLogging(debug=arguments.debug)
    # Run the appropriate function for each sub-parser.
    if hasattr(arguments, 'func'):
        # Run the appropriate function
        arguments.func(arguments)
    # If the 'func' attribute doesn't exist, display the basic help for the parental parser
    else:
        parser.parse_args(['-h'])
    logging.info('Mutation complete')
    # Prevent the arguments being printed to the console (they are returned in order for the tests to work)
    sys.stderr = open(os.devnull, 'w')
    return arguments


if __name__ == '__main__':
    cli()
