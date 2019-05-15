#!/usr/bin/env python3
from accessoryFunctions.accessoryFunctions import GenObject, run_subprocess, write_to_logfile
import logging
import pandas
import shutil
import os
__author__ = 'adamkoziol'


class ECTyper(object):

    def main(self):
        self.run_ec_typer()
        self.move_report()
        self.populate_dictionary()
        self.populate_metadata()

    def run_ec_typer(self):
        """
        Run the ECTyper system call
        """
        logging.info('Running ECTyper')
        ec_cmd = 'ectyper  -i {input_folder} -o {output_folder}'.format(input_folder=self.assembly_path,
                                                                        output_folder=self.report_path)
        if not os.path.isfile(self.report_output) and not os.path.isfile(self.report_final):
            out, err = run_subprocess(ec_cmd)
            write_to_logfile(out=out,
                             err=err,
                             logfile=self.logfile)

    def move_report(self):
        """
        Rename the report
        """
        if os.path.isfile(self.report_output):
            shutil.move(src=self.report_output,
                        dst=self.report_final)

    def populate_dictionary(self):
        """
        Parse the report, and populate a dictionary with the extracted values
        """
        with open(self.report_final, 'r') as report:
            next(report)
            for line in report:
                sample_name = line.split('\t')[0]
                if line.split('\t')[1].rstrip() == 'No serotyping-specific genes found':
                    self.nesteddictionary[sample_name] = {
                        'o_type': 'ND',
                        'h_type': 'ND'
                    }
                else:
                    self.nesteddictionary[sample_name] = {
                        'o_type': line.split('\t')[1].rstrip(),
                        'h_type': line.split('\t')[2].rstrip()
                    }

    def populate_metadata(self):
        """
        Populate the metadata object using the dictionary with values extracted from the report
        """
        for sample in self.metadata:
            # Initialise a GenObject to store the ECTyper results
            sample.ectyper = GenObject()
            # Populate the GenObject with the header: value pairs from the report
            try:
                for header, value in self.nesteddictionary[sample.name].items():
                    sample.ectyper[header] = value
            except KeyError:
                sample.ectyper.o_type = 'ND'
                sample.ectyper.h_type = 'ND'

    def __init__(self, metadata, report_path, assembly_path, threads, logfile):
        self.metadata = metadata.samples
        self.report_path = report_path
        self.report_output = os.path.join(self.report_path, 'output.tsv')
        self.report_final = os.path.join(self.report_path, 'ec_report.tsv')
        self.assembly_path = assembly_path
        self.logfile = logfile
        self.threads = threads
        self.nesteddictionary = dict()
