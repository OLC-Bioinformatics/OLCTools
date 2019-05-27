#!/usr/bin/env python 3
from accessoryFunctions.accessoryFunctions import run_subprocess
import logging
import os

__author__ = 'adamkoziol'


class PhiX(object):

    def main(self):
        logging.info('Attempting to extract PhiX mapping data')
        interop_folder = os.path.join(self.path, 'InterOp')
        # Determine if the InterOp folder is present
        if os.path.isdir(interop_folder):
            # Try to extract the relevant information
            try:
                self.run_interop_summary()
                self.interop_parse()
            # interop.py_interop_comm.file_not_found_exception: RunParameters.xml required for legacy run folders with
            # missing channel names
            except:
                for sample in self.metadata:
                    sample.run.actual_yield = 'ND'
                    sample.run.projected_yield = 'ND'
                    sample.run.phix_aligned = 'ND'
                    sample.run.error_rate = 'ND'
                    sample.run.over_q30 = 'ND'
        else:
            # Create attributes reflecting the lack of the InterOp folder
            for sample in self.metadata:
                sample.run.actual_yield = 'ND'
                sample.run.projected_yield = 'ND'
                sample.run.phix_aligned = 'ND'
                sample.run.error_rate = 'ND'
                sample.run.over_q30 = 'ND'

    def run_interop_summary(self):
        """
        Run the interop_summary script to create a .csv file summarising the run stats
        """
        interop_summary_cmd = 'interop_summary {path} > {report}'.format(path=self.path,
                                                                         report=self.interop_summary_report)
        if not os.path.isfile(self.interop_summary_report):
            run_subprocess(command=interop_summary_cmd)

    def interop_parse(self):
        """
        Parse the .csv summary file to extract the percent PhiX aligned and the error rate
        """
        if os.path.isfile(self.interop_summary_report):
            with open(self.interop_summary_report, 'r') as summary:
                for line in summary:
                    # Only the line starting with 'Total' is required
                    # Level        Yield    Projected Yield	  Aligned  	Error Rate  Intensity C1   	%>=Q30
                    # Total        17.23	17.23	          1.87	    2.75	    246	            65.19
                    if line.startswith('Total'):
                        total, run_yield, projected_yield, phix_aligned, error_rate, intensity_c1, over_q30 \
                            = line.replace(' ', '').rstrip().split(',')
                        for sample in self.metadata:
                            sample.run.actual_yield = run_yield
                            sample.run.projected_yield = projected_yield
                            sample.run.phix_aligned = phix_aligned
                            sample.run.error_rate = error_rate
                            sample.run.over_q30 = over_q30
        else:
            for sample in self.metadata:
                sample.run.actual_yield = 'ND'
                sample.run.projected_yield = 'ND'
                sample.run.phix_aligned = 'ND'
                sample.run.error_rate = 'ND'
                sample.run.over_q30 = 'ND'

    def __init__(self, inputobject):
        self.path = inputobject.path
        self.metadata = inputobject.runmetadata.samples
        self.reportpath = inputobject.reportpath
        self.interop_summary_report = os.path.join(self.reportpath, 'interop_summary_report.csv')
