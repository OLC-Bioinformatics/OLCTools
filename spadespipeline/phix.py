#!/usr/bin/env python 3
from accessoryFunctions.accessoryFunctions import printtime
from interop import py_interop_run_metrics, py_interop_run, py_interop_summary, py_interop_comm
import os

__author__ = 'adamkoziol'


class PhiX(object):

    def main(self):
        printtime('Attempting to extract PhiX mapping data', self.start)
        interop_folder = os.path.join(self.path, 'InterOp')
        # Determine if the InterOp folder is present
        if os.path.isdir(interop_folder):
            # Try to extract the relevant information
            try:
                self.interop_parse()
            # interop.py_interop_comm.file_not_found_exception: RunParameters.xml required for legacy run folders with
            # missing channel names
            except py_interop_comm.file_not_found_exception:
                for sample in self.metadata:
                    sample.run.error_rate = 'NA'
                    sample.run.phix_aligned = 'NA'
        else:
            # Create attributes reflecting the lack of the InterOp folder
            for sample in self.metadata:
                sample.run.error_rate = 'NA'
                sample.run.phix_aligned = 'NA'

    def interop_parse(self):
        """
        Use interop to parse the files in the InterOp folder to extract the number of reads mapping to PhiX as well as
        the error rate
        """
        # Parse the files and load the data
        run_metrics = py_interop_run_metrics.run_metrics()
        valid_to_load = py_interop_run.uchar_vector(py_interop_run.MetricCount, 0)
        py_interop_run_metrics.list_summary_metrics_to_load(valid_to_load)
        run_metrics.read(self.path, valid_to_load)
        summary = py_interop_summary.run_summary()
        py_interop_summary.summarize_run_metrics(run_metrics, summary)
        # PhiX error rate for run over all "usable cycles"
        errorrate = summary.total_summary().error_rate()
        # Percent aligned PhiX
        pctaligned = summary.total_summary().percent_aligned()
        # Add the error rate and the percent of reads that align to PhiX to the metadata object
        for sample in self.metadata:
            sample.run.error_rate = '{:.2f}'.format(errorrate)
            sample.run.phix_aligned = '{:.2f}'.format(pctaligned)

    def __init__(self, inputobject):
        self.path = inputobject.path
        self.start = inputobject.starttime
        self.metadata = inputobject.runmetadata.samples
