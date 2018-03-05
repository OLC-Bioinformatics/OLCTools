#!/usr/bin/env python
from threading import Thread, Lock
from queue import Queue
from accessoryFunctions.accessoryFunctions import printtime, GenObject, make_path, run_subprocess, write_to_logfile
import os
__author__ = 'adamkoziol'

threadlock = Lock()


class Prodigal(object):

    def predictthreads(self):
        printtime('Performing gene predictions', self.start)
        # Create the threads for the analyses
        for sample in self.metadata:
            if sample.general.bestassemblyfile != 'NA':
                threads = Thread(target=self.predict, args=())
                threads.setDaemon(True)
                threads.start()
        for sample in self.metadata:
            # Create the .prodigal attribute
            sample.prodigal = GenObject()
            if sample.general.bestassemblyfile != 'NA':
                self.predictqueue.put(sample)
        self.predictqueue.join()

    def predict(self):
        while True:
            sample = self.predictqueue.get()
            # Populate attributes
            sample.prodigal.reportdir = os.path.join(sample.general.outputdirectory, 'prodigal')
            sample.prodigal.results_file = os.path.join(sample.prodigal.reportdir,
                                                        '{}_prodigalresults.sco'.format(sample.name))
            sample.prodigal.results = sample.prodigal.results_file
            sample.commands.prodigal = 'prodigal -i {in1} -o {out1} -f sco -d {genes}'\
                .format(in1=sample.general.bestassemblyfile,
                        out1=sample.prodigal.results_file,
                        genes=os.path.join(sample.prodigal.reportdir, '{}_genes.fa'.format(sample.name)))
            # Create the folder to store the reports
            make_path(sample.prodigal.reportdir)
            # Determine if the report already exists, and that it is not empty
            size = 0
            if os.path.isfile(sample.prodigal.results_file):
                size = os.stat(sample.prodigal.results_file).st_size
            if not os.path.isfile(sample.prodigal.results_file) or size == 0:
                # Run the command
                out, err = run_subprocess(sample.commands.prodigal)
                threadlock.acquire()
                write_to_logfile(sample.commands.prodigal, sample.commands.prodigal, self.logfile,
                                 sample.general.logout, sample.general.logerr, None,
                                 None)
                write_to_logfile(out, err, self.logfile, sample.general.logout, sample.general.logerr, None, None)
                threadlock.release()
            self.predictqueue.task_done()

    def prodigalparse(self):
        printtime('Parsing gene predictions', self.start)
        for sample in self.metadata:
            sample.prodigal.predictedgenestotal = 0
            sample.prodigal.predictedgenesover3000bp = 0
            sample.prodigal.predictedgenesover1000bp = 0
            sample.prodigal.predictedgenesover500bp = 0
            sample.prodigal.predictedgenesunder500bp = 0
            if sample.general.bestassemblyfile != 'NA':
                with open(sample.prodigal.results, 'r') as results:
                    for line in results:
                        if line.startswith('>'):
                            start = int(line.split('_')[1])
                            end = int(line.split('_')[2])
                            length = abs(start - end)
                            sample.prodigal.predictedgenestotal += 1
                            if length > 3000:
                                sample.prodigal.predictedgenesover3000bp += 1
                            elif length > 1000:
                                sample.prodigal.predictedgenesover1000bp += 1
                            elif length > 500:
                                sample.prodigal.predictedgenesover500bp += 1
                            else:
                                sample.prodigal.predictedgenesunder500bp += 1

    def __init__(self, inputobject):
        self.metadata = inputobject.runmetadata.samples
        self.start = inputobject.starttime
        self.logfile = inputobject.logfile
        self.predictqueue = Queue()
        self.parsequeue = Queue()
        self.predictthreads()
        self.prodigalparse()
