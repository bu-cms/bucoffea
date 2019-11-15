#!/usr/bin/env python
import os
import pickle
import socket
import subprocess
import htcondor

pjoin = os.path.join
def condor_submit(jobfile):
    host = socket.gethostname()
    if 'lxplus' in host:
        cmd = ["condor_submit", jobfile]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)

        stdout, stderr = proc.communicate()
        if proc.returncode != 0:
            raise RuntimeError(f"Condor submission failed. Stderr:\n {stderr}.")
        jobid = stdout.split()[-1].decode('utf-8').replace('.','')
    elif 'fnal' in host:
        cmd = ["bash","/usr/local/bin/condor_submit", jobfile]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)

        stdout, stderr = proc.communicate()
        if proc.returncode != 0:
            raise RuntimeError(f"Condor submission failed. Stderr:\n {stderr}.")
        jobid = stdout.split()[-1].decode('utf-8').replace('.','')
    return jobid



def read_logs(directories):
    logs = []
    for directory in directories:
        for path, _, files in os.walk(directory):
            logs.extend([os.path.join(path, x) for x in files if x.startswith("log_")])
    return logs

class ConJob():
    '''Condor Job'''
    def __init__(self,log):
        self.log = log
        self.status = None
        self.name = os.path.basename(self.log).replace("log_","").replace(".txt","")
        self.resubcount = 0
        self.update()

    @property
    def log(self):
        return self._log

    @log.setter
    def log(self,log):
        log = os.path.abspath(log)
        if not os.path.exists(log):
            raise FileNotFoundError(f'Cannot find log file: {log}')
        self._log = log

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self,name):
        self._name = name

    @property
    def status(self):
        return self._status

    @status.setter
    def status(self,status):
        self._status = status

    @property
    def code(self):
        return self._code

    @code.setter
    def code(self,code):
        self._code = code

    @property
    def cluster(self):
        return self._cluster

    @cluster.setter
    def cluster(self,cluster):
        self._cluster = cluster

    @property
    def resubcount(self):
        return self._resubcount

    @resubcount.setter
    def resubcount(self,resubcount):
        self._resubcount = resubcount

    @property
    def runtime(self):
        return self._runtime

    @runtime.setter
    def runtime(self,runtime):
        self._runtime = runtime

    def update(self):
        # What's done is done
        if self.status == 'JOB_TERMINATED':
            return

        # Updating works by opening the log file,
        # Looping through and only keeping the last event,
        # which really tells us what's going on.
        # This is not very efficient, so an alternative
        # implementation would be welcome
        jel = htcondor.JobEventLog(self._log)

        first = None
        try:
            for event in jel.events(stop_after=0):
                if not first:
                    first = event
                latest = event
            try:
                self._code = latest["ReturnValue"]
            except KeyError:
                self._code = "-"
            self.status = str(htcondor.JobEventType.values[latest.type])
            self.cluster = latest.cluster
            self.runtime = latest.timestamp - first.timestamp
        except OSError:
            self.code = "-"
            self.status = "NOPARSE"
            self.cluster = "-"
            self.runtime = -1
        finally:
            jel.close()

    def jdl(self):
        return self.log.replace('log_','job_').replace('.txt','.jdl')

    def resubmit(self):
        self.status = None
        condor_submit(self.jdl())

        self.resubcount =  self.resubcount + 1

class ConMan():
    '''
    Condor manager.

    Wraps a list of condor jobs.
    '''
    def __init__(self,directory):
        self.directory = directory
        self.pkl = pjoin(directory, "jobs.pkl")
        self.autoresub = False

        if os.path.exists(self.pkl):
            self.jobs = pickle.load( open( self.pkl, "rb" ) )
        else:
            self.init_jobs()

    @property
    def directory(self):
        return self._directory

    @directory.setter
    def directory(self,directory):
        directory = os.path.abspath(directory)
        if not os.path.exists(directory):
            raise FileNotFoundError(f'Cannot find directory: {directory}')
        self._directory = directory

    @property
    def jobs(self):
        return self._jobs

    @jobs.setter
    def jobs(self,jobs):
        self._jobs = sorted(jobs, key=lambda x: x.name)

    @property
    def autoresub(self):
        return self._autoresub

    @autoresub.setter
    def autoresub(self,autoresub):
        self._autoresub = autoresub

    def init_jobs(self):
        logs = read_logs([self.directory])
        jobs = list(map(ConJob, logs))
        self.jobs = jobs

    def update(self):
        for j in self.jobs:
            j.update()
        if self.autoresub:
            self.resubmit_failed(max_resub=3)

    def resubmit_failed(self, max_resub=None):
        count = 0
        for j in self.jobs:
            if (not j.code  in [0, '-']) or (j.status == 'JOB_ABORTED'):
                if max_resub and max_resub <= j.resubcount:
                    continue
                j.resubmit()
                count += 1
        return count

    def save(self):
        pickle.dump( self.jobs, open( self.pkl, "wb" ) )
