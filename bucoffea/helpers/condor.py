#!/usr/bin/env python

import subprocess
import socket
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


