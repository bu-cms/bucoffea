#!/usr/bin/env python

import subprocess

def condor_submit(jobfile):
    cmd = ["condor_submit", jobfile]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)

    stdout, stderr = proc.communicate()
    if proc.returncode != 0:
        raise RuntimeError(f"Condor submission failed. Stderr:\n {stderr}.")

    return stdout


