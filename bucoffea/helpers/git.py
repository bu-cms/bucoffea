#!/usr/bin/env python

import subprocess

def git_rev_parse():
    return subprocess.check_output(['git', 'rev-parse', 'HEAD']).decode('utf-8')
    
def git_diff():
    return subprocess.check_output(['git', 'diff']).decode('utf-8')
