#!/usr/bin/env python

import os
import sys
import re
import csv
import uproot
import numpy as np
from glob import glob
from tqdm import tqdm
from pprint import pprint

pjoin = os.path.join

def dump_run_event_lumi(infiles, outfile):
    '''Dump run,lumi,event information for events passing VBF SR selection & failing HF shape cuts.'''
    with open(outfile, 'w+') as f:
        writer = csv.writer(f)
        writer.writerow(['Run','Event','Lumi'])
        for infile in tqdm(infiles):
            tree = uproot.open(infile)['sr_vbf_fail_hfcuts']
            events = tree['event'].array()
            runs = tree['run'].array()
            lumis = tree['lumi'].array()

            for event, run, lumi in zip(events, runs, lumis):
                writer.writerow([run, event, lumi])


def main():
    indir = sys.argv[1]
    infiles = glob(pjoin(indir, 'tree*.root'))

    # Set output directory based on the input directory
    outdir = indir.replace('input', 'output')
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outfile = pjoin(outdir, 'run_event_lumi.csv')

    dump_run_event_lumi(infiles, outfile)

if __name__ == '__main__':
    main()