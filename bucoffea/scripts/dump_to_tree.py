#!/usr/bin/env python

import os
import sys
from coffea.util import load
import uproot
import argparse
import numpy as np

pjoin = os.path.join

def parse_cli():
    parser = argparse.ArgumentParser()
    parser.add_argument('inpath', help='The path to directory containing input coffea files.')
    parser.add_argument('--year', help='The year to look at.')
    args = parser.parse_args()
    return args

args = parse_cli()
inpath = args.inpath
year = args.year

region_labels = ['2m', '2e']

for region_label in region_labels:
    tree_name = f'tree_{region_label}'
    
    # Find region and branch names
    variables = []
    regions = []
    for fname in os.listdir(inpath):
        # Only look at the requested year
        if year not in fname:
            continue
        path = pjoin(inpath, fname)
        acc = load(path)
        
        for region in acc[tree_name].keys():
            regions.append(region)
            variables.extend(acc[tree_name][region].keys())
    
    # Combine
    outdir = f'./output/trees_for_uw/{os.path.basename(inpath)}'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    outfile = pjoin(outdir, f'tree_{region_label}_{year}.root')
    with uproot.recreate(outfile) as f:
        for region in set(regions):
            f[region] = uproot.newtree({x:np.float64 for x in variables})
            for fname in os.listdir(inpath):
                # Only look at the requested year
                if year not in fname:
                    continue
                path = pjoin(inpath, fname)
                print(f'Processing file: {path}')
                acc = load(path)
                d = {x: acc[tree_name][region][x].value for x in variables}
                f[region].extend(d)
            
