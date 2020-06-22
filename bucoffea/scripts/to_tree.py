#!/usr/bin/env python
import argparse
import os
import re
import sys
from collections import defaultdict
from pprint import pprint

import numpy as np
import uproot
from coffea.util import load

pjoin = os.path.join

def files_by_dataset(filelist):
    # Split input files by dataset by dataset
    datasets = defaultdict(list)
    for ifile in filelist:
        m = re.match("(vbfhinv|monojet)_(.*_201(6|7|8)([A-Z])?)(_\d*)?.coffea", os.path.basename(ifile))

        if not m:
            print(f"Skipping file {ifile}")
            continue

        dataset = m.groups()[1]

        datasets[dataset].append(ifile)

    return datasets

def make_trees(args):



    filelists = files_by_dataset(args.files)
    # The output for each dataset will be written into a separate file
    for dataset, files in filelists.items():
        # Find region and branch names
        datatypes = {}
        tree_by_variable = {}
        variables = []
        regions = []

        # Scout out what branches there are
        for fname in files:
            acc = load(fname)
            
            treenames = [x for x in map(str,acc.keys()) if x.startswith("tree")]

            for tn in treenames:
                datatype = tn.split("_")[-1]
                for region in acc[tn].keys():
                    vars = acc[tn][region].keys()
                    regions.append(region)
                    variables.extend(vars)
                    for v in vars:
                        datatypes[v] = np.float64 #getattr(np, datatype)
                        tree_by_variable[v] = tn

        # Combine
        with uproot.recreate(pjoin(args.outdir, f"tree_{dataset}.root"),compression=uproot.ZLIB(4)) as f:
            for region in set(regions):
                for fname in files:
                    acc = load(fname)
                    d = {x: acc[tree_by_variable[x]][region][x].value for x in variables}

                    # Remove empty entries
                    to_remove = []
                    for k, v in d.items():
                        if not len(v):
                            to_remove.append(k)
                    for k in to_remove:
                        d.pop(k)

                    if not len(d):
                        continue
                    if not (region in [re.sub(";.*","",x.decode("utf-8")) for x in f.keys()]):
                        f[region] = uproot.newtree({x : np.float64 for x in d.keys()})

                    lengths = set()
                    for k,v in d.items():
                        lengths.add(len(v))
                    assert(len(lengths) == 1)
                    # write
                    f[region].extend(d)

def commandline():
    parser = argparse.ArgumentParser(prog='Convert coffea files to TTrees.')
    parser.add_argument('files', type=str, nargs='+', help='Input folder to use.')
    parser.add_argument('--outdir', type=str, default='./trees', help='Output directory')


    args = parser.parse_args()
    args.outdir = os.path.abspath(args.outdir)
    return args


def main():
    args = commandline()
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    make_trees(args)

if __name__ == "__main__":
    main()