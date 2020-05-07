#!/usr/bin/env python
import re
import sys
from coffea.util import load
import uproot
import numpy as np
from pprint import pprint
infiles = sys.argv[1:]

# Find region and branch names
datatypes = {}
tree_by_variable = {}
variables = []
regions = []
for fname in infiles:
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
with uproot.recreate("tree.root",compression=uproot.ZLIB(4)) as f:
    for region in set(regions):
        for fname in infiles:
            acc = load(fname)
            d = {x: acc[tree_by_variable[x]][region][x].value for x in variables}

            # Remove empty entries
            to_remove = []
            for k, v in d.items():
                if not len(v):
                    to_remove.append(k)
            for k in to_remove:
                d.pop(k)

            print(f.keys())
            # d.pop('two_electrons')
            if not (region in [re.sub(";.*","",x.decode("utf-8")) for x in f.keys()]):
                f[region] = uproot.newtree({x : np.float64 for x in d.keys()})

            lengths = set()
            for k,v in d.items():
                print(k, len(v))
                lengths.add(len(v))
            print(lengths)
            # write
            f[region].extend(d)
