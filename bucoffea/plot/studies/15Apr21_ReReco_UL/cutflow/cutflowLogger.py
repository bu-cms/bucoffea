#!/usr/bin/env python

import os
import sys
import re
import argparse
import numpy as np

from klepto.archives import dir_archive
from pprint import pprint
from coffea import processor
from tabulate import tabulate

pjoin = os.path.join

def parse_cli():
    parser = argparse.ArgumentParser()
    parser.add_argument('inpath', help='Path to merged coffea files.')
    parser.add_argument('--region', help='The region to look at.', default='cr_1e_vbf')
    parser.add_argument('--dregex', help='Dataset regex to look at.', default='Single(Electron|Photon)_.*_2017')
    args = parser.parse_args()
    return args

class CutflowLogger():
    def __init__(self, region, inpath, dataset_regex):
        self.inpath = inpath
        self.region = region
        self.dataset_regex = dataset_regex
        self._setup_acc()

    def _setup_acc(self):
        self.accumulator = dir_archive(self.inpath)
        self.accumulator.load(f'cutflow_{self.region}')
    
    def dump_cutflow(self):
        init_cf = processor.defaultdict_accumulator(int)
        
        cf = self.accumulator[f'cutflow_{self.region}']
        all_datasets = cf.keys()
        
        for dataset in all_datasets:
            if not re.match(self.dataset_regex, dataset):
                continue
            
            init_cf.add(cf[dataset])

        print(tabulate(
            init_cf.items(),
            headers=['Cut','# of Events'],
            disable_numparse=True
            ))

def main():
    args = parse_cli()
    cl = CutflowLogger(region=args.region, inpath=args.inpath, dataset_regex=args.dregex)
    cl.dump_cutflow()

if __name__ == '__main__':
    main()