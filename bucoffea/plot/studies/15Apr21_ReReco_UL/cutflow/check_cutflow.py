#!/usr/bin/env python
import os
import re
import sys
import csv
import tabulate
import numpy as np
import mplhep as hep

from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
from coffea import hist
from bucoffea.helpers.paths import bucoffea_path
from klepto.archives import dir_archive
from pprint import pprint

pjoin = os.path.join


def dump_cutflows(acc_dict, dataset_regex):
    cutflows = {}
    for key, acc in acc_dict.items():
        cutflows[key] = acc['cutflow_sr_vbf']

    outdir = './output'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # Get the dataset list that we're interested in
    datasetlist_eoy = [d for d in list(cutflows['EOY'].keys()) if re.match(dataset_regex, d)]
    datasetlist_ul = [d for d in list(cutflows['UL'].keys()) if re.match(dataset_regex, d)]

    for dataset_eoy in datasetlist_eoy:
        c_EOY = cutflows['EOY'][dataset_eoy]
        outpath = pjoin(outdir, f'eoy_{dataset_eoy}.csv')
        with open(outpath, 'w+') as f:
            writer = csv.writer(f)
            writer.writerow(['Cut', 'Count'])
            for cut, count in c_EOY.items():
                writer.writerow([cut, count])

        print(f'File saved: {outpath}')

    for dataset_ul in datasetlist_ul:
        c_UL = cutflows['UL'][dataset_ul]
        outpath = pjoin(outdir, f'ul_{dataset_ul}.csv')
        with open(outpath, 'w+') as f:
            writer = csv.writer(f)
            writer.writerow(['Cut', 'Count'])
            for cut, count in c_UL.items():
                writer.writerow([cut, count])

        print(f'File saved: {outpath}')

def main():
    acc_dict = {
        'EOY' : dir_archive( bucoffea_path('submission/merged_2021-03-19_vbfhinv_03Sep20v7_leadak4_notinHF')),
        'UL' : dir_archive( bucoffea_path('submission/merged_2021-05-03_vbfhinv_ULv8_05Feb21_trigSF_update')),
    }

    distribution = 'cutflow_sr_vbf'
    for acc in acc_dict.values():
        acc.load(distribution)

    dump_cutflows(acc_dict, dataset_regex='ZJetsToNuNu.*HT.*')

if __name__ == '__main__':
    main()