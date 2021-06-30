#!/usr/bin/env python

import os
import sys
import re
import uproot
import numpy as np

from matplotlib import pyplot as plt
from tabulate import tabulate
from pprint import pprint

pjoin = os.path.join

def main():
    inpath = 'output/merged_2021-06-28_vbfhinv_ULv8_05Feb21_lepton_uncs/lepton_sf_uncs.root'
    f = uproot.open(inpath)

    variations = [k.decode('utf-8').replace(';1','') for k in f.keys()]
    avg_variations = []
    for v in variations:
        unc = np.abs(f[v].values.mean() - 1) * 100
        avg_variations.append((v, unc))

    outtag = inpath.split('/')[-2]
    outfile = pjoin(f'./output/{outtag}', 'lepton_sf_uncs_table.txt')
    with open(outfile, 'w+') as f:
        f.write('Lepton SF uncertainties')
        f.write('\n')
        f.write(tabulate(avg_variations, headers=['Uncertainty', 'Average Unc (%)']))

    print(f'Fiie saved: {outfile}')

if __name__ == '__main__':
    main()