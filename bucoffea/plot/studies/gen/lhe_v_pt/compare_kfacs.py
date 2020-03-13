#!/usr/bin/env python

import os
import sys
import re
import pickle
import numpy as np
from matplotlib import pyplot as plt
from pprint import pprint

pjoin = os.path.join

def compare_kfactors(cases):
    '''
    For two cases, corresponding to two different LO sample usage for k-factor derivation,
    construct a comparison plot between the two k-factors.
    '''
    assert len(cases) == 2

    pickledir = './pkl_kfacs'

    picklefiles = [pjoin(pickledir, f'{case}_kfac.pkl') for case in cases]

    kfacs = {}
    uncs = {}

    # From the pickle files, get k-factors and 
    # associated stat uncertainties for each case
    # These are stored in two dictionaries
    for case, picklefile in zip(cases, picklefiles):
        with open(picklefile, 'rb') as pf:
            kfacs[case], uncs[case] = pickle.load(pf)

    # Calculate differences
    old, new = cases
    diff = kfacs[new] - kfacs[old]
    diff /= uncs[old]
    
    # Flatten to 1D array
    diff = diff.flatten()

    # Remove possible inf values (value will be outside binning range)
    diff[np.isinf(diff)] = -10

    # Plot the diff quantity as a histogram
    fig, ax = plt.subplots(1,1)
    ax.hist(diff, bins=np.arange(-6,6,0.5))

    pretty_title =  {
        'gjets_dr_16' : 'GJets_DR-0p4_HT_2016',
        'gjets_dr_17' : 'GJets_DR-0p4_HT_2017',
        'gjets_ht_16' : 'GJets_HT_2016',
        'gjets_ht_17' : 'GJets_HT_2017'
    }

    title = ' VS '.join([pretty_title[case] for case in cases]) 
    ax.set_title(title)

    ax.set_xlabel('k-fac diff / stat unc')
    ax.set_ylabel('Counts')

    outdir = './output/kfactor_comparison'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    outpath = pjoin(outdir, f'{old}_vs_{new}.pdf')
    fig.savefig(outpath)

def main():
    cases = ['gjets_dr_17', 'gjets_ht_17']
    compare_kfactors(cases)

if __name__ == '__main__':
    main()