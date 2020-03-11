#!/usr/bin/env python

import os
import sys
import re
from bucoffea.plot.util import merge_datasets, merge_extensions, scale_xs_lumi
from coffea import hist
from matplotlib import pyplot as plt
from klepto.archives import dir_archive
from pprint import pprint

# Plot aesthetics
colors = plt.rcParams['axes.prop_cycle'].by_key()['color'][:5]

recoil_bins_2016 = [250,  280,  310,  340,  370,  400,  430,  470,  510, 550,  590,  640,  690,  740,  790,  840,  900]

REBIN = {
    'mjj' : hist.Bin('mjj', r'$M_{jj}$ (GeV)', list(range(200,800,300)) + list(range(800,2000,400)) + [2000, 2750, 3500]),
    'recoil' : hist.Bin('recoil','Recoil (GeV)', recoil_bins_2016),
    'ak4_pt0' : hist.Bin('jetpt',r'Leading AK4 jet $p_{T}$ (GeV)',list(range(80,600,20)) + list(range(600,1000,20)) ),
    'ak4_pt1' : hist.Bin('jetpt',r'Trailing AK4 jet $p_{T}$ (GeV)',list(range(40,600,20)) + list(range(600,1000,20)) ),
    'ak4_pt' : hist.Bin('jetpt',r'All AK4 jet $p_{T}$ (GeV)',list(range(100,600,20)) + list(range(600,1000,20)) )
}

def plot_dist(acc, param, regex, region, outtag):
    '''
    Given the input accumulator, plot distribution of param for ALL variations for given dataset and region.
    =================
    PARAMETERS:
    =================
    acc:    Input accumulator containing all histograms.
    param:  The parameter to be plotted.
    regex:  Regular expression matching the dataset name.
    region: The region from which data will be taken.
    outtag: Tag for naming output directory.
    '''
    acc.load(param)
    h = acc[param]

    # Rebin, if neccessary
    if param in REBIN.keys():
        new_bin = REBIN[param]
        h = h.rebin(h.axis(param), new_bin)
    
    # Merge extensions/datasets, scale w.r.t. xs and lumi
    h = merge_extensions(h, acc, reweight_pu=False)
    scale_xs_lumi(h)
    h = merge_datasets(h)

    # Pick the relevant region and dataset
    h = h[re.compile(regex)].integrate('dataset')
    h = h.integrate('region', re.compile(f'{region}.*'))

    # Plot for each variation
    var_to_label = {
        '' : 'Nominal',
        '_jerup' : 'JER up',
        '_jerdown' : 'JER down',
        '_jesup' : 'JES up',
        '_jesdown' : 'JES down',
    }

    fig, ax = plt.subplots(1,1)
    # Rename variations so that they show up in legend 
    # when plotted by plot1d func
    for var in h.identifiers('var'):
        var.label = var_to_label[var.name]

    hist.plot1d(h, ax=ax, overlay='var')
    ax.set_ylabel('Counts')        

    fig.savefig('test.pdf')

def main():
    inpath = sys.argv[1]
    
    acc = dir_archive(
        inpath,
        memsize=1e3,
        compression=0,
        serialized=True
    )

    acc.load('sumw')
    acc.load('sumw2')

    plot_dist(acc, param='recoil', regex='WJetsToLNu.*2017', region='sr_vbf', outtag=None)

if __name__ == '__main__':
    main()

