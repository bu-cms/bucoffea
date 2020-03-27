#!/usr/bin/env python

## Compare nominal sumw from 06Jan20 and 19Feb20 inputs

import os
import sys
import re
import argparse
import numpy as np
from bucoffea.plot.util import merge_datasets, scale_xs_lumi, merge_extensions
from coffea import hist
from matplotlib import pyplot as plt
from klepto.archives import dir_archive
from pprint import pprint

pjoin = os.path.join

recoil_bins_2016 = [250,  280,  310,  340,  370,  400,  430,  470,  510, 550,  590,  640,  690,  740,  790,  840,  900]

REBIN = {
    'mjj' : hist.Bin('mjj', r'$M_{jj}$ (GeV)', list(range(200,800,300)) + list(range(800,2000,400)) + [2000, 2750, 3500]),
    'recoil' : hist.Bin('recoil','Recoil (GeV)', recoil_bins_2016),
    'ak4_pt0' : hist.Bin('jetpt',r'Leading AK4 jet $p_{T}$ (GeV)',list(range(80,700,20)) ),
    'ak4_pt1' : hist.Bin('jetpt',r'Trailing AK4 jet $p_{T}$ (GeV)',list(range(40,400,20)) ),
    'ak4_pt' : hist.Bin('jetpt',r'All AK4 jet $p_{T}$ (GeV)',list(range(100,700,20)) )
}

def parse_commandline():
    parser = argparse.ArgumentParser()
    parser.add_argument('--path_to_06Jan', help='Path to coffea files for 06Jan20 inputs.')
    parser.add_argument('--path_to_19Feb', help='Path to coffea files for 19Feb20 inputs.')
    parser.add_argument('--distribution', help='Distribution to plot.')
    parser.add_argument('--variation', help='JES/JER variation to consider.')
    args = parser.parse_args()
    return args

def plot_comparison(acc_06Jan, acc_19Feb, ptag, regex, region, distribution, var=''):
    '''
    Given two accumulators containing the coffea files with 06Jan20 and 19Feb20 inputs,
    plot the specified distributions for comparison.
    ==============
    PARAMETERS:
    ==============
    acc_06Jan    : Accumulator containing coffea files with 06Jan20 inputs. 
    acc_19Feb    : Accumulator containing coffea files with 19Feb20 inputs. 
    ptag         : Tag for the physics process.
    regex        : Regular expression matching dataset name.
    region       : The region from where the data will be taken.
    distribution : Distribution to plot.
    var          : Variation for the comparison to be made, by default nominal case will be compared.
    '''
    acc_06Jan.load(distribution)
    h_06Jan = acc_06Jan[distribution]

    acc_19Feb.load(distribution)
    h_19Feb = acc_19Feb[distribution]

    def merge_and_scale(h, acc):
        h = merge_extensions(h, acc, reweight_pu=False)
        scale_xs_lumi(h)
        h = merge_datasets(h)
        return h
    
    h_06Jan = merge_and_scale(h_06Jan, acc_06Jan)
    h_19Feb = merge_and_scale(h_19Feb, acc_19Feb)

    if distribution in REBIN.keys():
        new_bin = REBIN[distribution]
        axis_name = 'jetpt' if 'ak4_pt' in distribution else distribution
        h_06Jan = h_06Jan.rebin(h_06Jan.axis(axis_name), new_bin)
        h_19Feb = h_19Feb.rebin(h_19Feb.axis(axis_name), new_bin)
    
    histograms = {
        '06Jan' : h_06Jan,
        '19Feb' : h_19Feb
    }

    sumw = {}
    sumw2 = {}
    # Store bin edges and centers for histogram plotting
    edges = h_06Jan.axis(axis_name).edges(overflow='over')
    centers = h_06Jan.axis(axis_name).centers(overflow='over')

    # Pick the relevant dataset and region
    for tag, h in histograms.items():
        h = h[re.compile(regex)].integrate('dataset')
        h = h.integrate('region', re.compile(f'{region}.*'))
        # Get the relevant variation (by default, it is nominal)
        h = h.integrate('var', var)

        # Store the values in new dict
        sumw[tag], sumw2[tag] = h.values(overflow='over', sumw2=True)[()]
    
    # Calculate and store the ratios
    ratios = sumw['19Feb'] / sumw['06Jan']

    # Gaussian error propagation on the ratios
    uncs = np.hypot(
        np.sqrt(sumw2['19Feb'])/sumw['06Jan'],
        np.sqrt(sumw2['06Jan'])/sumw['19Feb']
        )
    
    # Plot the comparison
    fig, (ax, rax) = plt.subplots(2, 1, figsize=(7,7), gridspec_kw={"height_ratios": (2,1)}, sharex=True)
    for tag, val in sumw.items():
        ax.step(edges[:-1], val, label=tag, where='post')
    
    ax.legend()
    ax.set_xlim(edges[0], edges[-2]) # Do not include overflow bin

    dataset_name = regex.replace('.*', '_HT_')
    ax.set_title(dataset_name)

    # Ratio pad
    rax.errorbar(x=centers, y=ratios, yerr=uncs, ls='', marker='o', color='k')
    rax.set_xlabel(h_06Jan.axis(axis_name).label)
    rax.set_ylabel('19 Feb / 06 Jan')
    rax.set_ylim(0.8,1.2)
    rax.grid(True)

    rax.plot(rax.get_xlim(), [1.0, 1.0], 'r--')

    # Save plot
    if var == '':
        outdir = f'./output/input_comparisons'
    else:
        var_label = var.replace('_', '', 1)
        outdir = f'./output/input_comparisons/{var_label}'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    outpath = pjoin(outdir, f'{ptag}_{distribution}_comparison.pdf')
    fig.savefig(outpath)
    plt.close()
    print(f'File saved: {outpath}')

def main():
    args = parse_commandline()
    # Get the specific distribution, if requested
    # Otherwise, just loop over all distributions we are interested in
    if args.distribution:
        distributions = [args.distribution]
    else:
        distributions = ['mjj', 'recoil', 'ak4_pt', 'ak4_pt0', 'ak4_pt1']

    # List of dataset regular expressions and tags
    # No need for 2018 samples as the nanoaod-tools fix was about 2017 only
    dataset_info = {
        'wjets17' : {'regex': 'WJetsToLNu.*2017', 'region': 'sr_vbf'},
        'zjets17' : {'regex': 'ZJetsToNuNu.*2017', 'region': 'sr_vbf'},
        'gjets17' : {'regex': 'GJets_DR-0p4.*2017', 'region': 'cr_g_vbf'}
    }

    # Use the specific variation if specified
    # Otherwise, just look at the nominal case
    var = ''
    if args.variation:
        var = args.variation

    acc_06Jan = dir_archive(
        args.path_to_06Jan,
        serialized=True,
        memsize=1e3,
        compression=0
    )

    acc_06Jan.load('sumw')
    acc_06Jan.load('sumw2')

    acc_19Feb = dir_archive(
        args.path_to_19Feb,
        serialized=True,
        memsize=1e3,
        compression=0
    )

    acc_19Feb.load('sumw')
    acc_19Feb.load('sumw2')

    for dist in distributions:
        for ptag, info in dataset_info.items():
            regex, region = info.values()
            plot_comparison(acc_06Jan, acc_19Feb, ptag=ptag, regex=regex, region=region, distribution=dist, var=var)

if __name__ == '__main__':
    main()

