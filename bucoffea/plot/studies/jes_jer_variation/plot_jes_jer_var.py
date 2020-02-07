#!/usr/bin/env python

import os
import sys
import re
from bucoffea.plot.util import merge_datasets, merge_extensions, scale_xs_lumi
from klepto.archives import dir_archive
from coffea import hist
from matplotlib import pyplot as plt
import mplhep as hep
import numpy as np

pjoin = os.path.join

def dict_to_arr(d):
    '''Given a dictionary containing different weights as
       its values, concatenate the weights as a 2D numpy array.'''
    shape = len(d), len(list(d.values())[0])
    arr = np.zeros(shape)
    for idx, weight_arr in enumerate(d.values()):
        arr[idx] = weight_arr
    return arr

def get_unc(d, edges):
    '''Given a dictionary containing different weights as
       its values, calculate the uncertainty in each bin.'''
    # Transform to 2D array
    arr = dict_to_arr(d)
    nom = arr[0] 
    # Calculate uncertainty in each mjj bin
    unc = np.zeros_like(nom)
    for idx, entry in enumerate(nom):
        bin_ = arr[:, idx]
        rng = bin_.max() - bin_.min()
        unc[idx] = rng/nom[idx]
        
    # Print the results
    print('*'*20)
    print('Combined Uncertainties')
    print('*'*20)
    for idx, sigma in enumerate(unc):
        print(f'{edges[idx]:.0f} < mjj < {edges[idx+1]:.0f}: {sigma*100:.2f}%')

def plot_jes_jer_var(acc, regex, tag, out_tag, dataset_name):
    '''Given the input accumulator and the regex
       describing the dataset, plot the mjj distribution
       with different JES/JER variations in the same canvas.'''
    acc.load('mjj')
    h = acc['mjj']

    print(f'Working on: {tag}')

    h = merge_extensions(h, acc, reweight_pu=False)
    scale_xs_lumi(h)
    h = merge_datasets(h)
    
    # Pick the relevant dataset
    h = h[re.compile(regex)].integrate('dataset')

    # Rebin mjj
    mjj_bins = hist.Bin('mjj', r'$M_{jj}$ (GeV)', list(range(200,800,300)) + list(range(800,2000,400)) + [2000, 2750, 3500])
    h = h.rebin('mjj', mjj_bins)

    h = h[re.compile('sr_vbf.*')].integrate('var')

    # Plot the variation
    fig, ax = plt.subplots(1,1)
    hist.plot1d(h, ax=ax, overlay='region', binwnorm=True)
    ax.set_ylabel('Counts / Bin Width')
    ax.set_title(dataset_name)

    # Save figure
    outdir = f'./output/{out_tag}'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    outpath = pjoin(outdir, f'{tag}_jes_jer_variations.pdf')
    fig.savefig(outpath)
    
    print(f'Histogram saved in: {outpath}')

def plot_jes_jer_var_ratio(acc, regex1, regex2, region1, region2, tag, out_tag):
    '''Given the input accumulator, plot ratio of two datasets
       for each JES/JER variation, on the same canvas.'''
    acc.load('mjj')
    h = acc['mjj']

    print(f'Working on: {tag}')

    h = merge_extensions(h, acc, reweight_pu=False)
    scale_xs_lumi(h)
    h = merge_datasets(h)
    
    # Pick the relevant dataset
    h = h[re.compile(f'{regex1}|{regex2}')]
    
    # Rebin mjj
    mjj_bins = hist.Bin('mjj', r'$M_{jj}$ (GeV)', list(range(200,800,300)) + list(range(800,2000,400)) + [2000, 2750, 3500])
    h = h.rebin('mjj', mjj_bins)

    h1 = h[re.compile(regex1)].integrate('dataset')
    h2 = h[re.compile(regex2)].integrate('dataset')

    h1 = h1[re.compile(f'{region1}.*')].integrate('region')
    h2 = h2[re.compile(f'{region2}.*')].integrate('region')

    # Calculate the ratios for each JES/JER variation
    ratios = {}
    h1_vals = h1.values(overflow='over')
    h2_vals = h2.values(overflow='over')
    
    for var in h1_vals.keys():
        ratios[var] = h1_vals[var] / h2_vals[var]
    
    # Plot the ratios for each variation
    var_to_legend_label = {
        ''         : 'Nominal',
        '_jerup'   : 'JER up',
        '_jerdown' : 'JER down',
        '_jesup'   : 'JES up',
        '_jesdown' : 'JES down'
    }

    # The y-axis labels for each tag
    tag_to_ylabel = {
        'z_over_w17' : r'$Z\rightarrow \nu \nu$ SR / $W\rightarrow \ell \nu$ SR (2017)',
        'z_over_w18' : r'$Z\rightarrow \nu \nu$ SR / $W\rightarrow \ell \nu$ SR (2018)'
    }
    
    # Upper y-limits for each tag
    tag_to_ylim = {
        'z_over_w17' : 2.5, 
        'z_over_w18' : 2.5, 
    }
   
    mjj_edges = h1.axes()[0].edges(overflow='over')
    
    fig, ax = plt.subplots(1,1)
    for var, ratio_arr in ratios.items():
        hep.histplot(ratio_arr, 
                     mjj_edges, 
                     label=var_to_legend_label[var[0]],
                     ax=ax,
                     stack=True,
                     histtype='step'
                     )

    ax.set_xlim(200,4000)
    ax.set_ylim(0,tag_to_ylim[tag])
    ax.set_xlabel(r'$M_{jj} \ (GeV)$')
    ax.set_ylabel(tag_to_ylabel[tag])

    plt.legend()

    # Save the figure
    outdir = f'./output/{out_tag}'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    outpath = pjoin(outdir, f'{tag}_jes_jer_variations.pdf')
    fig.savefig(outpath)

    print(f'Histogram saved in: {outpath}')

    # Calculate and print the uncertainties
    # for each mjj bin
    get_unc(ratios, mjj_edges)

def main():
    inpath = sys.argv[1]

    acc = dir_archive(
                         inpath,
                         serialized=True,
                         compression=0,
                         memsize=1e3
                      )

    acc.load('sumw')
    acc.load('sumw2')

    if inpath.endswith('/'):
        out_tag = inpath.split('/')[-2]
    else:
        out_tag = inpath.split('/')[-1]
    
    # Dict mappping dataset names to 
    # corresponding regexps
    dataset_regex = {
        ('wjet17', 'WJetsToLNu_HT_2017') : 'WJetsToLNu.*2017',
        ('wjet18', 'WJetsToLNu_HT_2018') : 'WJetsToLNu.*2018',
        ('dy17', 'DYJetsToLL_HT_2017') : 'DYJetsToLL.*2017',
        ('dy18', 'DYJetsToLL_HT_2018') : 'DYJetsToLL.*2018',
        ('znunu17', 'ZJetsToNuNu_HT_2017') : 'ZJetsToNuNu.*2017',
        ('znunu18', 'ZJetsToNuNu_HT_2018') : 'ZJetsToNuNu.*2018'
    }

   # for dataset, regex in dataset_regex.items():
   #     dataset_tag, dataset_name = dataset
   #     plot_jes_jer_var(acc, regex=regex, dataset_name=dataset_name, tag=dataset_tag, out_tag=out_tag)

    plot_jes_jer_var_ratio(acc, regex1='ZJetsToNuNu.*2017', regex2='WJetsToLNu.*2017', region1='sr_vbf', region2='sr_vbf', tag='z_over_w17', out_tag=out_tag)
    plot_jes_jer_var_ratio(acc, regex1='ZJetsToNuNu.*2018', regex2='WJetsToLNu.*2018', region1='sr_vbf', region2='sr_vbf', tag='z_over_w18', out_tag=out_tag)

if __name__ == '__main__':
    main()
