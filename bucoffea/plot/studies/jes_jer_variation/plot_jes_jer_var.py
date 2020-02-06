#!/usr/bin/env python

import os
import sys
import re
from bucoffea.plot.util import merge_datasets, merge_extensions, scale_xs_lumi
from klepto.archives import dir_archive
from coffea import hist
from matplotlib import pyplot as plt
import mplhep as hep

pjoin = os.path.join

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

def plot_jes_jer_var_ratio(acc, regex1, regex2, tag, out_tag):
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

    h1 = h1[re.compile('sr_vbf.*')].integrate('var')
    h2 = h2[re.compile('sr_vbf.*')].integrate('var')

    # Calculate the ratios for each JES/JER variation
    ratios = {}
    h1_vals = h1.values(overflow='over')
    h2_vals = h2.values(overflow='over')
    for region in h1_vals.keys():
        h1_weights = h1_vals[region]
        h2_weights = h2_vals[region]
        region_name = region[0]
        ratios[region_name] = h1_weights / h2_weights
    
    # Plot the ratios for each variation
    region_to_legend_label = {
        'sr_vbf'         : 'Nominal',
        'sr_vbf_jerup'   : 'JER up',
        'sr_vbf_jerdown' : 'JER down',
        'sr_vbf_jesup'   : 'JES up',
        'sr_vbf_jesdown' : 'JES down'
    }

    # The y-axis labels for each tag
    tag_to_ylabel = {
        'z_over_w' : r'$Z\rightarrow \ell \ell$ / $W\rightarrow \ell \nu$'
    }
    
    mjj_edges = h1.axes()[1].edges(overflow='over')
    
    fig, ax = plt.subplots(1,1)
    for region, ratio_arr in ratios.items():
        hep.histplot(ratio_arr, 
                     mjj_edges, 
                     label=region_to_legend_label[region],
                     ax=ax,
                     stack=True,
                     histtype='step'
                     )

    ax.set_xlim(200,4000)
    ax.set_ylim(0,0.04)
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

    plot_jes_jer_var_ratio(acc, regex1='DYJetsToLL.*2017', regex2='WJetsToLNu.*2017', tag='z_over_w', out_tag=out_tag)

if __name__ == '__main__':
    main()
