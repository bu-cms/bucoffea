#!/usr/bin/env python

import os
import sys
import re
from bucoffea.plot.util import merge_datasets, merge_extensions, scale_xs_lumi, fig_ratio
from klepto.archives import dir_archive
from coffea import hist
from matplotlib import pyplot as plt
import matplotlib.ticker
import mplhep as hep
import numpy as np
from data import tag_to_dataset_pairs, dataset_regex

pjoin = os.path.join

# Plot aesthetics
colors = plt.rcParams['axes.prop_cycle'].by_key()['color'][:5]

# Legend labels for each variation
var_to_legend_label = {
    ''         : 'Nominal',
    '_jerup'   : 'JER up',
    '_jerdown' : 'JER down',
    '_jesup'   : 'JES up',
    '_jesdown' : 'JES down'
}

def dict_to_arr(d):
    '''Given a dictionary containing different weights as
       its values, concatenate the weights as a 2D numpy array.'''
    shape = len(d), len(list(d.values())[0])
    arr = np.zeros(shape)
    for idx, weight_arr in enumerate(d.values()):
        arr[idx] = weight_arr
    return arr

def get_unc(d, edges, out_tag, tag, sample_type):
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
        unc[idx] = rng/(2*nom[idx])
        
    # Dump the results to a .txt file
    outdir = f'./output/{out_tag}/txt'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outpath = pjoin(outdir, f'{tag}_{sample_type}_jes_jer_unc.txt')
    with open(outpath, 'w+') as f:
        f.write('*'*20 + '\n')
        f.write('Combined Uncertainties' + '\n')
        f.write('*'*20 + '\n')
        for idx, sigma in enumerate(unc):
            f.write(f'{edges[idx]:.0f} < mjj < {edges[idx+1]:.0f}: {sigma*100:.2f}%\n')

    print(f'File saved: {outpath}')

def plot_jes_jer_var(acc, regex, region, tag, out_tag, title, sample_type):
    '''Given the input accumulator and the regex
    describing the dataset, plot the mjj distribution
    with different JES/JER variations in the same canvas.
    =============
    PARAMETERS:
    =============
    acc         : The input accumulator containing the histograms.
    regex       : Regular expression matching the dataset name.
    region      : The region from which the event yields will be taken.
    tag         : Tag representing the process. (e.g. "wjet")
    out_tag     : Out-tag for output directory naming. The output files are going to be saved
                  under this directory.
    title       : Histogram title for plotting.
    sample_type : QCD ("qcd") or EWK ("ewk") sample.
    '''
    acc.load('mjj')
    h = acc['mjj']

    print(f'Working on: {tag}, {sample_type}')

    h = merge_extensions(h, acc, reweight_pu=False)
    scale_xs_lumi(h)
    h = merge_datasets(h)
    
    # Pick the relevant dataset
    h = h[re.compile(regex)].integrate('dataset')

    # Rebin mjj
    mjj_bins = hist.Bin('mjj', r'$M_{jj}$ (GeV)', list(range(200,800,300)) + list(range(800,2000,400)) + [2000, 2750, 3500])
    mjj_bins_coarse = hist.Bin('mjj', r'$M_{jj}$ (GeV)', list(range(0,4000,1000))) 
    mjj_bins_very_coarse = hist.Bin('mjj', r'$M_{jj}$ (GeV)', [0,4000]) 
    h = h.rebin('mjj', mjj_bins_very_coarse)

    h = h[re.compile(f'{region}')].integrate('region')

    # Calculate the ratios of each variation
    # with respect to nominal counts
    h_nom = h.integrate('var', '').values()[()]
    ratios = {}
    for variation in h.identifiers('var'):
        ratios[variation.name] = h.integrate('var', variation).values()[()] / h_nom

    mjj_edges = h.axes()[0].edges()
    mjj_centers = ((mjj_edges + np.roll(mjj_edges, -1))/2)[:-1]
    
    # Store counts with all variations
    counts = []

    # Plot the variation + ratio pad
    fig, (ax, rax) = plt.subplots(2, 1, figsize=(7,7), gridspec_kw={"height_ratios": (3, 2)}, sharex=True)
    for idx, (var, ratio_arr) in enumerate(ratios.items()):
        h_var = h.integrate('var', var).values()[()]
        hep.histplot(h_var, 
                     mjj_edges, 
                     label=var_to_legend_label[var],
                     ax=ax,
                     stack=True,
                     histtype='step'
                     )

        counts.append(h_var)

        if var != '':
            rax.plot(mjj_centers, ratio_arr, 'o', label=var_to_legend_label[var], c=colors[idx])

    # Aesthetics
    ax.set_xlim(200,4000)
    ax.set_ylabel('Counts / Bin Width')
    ax.set_title(title)
    ax.legend()
    
    # Handle y-limit dynamically
    min_count, max_count = min(counts), max(counts)
    lower_ylim = 0.8 * min_count
    upper_ylim = 1.2 * max_count
    ax.set_ylim(lower_ylim, upper_ylim)

    if '17' in tag:
        ax.text(1., 1., '2017',
                fontsize=12,
                horizontalalignment='right',
                verticalalignment='bottom',
                transform=ax.transAxes
                )
    
    elif '18' in tag:
        ax.text(1., 1., '2018',
                fontsize=12,
                horizontalalignment='right',
                verticalalignment='bottom',
                transform=ax.transAxes
                )

    rax.set_ylim(0., 2.)
    loc = matplotlib.ticker.MultipleLocator(base=0.2)
    rax.yaxis.set_major_locator(loc)
    rax.set_ylabel('Varied / Nominal')
    rax.set_xlabel(r'$M_{jj} \ (GeV)$')
    rax.legend(ncol=2)
    rax.grid(True)

    # Save figure
    outdir = f'./output/{out_tag}'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    outpath = pjoin(outdir, f'{tag}_{sample_type}_jes_jer_variations.pdf')
    fig.savefig(outpath)
    plt.close()
    
    print(f'Histogram saved in: {outpath}')

def plot_jes_jer_var_ratio(acc, regex1, regex2, region1, region2, tag, out_tag, sample_type):
    '''Given the input accumulator, plot ratio of two datasets
    for each JES/JER variation, on the same canvas.
    ==============
    PARAMETERS:
    ==============
    acc         : Input accumulator containing the histograms.
    regex1      : Regular expression matching the dataset name in the numerator of the ratio.
    regex2      : Regular expression matching the dataset name in the denominator of the ratio.
    region1     : The region from which the data for the numerator is going to be taken from.
    region2     : The region from which the data for the denominator is going to be taken from.
    tag         : Tag for the process name. (e.g "wjet")
    out_tag     : Out-tag for naming output directory, output files are going to be saved under this directory.
    sample_type : QCD ("qcd") or EWK ("ewk") sample. 
    '''
    acc.load('mjj')
    h = acc['mjj']

    print(f'Working on: {tag}, {sample_type}')

    h = merge_extensions(h, acc, reweight_pu=False)
    scale_xs_lumi(h)
    h = merge_datasets(h)
    
    # Pick the relevant datasets 
    # Regex 1 matches the dataset for the numerator
    # Regex 2 matches the dataset for the denominator
    h = h[re.compile(f'{regex1}|{regex2}')]
    
    # Rebin mjj
    mjj_bins = hist.Bin('mjj', r'$M_{jj}$ (GeV)', list(range(200,800,300)) + list(range(800,2000,400)) + [2000, 2750, 3500])
    mjj_bins_coarse = hist.Bin('mjj', r'$M_{jj}$ (GeV)', list(range(0,4000,1000))) 
    mjj_bins_very_coarse = hist.Bin('mjj', r'$M_{jj}$ (GeV)', [0,4000]) 
    h = h.rebin('mjj', mjj_bins_very_coarse)

    # h1: Histogram for the numerator
    # h2: Histogram for the denominator
    h1 = h[re.compile(regex1)].integrate('dataset')
    h2 = h[re.compile(regex2)].integrate('dataset')

    h1 = h1[re.compile(f'{region1}.*')].integrate('region')
    h2 = h2[re.compile(f'{region2}.*')].integrate('region')

    # Calculate the ratios and errors in ratios 
    # for each JES/JER variation
    ratios = {}
    err = {}
    h1_vals = h1.values(overflow='over', sumw2=True)
    h2_vals = h2.values(overflow='over', sumw2=True)
    
    for var in h1_vals.keys():
        h1_sumw, h1_sumw2 = h1_vals[var]
        h2_sumw, h2_sumw2 = h2_vals[var]
        # Name of the variation is contained in the first element of the tuple
        varname = var[0]
        ratios[varname] = h1_sumw / h2_sumw 
        # Gaussian error propagation
        gaus_error = np.sqrt((h2_sumw*np.sqrt(h1_sumw2))**2 + (h1_sumw*np.sqrt(h2_sumw2))**2)/h2_sumw**2
        err[varname] = gaus_error

    # Set y-label for either QCD or EWK samples
    sample_label = sample_type.upper()

    # The y-axis labels for each tag
    tag_to_ylabel = {
        'znunu_over_wlnu17' : r'{} $Z\rightarrow \nu \nu$ SR / {} $W\rightarrow \ell \nu$ SR'.format(sample_label, sample_label),
        'znunu_over_wlnu18' : r'{} $Z\rightarrow \nu \nu$ SR / {} $W\rightarrow \ell \nu$ SR'.format(sample_label, sample_label),
        'znunu_over_zmumu17' : r'{} $Z\rightarrow \nu \nu$ SR / {} $Z\rightarrow \mu \mu$ CR'.format(sample_label, sample_label),
        'znunu_over_zmumu18' : r'{} $Z\rightarrow \nu \nu$ SR / {} $Z\rightarrow \mu \mu$ CR'.format(sample_label, sample_label),
        'znunu_over_zee17' : r'{} $Z\rightarrow \nu \nu$ SR / {} $Z\rightarrow ee$ CR'.format(sample_label, sample_label),
        'znunu_over_zee18' : r'{} $Z\rightarrow \nu \nu$ SR / {} $Z\rightarrow ee$ CR'.format(sample_label, sample_label),
        'znunu_over_gjets17' : r'{} $Z\rightarrow \nu \nu$ SR / {} $\gamma$ + jets CR'.format(sample_label, sample_label),
        'znunu_over_gjets18' : r'{} $Z\rightarrow \nu \nu$ SR / {} $\gamma$ + jets CR'.format(sample_label, sample_label),
        'wlnu_over_wenu17' : r'{} $W\rightarrow \ell \nu$ SR / {} $W\rightarrow e\nu$ CR'.format(sample_label, sample_label),
        'wlnu_over_wenu18' : r'{} $W\rightarrow \ell \nu$ SR / {} $W\rightarrow e\nu$ CR'.format(sample_label, sample_label),
        'wlnu_over_wmunu17' : r'{} $W\rightarrow \ell \nu$ SR / {} $W\rightarrow \mu \nu$ CR'.format(sample_label, sample_label),
        'wlnu_over_wmunu18' : r'{} $W\rightarrow \ell \nu$ SR / {} $W\rightarrow \mu \nu$ CR'.format(sample_label, sample_label),
        'wlnu_over_gjets17' : r'{} $W\rightarrow \ell \nu$ SR / {} $\gamma$ + jets CR'.format(sample_label, sample_label),
        'wlnu_over_gjets18' : r'{} $W\rightarrow \ell \nu$ SR / {} $\gamma$ + jets CR'.format(sample_label, sample_label),
    }
    
    # Lower + upper y-limits for each tag
    tag_to_ylim = {
        'znunu_over_wlnu17' : {'qcd' : (1.9, 2.1), 'ewk' : (1.95, 2.15)}, 
        'znunu_over_wlnu18' : {'qcd' : (1.95, 2.15), 'ewk' : (1.95, 2.15)}, 
        'znunu_over_zmumu17' : {'qcd' : (8, 9.5), 'ewk' : (7.6, 8.6)},
        'znunu_over_zmumu18' : {'qcd' : (8, 9.5), 'ewk' : (7.4, 8.8)},
        'znunu_over_zee17' : {'qcd' : (10.5, 11.1), 'ewk' : (0,10)},
        'znunu_over_zee18' : {'qcd' : (10, 10.8), 'ewk' : (0,10)},
        'znunu_over_gjets17' : {'qcd' : (0.65,0.75), 'ewk' : (0,10)},
        'znunu_over_gjets18' : {'qcd' : (0.6,0.68), 'ewk' : (0,10)},
        'wlnu_over_wenu17' : {'qcd' : (0.75,0.85), 'ewk' : (0.15, 0.35)},
        'wlnu_over_wenu18' : {'qcd' : (0.65,0.75), 'ewk' : (0.43, 0.47)},
        'wlnu_over_wmunu17' : {'qcd' : (0.45,0.55), 'ewk' : (0.27, 0.32)},
        'wlnu_over_wmunu18' : {'qcd' : (0.45,0.55), 'ewk' : (0.28, 0.33)},
        'wlnu_over_gjets17' : {'qcd' : (0.3,0.4), 'ewk' : (0.6, 0.8)},
        'wlnu_over_gjets18' : {'qcd' : (0.3,0.4), 'ewk' : (0.95, 1.15)}
    }
   
    mjj_edges = h1.axes()[0].edges(overflow='over')
    mjj_centers = ((mjj_edges + np.roll(mjj_edges, -1))/2)[:-1]

    # Plot the ratios for each variation
    fig, (ax, rax) = plt.subplots(2, 1, figsize=(7,7), gridspec_kw={"height_ratios": (3, 2)}, sharex=True)
    for idx, (var, ratio_arr) in enumerate(ratios.items()):
        hep.histplot(ratio_arr, 
                     mjj_edges, 
                     label=var_to_legend_label[var],
                     ax=ax,
                     stack=True,
                     histtype='step',
                     yerr=err[var]
                     )

        if var != '':
            r = ratios[var] / ratios['']
            rax.plot(mjj_centers, r, 'o', label=var_to_legend_label[var], c=colors[idx])

    # Aesthetics
    ax.set_xlim(200,4000)
    ax.set_ylim(tag_to_ylim[tag][sample_type][0], tag_to_ylim[tag][sample_type][1])
    ax.set_ylabel(tag_to_ylabel[tag])
    ax.legend()

    if '17' in tag:
        ax.text(1., 1., '2017',
                fontsize=12,
                horizontalalignment='right',
                verticalalignment='bottom',
                transform=ax.transAxes
                )
    
    elif '18' in tag:
        ax.text(1., 1., '2018',
                fontsize=12,
                horizontalalignment='right',
                verticalalignment='bottom',
                transform=ax.transAxes
                )

    rax.set_ylim(0.94, 1.06)
    loc = matplotlib.ticker.MultipleLocator(base=0.02)
    rax.yaxis.set_major_locator(loc)
    rax.set_ylabel('Varied / Nominal')
    rax.set_xlabel(r'$M_{jj} \ (GeV)$')
    rax.legend(ncol=2)
    rax.grid(True)

    # Save the figure
    outdir = f'./output/{out_tag}'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    outpath = pjoin(outdir, f'{tag}_{sample_type}_jes_jer_variations.pdf')
    fig.savefig(outpath)
    plt.close()

    print(f'Histogram saved in: {outpath}')

    # Calculate and print the uncertainties
    # for each mjj bin
    get_unc(ratios, mjj_edges, out_tag, tag, sample_type)

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

    sample_types = ['qcd', 'ewk']

    for tag, data_dict in dataset_regex.items():
        for sample_type in sample_types:
            try:
                title, regex, region = data_dict[sample_type].values()
                plot_jes_jer_var(acc, regex=regex, title=title, tag=tag, out_tag=out_tag, region=region, sample_type=sample_type)
            except KeyError:
                continue

    for tag, data_dict in tag_to_dataset_pairs.items():
        for sample_type in sample_types:
            try:
                datapair_dict = data_dict[sample_type] 
                data1_info = datapair_dict['dataset1']
                data2_info = datapair_dict['dataset2']
                plot_jes_jer_var_ratio( acc, 
                                        regex1=data1_info['regex'], 
                                        regex2=data2_info['regex'], 
                                        region1=data1_info['region'], 
                                        region2=data2_info['region'], 
                                        tag=tag, 
                                        out_tag=out_tag,
                                        sample_type=sample_type)
            except KeyError:
                continue

if __name__ == '__main__':
    main()
