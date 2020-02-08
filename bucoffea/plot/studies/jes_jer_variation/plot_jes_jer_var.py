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

# Dict mapping tags to dataset pairs
# and corresponding regexps
tag_to_dataset_pairs = {
    'znunu_over_wlnu17' : {
        'qcd' : {
            'dataset1' : {'regex' : 'ZJetsToNuNu.*2017', 'region' : 'sr_vbf'},
            'dataset2' : {'regex' : 'WJetsToLNu.*2017', 'region' : 'sr_vbf'}
        },
        'ewk' : {
            'dataset1' : {'regex' : 'EWKZ2Jets_ZToNuNu.*2017', 'region' : 'sr_vbf'},
            'dataset2' : {'regex' : 'EWKW2Jets_WToLNu.*2017', 'region' : 'sr_vbf'}
        },
    },
    'znunu_over_wlnu18' : {
        'qcd' : {
            'dataset1' : {'regex' : 'ZJetsToNuNu.*2018', 'region' : 'sr_vbf'},
            'dataset2' : {'regex' : 'WJetsToLNu.*2018', 'region' : 'sr_vbf'}
        },
        'ewk' : {
            'dataset1' : {'regex' : 'EWKZ2Jets_ZToNuNu.*2018', 'region' : 'sr_vbf'},
            'dataset2' : {'regex' : 'EWKW2Jets_WToLNu.*2018', 'region' : 'sr_vbf'}
        },
    },
    'znunu_over_zmumu17' : {
        'qcd': {
            'dataset1' : {'regex' : 'ZJetsToNuNu.*2017', 'region' : 'sr_vbf'},
            'dataset2' : {'regex' : 'DYJetsToLL.*2017', 'region' : 'cr_2m_vbf'},
        },
        'ewk': {
            'dataset1' : {'regex' : 'EWKZ2Jets_ZToNuNu.*2017', 'region' : 'sr_vbf'},
            'dataset2' : {'regex' : 'EWKZ2Jets_ZToLL.*2017', 'region' : 'cr_2m_vbf'},
        },
    },
    'znunu_over_zmumu18' : {
        'qcd': {
            'dataset1' : {'regex' : 'ZJetsToNuNu.*2018', 'region' : 'sr_vbf'},
            'dataset2' : {'regex' : 'DYJetsToLL.*2018', 'region' : 'cr_2m_vbf'},
        },
        'ewk': {
            'dataset1' : {'regex' : 'EWKZ2Jets_ZToNuNu.*2018', 'region' : 'sr_vbf'},
            'dataset2' : {'regex' : 'EWKZ2Jets_ZToLL.*2018', 'region' : 'cr_2m_vbf'},
        },
    },
    'znunu_over_zee17' : {
        'qcd': {
            'dataset1' : {'regex' : 'ZJetsToNuNu.*2017', 'region' : 'sr_vbf'},
            'dataset2' : {'regex' : 'DYJetsToLL.*2017', 'region' : 'cr_2e_vbf'},
        },
        'ewk': {
            'dataset1' : {'regex' : 'EWKZ2Jets_ZToNuNu.*2017', 'region' : 'sr_vbf'},
            'dataset2' : {'regex' : 'EWKZ2Jets_ZToLL.*2017', 'region' : 'cr_2e_vbf'},
        },
    },
    'znunu_over_zee18' : {
        'qcd': {
            'dataset1' : {'regex' : 'ZJetsToNuNu.*2018', 'region' : 'sr_vbf'},
            'dataset2' : {'regex' : 'DYJetsToLL.*2018', 'region' : 'cr_2e_vbf'},
        },
        'ewk': {
            'dataset1' : {'regex' : 'EWKZ2Jets_ZToNuNu.*2018', 'region' : 'sr_vbf'},
            'dataset2' : {'regex' : 'EWKZ2Jets_ZToLL.*2018', 'region' : 'cr_2e_vbf'},
        },
    },
    'wlnu_over_wenu17' : {
        'qcd': {
            'dataset1' : {'regex' : 'WJetsToLNu.*2017', 'region' : 'sr_vbf'},
            'dataset2' : {'regex' : 'WJetsToLNu.*2017', 'region' : 'cr_1e_vbf'},
        },
        'ewk': {
            'dataset1' : {'regex' : 'EWKW2Jets_WToLNu.*2017', 'region' : 'sr_vbf'},
            'dataset2' : {'regex' : 'EWKW2Jets_WToLNu.*2017', 'region' : 'cr_1e_vbf'},
        },
    },
    'wlnu_over_wenu18' : {
        'qcd': {
            'dataset1' : {'regex' : 'WJetsToLNu.*2018', 'region' : 'sr_vbf'},
            'dataset2' : {'regex' : 'WJetsToLNu.*2018', 'region' : 'cr_1e_vbf'},
        },
        'ewk': {
            'dataset1' : {'regex' : 'EWKW2Jets_WToLNu.*2018', 'region' : 'sr_vbf'},
            'dataset2' : {'regex' : 'EWKW2Jets_WToLNu.*2018', 'region' : 'cr_1e_vbf'},
        },
    },
    'wlnu_over_wmunu17' : {
        'qcd': {
            'dataset1' : {'regex' : 'WJetsToLNu.*2017', 'region' : 'sr_vbf'},
            'dataset2' : {'regex' : 'WJetsToLNu.*2017', 'region' : 'cr_1m_vbf'},
        },
        'ewk': {
            'dataset1' : {'regex' : 'EWKW2Jets_WToLNu.*2017', 'region' : 'sr_vbf'},
            'dataset2' : {'regex' : 'EWKW2Jets_WToLNu.*2017', 'region' : 'cr_1m_vbf'},
        },
    },
    'wlnu_over_wmunu18' : {
        'qcd': {
            'dataset1' : {'regex' : 'WJetsToLNu.*2018', 'region' : 'sr_vbf'},
            'dataset2' : {'regex' : 'WJetsToLNu.*2018', 'region' : 'cr_1m_vbf'},
        },
        'ewk': {
            'dataset1' : {'regex' : 'EWKW2Jets_WToLNu.*2018', 'region' : 'sr_vbf'},
            'dataset2' : {'regex' : 'EWKW2Jets_WToLNu.*2018', 'region' : 'cr_1m_vbf'},
        },
    },
}

# Dict mappping dataset names to 
# corresponding regexps
dataset_regex = {
    ### W processes 
    'wlnu17'  : {
        'qcd' : {'title': r'QCD $W\rightarrow \ell \nu$', 'regex': 'WJetsToLNu.*2017', 'region': 'sr_vbf'},
        'ewk' : {'title': r'EWK $W\rightarrow \ell \nu$', 'regex': 'EWKW2Jets_WToLNu.*2017', 'region': 'sr_vbf'}
    },
    'wlnu18'  : {
        'qcd' : {'title': r'QCD $W\rightarrow \ell \nu$', 'regex': 'WJetsToLNu.*2018', 'region': 'sr_vbf'},
        'ewk' : {'title': r'EWK $W\rightarrow \ell \nu$', 'regex': 'EWKW2Jets_WToLNu.*2018', 'region': 'sr_vbf'}
    },
    'wenu17'  : {
        'qcd' : {'title': r'QCD $W\rightarrow e \nu$', 'regex': 'WJetsToLNu.*2017', 'region': 'cr_1e_vbf'},
        'ewk' : {'title': r'EWK $W\rightarrow e \nu$', 'regex': 'EWKW2Jets_WToLNu.*2017', 'region': 'cr_1e_vbf'}
    },
    'wenu18'  : {
        'qcd' : {'title': r'QCD $W\rightarrow e \nu$', 'regex': 'WJetsToLNu.*2018', 'region': 'cr_1e_vbf'},
        'ewk' : {'title': r'EWK $W\rightarrow e \nu$', 'regex': 'EWKW2Jets_WToLNu.*2018', 'region': 'cr_1e_vbf'}
    },
    'wmunu17'  : {
        'qcd' : {'title': r'QCD $W\rightarrow \mu \nu$', 'regex': 'WJetsToLNu.*2017', 'region': 'cr_1m_vbf'},
        'ewk' : {'title': r'EWK $W\rightarrow \mu \nu$', 'regex': 'EWKW2Jets_WToLNu.*2017', 'region': 'cr_1m_vbf'}
    },
    'wmunu18'  : {
        'qcd' : {'title': r'QCD $W\rightarrow \mu \nu$', 'regex': 'WJetsToLNu.*2018', 'region': 'cr_1m_vbf'},
        'ewk' : {'title': r'EWK $W\rightarrow \mu \nu$', 'regex': 'EWKW2Jets_WToLNu.*2018', 'region': 'cr_1m_vbf'}
    },
    ### Z processes
    'zmumu17'  : {
        'qcd' : {'title': r'QCD $Z\rightarrow \mu \mu$', 'regex': 'DYJetsToLL.*2017', 'region': 'cr_2m_vbf'},
        'ewk' : {'title': r'EWK $Z\rightarrow \mu \mu$', 'regex': 'EWKZ2Jets_ZToLL.*2017', 'region': 'cr_2m_vbf'}
    },
    'zmumu18'  : {
        'qcd' : {'title': r'QCD $Z\rightarrow \mu \mu$', 'regex': 'DYJetsToLL.*2018', 'region': 'cr_2m_vbf'},
        'ewk' : {'title': r'EWK $Z\rightarrow \mu \mu$', 'regex': 'EWKZ2Jets_ZToLL.*2018', 'region': 'cr_2m_vbf'}
    },
    'zee17'  : {
        'qcd' : {'title': r'QCD $Z\rightarrow ee$', 'regex': 'DYJetsToLL.*2017', 'region': 'cr_2e_vbf'},
        'ewk' : {'title': r'EWK $Z\rightarrow ee$', 'regex': 'EWKZ2Jets_ZToLL.*2017', 'region': 'cr_2e_vbf'}
    },
    'zee18'  : {
        'qcd' : {'title': r'QCD $Z\rightarrow ee$', 'regex': 'DYJetsToLL.*2018', 'region': 'cr_2e_vbf'},
        'ewk' : {'title': r'EWK $Z\rightarrow ee$', 'regex': 'EWKZ2Jets_ZToLL.*2018', 'region': 'cr_2e_vbf'}
    },
    'znunu17'  : {
        'qcd' : {'title': r'QCD $Z\rightarrow \nu \nu$', 'regex': 'ZJetsToNuNu.*2017', 'region': 'sr_vbf'},
        'ewk' : {'title': r'EWK $Z\rightarrow \nu \nu$', 'regex': 'EWKZ2Jets_ZToNuNu.*2017', 'region': 'sr_vbf'}
    },
    'znunu18'  : {
        'qcd' : {'title': r'QCD $Z\rightarrow \nu \nu$', 'regex': 'ZJetsToNuNu.*2018', 'region': 'sr_vbf'},
        'ewk' : {'title': r'EWK $Z\rightarrow \nu \nu$', 'regex': 'EWKZ2Jets_ZToNuNu.*2018', 'region': 'sr_vbf'}
    },
    
    
}

def dict_to_arr(d):
    '''Given a dictionary containing different weights as
       its values, concatenate the weights as a 2D numpy array.'''
    shape = len(d), len(list(d.values())[0])
    arr = np.zeros(shape)
    for idx, weight_arr in enumerate(d.values()):
        arr[idx] = weight_arr
    return arr

def get_unc(d, edges, out_tag, tag):
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
        
    # Dump the results to a .txt file
    outdir = f'./output/{out_tag}/txt'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outpath = pjoin(outdir, f'{tag}_jes_jer_unc.txt')
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
       with different JES/JER variations in the same canvas.'''
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
    h = h.rebin('mjj', mjj_bins)

    h = h[re.compile(f'{region}')].integrate('region')

    # Calculate the ratios of each variation
    # with respect to nominal counts
    h_nom = h.integrate('var', '').values(overflow='over')[()]
    ratios = {}
    for variation in h.identifiers('var'):
        ratios[variation.name] = h.integrate('var', variation).values(overflow='over')[()] / h_nom

    mjj_edges = h.axes()[0].edges(overflow='over')
    mjj_centers = ((mjj_edges + np.roll(mjj_edges, -1))/2)[:-1]
    
    # Plot the variation + ratio pad
    fig, (ax, rax) = plt.subplots(2, 1, figsize=(7,7), gridspec_kw={"height_ratios": (3, 2)}, sharex=True)
    for idx, (var, ratio_arr) in enumerate(ratios.items()):
        h_var = h.integrate('var', var).values(overflow='over')[()]
        hep.histplot(h_var, 
                     mjj_edges, 
                     label=var_to_legend_label[var],
                     ax=ax,
                     stack=True,
                     histtype='step'
                     )

        if var != '':
            r = h_var / h_nom
            rax.plot(mjj_centers, r, 'o', label=var_to_legend_label[var], c=colors[idx])

    # Aesthetics
    ax.set_xlim(200,4000)
    ax.set_ylabel('Counts / Bin Width')
    ax.set_title(title)
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

    outpath = pjoin(outdir, f'{tag}_jes_jer_variations.pdf')
    fig.savefig(outpath)
    plt.close()
    
    print(f'Histogram saved in: {outpath}')

def plot_jes_jer_var_ratio(acc, regex1, regex2, region1, region2, tag, out_tag, sample_type):
    '''Given the input accumulator, plot ratio of two datasets
       for each JES/JER variation, on the same canvas.'''
    acc.load('mjj')
    h = acc['mjj']

    print(f'Working on: {tag}, {sample_type}')

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
        ratios[var[0]] = h1_vals[var] / h2_vals[var]

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
        'wlnu_over_wenu17' : r'{} $W\rightarrow \ell \nu$ SR / {} $W\rightarrow e\nu$ CR'.format(sample_label, sample_label),
        'wlnu_over_wenu18' : r'{} $W\rightarrow \ell \nu$ SR / {} $W\rightarrow e\nu$ CR'.format(sample_label, sample_label),
        'wlnu_over_wmunu17' : r'{} $W\rightarrow \ell \nu$ SR / {} $W\rightarrow \mu \nu$ CR'.format(sample_label, sample_label),
        'wlnu_over_wmunu18' : r'{} $W\rightarrow \ell \nu$ SR / {} $W\rightarrow \mu \nu$ CR'.format(sample_label, sample_label)
    }
    
    # Upper y-limits for each tag
    tag_to_ylim = {
        'znunu_over_wlnu17' : 2.5, 
        'znunu_over_wlnu18' : 2.5, 
        'znunu_over_zmumu17' : 20,
        'znunu_over_zmumu18' : 20,
        'znunu_over_zee17' : 20,
        'znunu_over_zee18' : 20,
        'wlnu_over_wenu17' : 2,
        'wlnu_over_wenu18' : 2,
        'wlnu_over_wmunu17' : 2,
        'wlnu_over_wmunu18' : 2,
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
                     histtype='step'
                     )

        h1_ = h1.integrate('var', var)
        h2_ = h2.integrate('var', var)

        if var != '':
            r = ratios[var] / ratios['']
            rax.plot(mjj_centers, r, 'o', label=var_to_legend_label[var], c=colors[idx])

    # Aesthetics
    ax.set_xlim(200,4000)
    ax.set_ylim(0,tag_to_ylim[tag])
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

    rax.set_ylim(0., 2.)
    loc = matplotlib.ticker.MultipleLocator(base=0.2)
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
    get_unc(ratios, mjj_edges, out_tag, tag)

def main():
    inpath, sample_type = sys.argv[1:]

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
    

    for tag, data_dict in dataset_regex.items():
        title, regex, region = data_dict[sample_type].values()
        plot_jes_jer_var(acc, regex=regex, title=title, tag=tag, out_tag=out_tag, region=region, sample_type=sample_type)


    for tag, data_dict in tag_to_dataset_pairs.items():
        # Choose QCD or EWK
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

if __name__ == '__main__':
    main()
