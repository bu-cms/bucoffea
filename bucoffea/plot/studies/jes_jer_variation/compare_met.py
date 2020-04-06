#!/usr/bin/env python

#####################
# Script to compare several types of MET in NanoAOD:
# JER smeared MET pt: met_pt_jer
# MET pt with no JER smearing applied: met_pt_nom
# MET pt with JES / JER variations applied

# INPUT: Use with 2020-03-28_vbf_jes_jer_var job
######################

import os
import sys
import re
import numpy as np
from bucoffea.plot.util import merge_extensions, merge_datasets, scale_xs_lumi
from coffea import hist
from klepto.archives import dir_archive
from matplotlib import pyplot as plt
from pprint import pprint

pjoin = os.path.join

def preprocess_histos(h, acc, regex, integrate_region=True):
    '''Merging, scaling, integrating.'''
    h = merge_extensions(h, acc, reweight_pu=False)
    scale_xs_lumi(h)
    h = merge_datasets(h)
    # Histograms filled only for signal region, integrate it out
    if integrate_region:
        h = h.integrate('region').integrate('dataset', re.compile(regex) )
    else:
        h = h.integrate('dataset', re.compile(regex) )

    # Rebin
    new_met_bin = hist.Bin('met',r'$p_{T}^{miss}$ (GeV)',list(range(0,500,50)) + list(range(500,1100,100)) )
    h = h.rebin('met', new_met_bin)
    return h

def compare_jer_nom_met(acc, regex, dataset_name, tag, outtag, inclusive=True):
    '''Plot JER smeared and non-JER smeared MET pt.'''
    inc_suffix = '_inc' if inclusive else ''
    acc.load(f'met_jer{inc_suffix}')
    acc.load(f'met_nom{inc_suffix}')

    # Get pre-processed histograms
    histograms = {
        'with JER' : preprocess_histos(acc[f'met_jer{inc_suffix}'], acc, regex),
        'without JER' : preprocess_histos(acc[f'met_nom{inc_suffix}'], acc, regex)
    }

    met_edges = histograms['with JER'].axes()[0].edges()

    # Get histogram values as numpy arrays
    values = {key : histograms[key].values()[()] for key in histograms.keys()}

    # Plot comparison between the two
    fig, ax = plt.subplots(1,1)
    for key, arr in values.items():
        ax.step(met_edges[:-1], arr, label=key, where='post')

    ax.legend(title='Signal Region VBF')
    ax.set_xlim(met_edges[0], met_edges[-2])
    ax.set_xlabel(r'$p_T^{miss}$ (GeV)')
    ax.set_ylabel('Counts')
    title = f'{dataset_name} (Inclusive)' if inclusive else f'{dataset_name} (VBF selection)' 
    ax.set_title(title)

    # Save figure
    outdir = f'./output/{outtag}/met_comparison'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    filename = f'{tag}_nom_jer_met_comparison{inc_suffix}.pdf'

    outpath = pjoin(outdir, filename)
    fig.savefig(outpath)
    print(f'Figure saved: {outpath}')

def plot_varied_met(acc, regex, region, dataset_name, tag, outtag, inclusive=True):
    '''Plot JES/JER varied MET.'''
    inc_suffix = '_inc' if inclusive else ''
    dist = f'met{inc_suffix}'
    acc.load(dist)
    h = preprocess_histos(acc[dist], acc, regex, integrate_region=False)
    h = h[re.compile(f'{region}.*')]

    fig, ax = plt.subplots(1,1)
    hist.plot1d(h, ax=ax, overlay='region')
    title = f'{dataset_name} (Inclusive)' if inclusive else f'{dataset_name} (VBF selection)' 
    ax.set_title(title)

    # Save figure
    outdir = f'./output/{outtag}/met_comparison'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    filename = f'{tag}_varied_met{inc_suffix}.pdf'
    outpath = pjoin(outdir, filename)
    fig.savefig(outpath)
    
    print(f'Figure saved: {outpath}')
    plt.close()

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
        outtag = inpath.split('/')[-2]
    else:
        outtag = inpath.split('/')[-1]

    data_info = {
        # W+jets processes
        'wjets_sr_2017' : {'dataset_name' : 'WJetsToLNu_HT_2017', 'regex' : 'WJetsToLNu.*2017', 'region' : 'sr_vbf'},
        'wjets_sr_2018' : {'dataset_name' : 'WJetsToLNu_HT_2018', 'regex' : 'WJetsToLNu.*2018', 'region' : 'sr_vbf'},
        'wtoenu_2017' : {'dataset_name' : 'WJetsToLNu_HT_2017', 'regex' : 'WJetsToLNu.*2017', 'region' : 'cr_1e_vbf'},
        'wtoenu_2018' : {'dataset_name' : 'WJetsToLNu_HT_2018', 'regex' : 'WJetsToLNu.*2018', 'region' : 'cr_1e_vbf'},
        'wtomunu_2017' : {'dataset_name' : 'WJetsToLNu_HT_2017', 'regex' : 'WJetsToLNu.*2017', 'region' : 'cr_1m_vbf'},
        'wtomunu_2018' : {'dataset_name' : 'WJetsToLNu_HT_2018', 'regex' : 'WJetsToLNu.*2018', 'region' : 'cr_1m_vbf'},

        # Z+jets processes
        'zjets_sr_2017' : {'dataset_name' : 'ZJetsToNuNu_HT_2017', 'regex' : 'ZJetsToNuNu.*2017', 'region' : 'sr_vbf'},
        'zjets_sr_2018' : {'dataset_name' : 'ZJetsToNuNu_HT_2018', 'regex' : 'ZJetsToNuNu.*2018', 'region' : 'sr_vbf'},
        'ztoee_2017' : {'dataset_name' : 'WJetsToLNu_HT_2017', 'regex' : 'ZJetsToNuNu.*2017', 'region' : 'cr_2e_vbf'},
        'ztoee_2018' : {'dataset_name' : 'WJetsToLNu_HT_2018', 'regex' : 'ZJetsToNuNu.*2018', 'region' : 'cr_2e_vbf'},
        'ztomumu_2017' : {'dataset_name' : 'WJetsToLNu_HT_2017', 'regex' : 'ZJetsToNuNu.*2017', 'region' : 'cr_2m_vbf'},
        'ztomumu_2018' : {'dataset_name' : 'WJetsToLNu_HT_2018', 'regex' : 'ZJetsToNuNu.*2018', 'region' : 'cr_2m_vbf'},

        # gamma+jets processes
        'gjets2017' : {'dataset_name' : 'GJets_DR-0p4_HT_2017', 'regex' : 'GJets_DR-0p4.*2017', 'region' : 'cr_g_vbf'},
        'gjets2018' : {'dataset_name' : 'GJets_DR-0p4_HT_2018', 'regex' : 'GJets_DR-0p4.*2018', 'region' : 'cr_g_vbf'}
    }

    for tag, info in data_info.items():
        dataset_name, regex, region = info.values()
        compare_jer_nom_met(acc, dataset_name=dataset_name, regex=regex, tag=tag, outtag=outtag, inclusive=True)
        compare_jer_nom_met(acc, dataset_name=dataset_name, regex=regex, tag=tag, outtag=outtag, inclusive=False)
    
        plot_varied_met(acc, dataset_name=dataset_name, region=region, regex=regex, tag=tag, outtag=outtag, inclusive=True)
        plot_varied_met(acc, dataset_name=dataset_name, region=region, regex=regex, tag=tag, outtag=outtag, inclusive=False)

if __name__ == '__main__':
    main()