#!/usr/bin/env python

import os
import sys
import re
import uproot
import warnings
import numpy as np

from matplotlib import pyplot as plt
from coffea import hist
from bucoffea.plot.util import merge_datasets, merge_extensions, scale_xs_lumi, fig_ratio
from klepto.archives import dir_archive
from pprint import pprint

pjoin = os.path.join

warnings.filterwarnings('ignore')

def get_pretty_dataset_tag(dataset):
    mapping = {
        'ZJetsToNuNu' : r'QCD $Z(\nu\nu)$',
        'WJetsToLNu' : r'QCD $W(\ell\nu)$',
        'DYJetsToLL' : r'QCD $Z(\ell\ell)$',
    }
    return mapping[dataset]

def preprocess(h, acc, dataset, region, year):
    h = merge_extensions(h, acc)
    scale_xs_lumi(h)
    h = merge_datasets(h)

    # Rebin mjj
    mjj_ax = hist.Bin('mjj', r'$M_{jj} \ (GeV)$', [200., 400., 600., 900., 1200., 1500., 2000., 2750., 3500., 5000.])
    h = h.rebin('mjj', mjj_ax)

    return h.integrate('region', region).integrate('dataset', re.compile(f'{dataset}.*{year}'))

def get_btag_variations(acc, outtag, dataset='ZJetsToNuNu', region='sr_vbf_no_veto_all', year=2017, outputrootfile=None):
    '''Get b-tag weight variations as a function of mjj.'''
    distributions = {'central' : 'mjj', 'up' : 'mjj_bveto_up', 'down' : 'mjj_bveto_down'}
    histograms = {}
    for k, d in distributions.items():
        acc.load(d)
        histograms[k] = preprocess(acc[d], acc, dataset=dataset, region=region, year=year)

    fig, ax, rax = fig_ratio()
    hist.plot1d(histograms['central'], ax=ax)
    hist.plot1d(histograms['up'], ax=ax, clear=False)
    hist.plot1d(histograms['down'], ax=ax, clear=False)

    ax.set_yscale('log')
    ax.set_ylim(1e-2,1e6)

    ax.legend(title='b-tag Weights', labels=['Central', 'Up', 'Down'])

    ax.text(0.,1.,get_pretty_dataset_tag(dataset),
        fontsize=14,
        ha='left',
        va='bottom',
        transform=ax.transAxes
    )

    ax.text(1.,1.,year,
        fontsize=14,
        ha='right',
        va='bottom',
        transform=ax.transAxes
    )

    data_err_opts = {
        'linestyle':'none',
        'marker': '.',
        'markersize': 10.,
    }

    hist.plotratio(
        histograms['up'],        
        histograms['central'],
        ax=rax,
        unc='num',
        label='Up',
        error_opts=data_err_opts        
    )

    hist.plotratio(
        histograms['down'],        
        histograms['central'],
        ax=rax,
        unc='num',
        label='Down',
        error_opts=data_err_opts,
        clear=False
    )

    rax.grid(True)
    rax.set_ylim(0.9,1.1)
    rax.set_ylabel('Variation / Nominal')
    rax.legend()

    outdir = f'./output/{outtag}'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outpath = pjoin(outdir, f'btag_variations_{dataset}_{year}.pdf')
    fig.savefig(outpath)
    plt.close(fig)

    print(f'File saved: {outpath}')

    # Save the variations to ROOT
    xedges = histograms['central'].axis('mjj').edges()
    r_up = histograms['up'].values()[()] / histograms['central'].values()[()]
    r_down = histograms['down'].values()[()] / histograms['central'].values()[()]

    outputrootfile[f'{dataset}_{year}_btagUp'] = (r_up, xedges)
    outputrootfile[f'{dataset}_{year}_btagDown'] = (r_down, xedges)

def main():
    inpath = sys.argv[1]
    acc = dir_archive(inpath)
    acc.load('sumw')
    acc.load('sumw_pileup')
    acc.load('nevents')

    outtag = re.findall('merged_.*', inpath)[0].replace('/','')

    datasets = [
        'ZJetsToNuNu',
        'WJetsToLNu',
        'DYJetsToLL',
    ]

    # Output ROOT file to save the variations per dataset
    outdir = f'./output/{outtag}'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    outputrootfile = uproot.recreate( pjoin(outdir, 'btag_variations.root') )


    for year in [2017, 2018]:
        for dataset in datasets:
            get_btag_variations(acc,
                outtag=outtag,
                dataset=dataset,
                year=year,
                outputrootfile=outputrootfile
            )

if __name__ == '__main__':
    main()