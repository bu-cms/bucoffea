#!/usr/bin/env python

import os
import sys
import re
import argparse
import warnings
import mplhep as hep
import numpy as np
import uproot

from bucoffea.plot.util import merge_datasets, merge_extensions, scale_xs_lumi, fig_ratio
from coffea import hist
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
from klepto.archives import dir_archive
from pprint import pprint

pjoin = os.path.join

def pretty_legend_label(dataset, region):
    if 'ZJetsToNuNu' in dataset:
        return r'QCD $Z(\nu\nu)$'
    elif 'WJetsToLNu' in dataset:
        if region == 'sr_vbf':
            return r'QCD $W(\ell\nu)$'
        elif region == 'cr_1m_vbf':
            return r'QCD $W(\mu\nu)$'
        else:
            return r'QCD $W(e\nu)$'
    elif 'DYJetsToLL' in dataset:
        if region == 'cr_2m_vbf':
            return r'QCD $Z(\mu\nu)$'
        else:
            return r'QCD $Z(ee)$'
    elif 'GJets' in dataset:
        return r'$\gamma$ + jets'

def plot_jer_uncertainties_in_singlebin(acc, outtag, dataset_regex, region, distribution='mjj'):
    '''Plot JER up and down variations in a single bin (stats!).'''
    acc.load(distribution)
    h = acc[distribution]

    h = merge_extensions(h, acc)
    scale_xs_lumi(h)
    h = merge_datasets(h)

    # Merge everything into a single bin
    mjj_ax = hist.Bin('mjj', r'$M_{jj} \ (GeV)$', [200., 5000.])
    h = h.rebin('mjj', mjj_ax)

    for year in [2017, 2018]:
        _h = h.integrate('dataset', re.compile(f'{dataset_regex}.*{year}'))[re.compile(f'{region}.*')]
        # Nominal histogram
        h_nom = _h.integrate('region', region)
        sumw_nom = h_nom.values()[()]
        xcenters = h_nom.axis('mjj').centers()

        fig, ax = plt.subplots()
        regions = [r for r in map(str, _h.identifiers('region')) if (r != region and 'jesTotal' not in r and 'no_veto_all' not in r)]
        labels = [re.findall('je.*', r)[0] for r in regions]
        ratios = {}
        for reg, label in zip(regions, labels):
            # Calculate the ratio w.r.t. nominal
            ratio = _h.integrate('region', reg).values()[()] / sumw_nom
            opts = {'marker': 'o', 'label': label}
            ax.plot(xcenters,ratio, **opts)

            ratios[label] = ratio

        ax.legend(ncol=1, fontsize='x-small')
        ax.set_xlabel(r'$M_{jj} \ (GeV)$')
        ax.set_ylabel('Ratio to Nominal')
        ax.grid(True)

        ax.text(0.,1.,pretty_legend_label(dataset_regex, region),
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

        outdir = f'./output/{outtag}'
        try:
            os.makedirs(outdir)
        except FileExistsError:
            pass
        
        outpath = pjoin(outdir, f'{dataset_regex.replace(".*", "")}_{region}_{year}.pdf')
        fig.savefig(outpath)
        plt.close(fig)

        print(f'File saved: {outpath}')

def main():
    inpath = sys.argv[1]
    acc = dir_archive(inpath)
    acc.load('sumw')
    acc.load('sumw_pileup')
    acc.load('nevents')

    outtag = re.findall('merged_.*', inpath)[0].replace('/','')

    datasets_regions = {
        'ZJetsToNuNu.*HT.*': ['sr_vbf'],
        'WJetsToLNu_HT.*': ['sr_vbf', 'cr_1m_vbf', 'cr_1e_vbf'],
        'DYJetsToLL_M-50_HT.*': ['cr_2m_vbf', 'cr_2e_vbf'],
        'GJets_HT.*': ['cr_g_vbf'],
    }

    for dataset, regions in datasets_regions.items():
        for region in regions:
            plot_jer_uncertainties_in_singlebin(acc, outtag, dataset_regex=dataset, region=region)

if __name__ == '__main__':
    main()