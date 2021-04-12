#!/usr/bin/env python

import os
import sys
import re
import warnings
import numpy as np

from coffea import hist
from bucoffea.plot.util import merge_datasets, merge_extensions, scale_xs_lumi, fig_ratio, lumi
from bucoffea.plot.style import matplotlib_rc
from bucoffea.helpers.paths import bucoffea_path
from matplotlib import pyplot as plt
from klepto.archives import dir_archive
from pprint import pprint

pjoin = os.path.join

warnings.filterwarnings('ignore')

matplotlib_rc()

def preprocess(h, acc, dataset, distribution, region):
    h = merge_extensions(h, acc, reweight_pu=False)
    scale_xs_lumi(h)
    h = merge_datasets(h)

    h = h.integrate('region', region).integrate('dataset', dataset)

    if distribution == 'mjj':
        mjj_ax = hist.Bin('mjj', r'$M_{jj} \ (GeV)$', [200., 400., 600., 900., 1200., 1500., 2000., 2750., 3500., 5000.])
        h = h.rebin('mjj', mjj_ax)

    return h

def plot_data_to_data_comparison(acc_dict, outtag, dataset='MET_2018', distribution='mjj', region='sr_vbf'):
    histos = {}
    for key, acc in acc_dict.items():
        acc.load(distribution)
        histos[key] = preprocess(acc[distribution], acc, dataset, distribution, region)

    fig, ax, rax = fig_ratio()
    hist.plot1d(histos['03Sep20v7'], ax=ax)
    hist.plot1d(histos['ULv8'], ax=ax, clear=False)

    ax.set_yscale('log')
    ax.set_ylim(1e0,1e6)

    ax.legend(title='Data Version', labels=[
        'ReReco v7',
        'UL v8',
    ])

    ax.text(0.,1.,f'{dataset.replace("_"," ")}, ${lumi(2018)} \\ fb^{{-1}}$',
        fontsize=14,
        ha='left',
        va='bottom',
        transform=ax.transAxes
    )

    ax.text(1.,1.,'VBF SR (ReReco Cleaning Cuts)',
        fontsize=14,
        ha='right',
        va='bottom',
        transform=ax.transAxes
    )

    data_err_opts = {
        'linestyle':'none',
        'marker': '.',
        'markersize': 10.,
        'color':'k',
    }

    hist.plotratio(
        histos['ULv8'],
        histos['03Sep20v7'],
        ax=rax,
        unc='num',
        error_opts=data_err_opts
    )

    rax.grid(True)
    rax.set_ylabel('v8 / v7')
    rax.set_ylim(0.5,1.5)

    new_xlabels = {
        'ak4_eta0': r'Leading Jet $\eta$',
        'ak4_eta1': r'Trailing Jet $\eta$',
    }

    if distribution in new_xlabels.keys():
        ax.set_xlabel(new_xlabels[distribution])
        rax.set_xlabel(new_xlabels[distribution])

    # Save figure
    outdir = f'./output/{outtag}'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outpath = pjoin(outdir, f'{dataset}_data_to_data_comp_{region}_{distribution}.pdf')
    fig.savefig(outpath)
    plt.close(fig)

    print(f'File saved: {outpath}')

def main():
    acc_dict = {
        '03Sep20v7' : dir_archive(bucoffea_path('submission/merged_2021-04-12_vbfhinv_03Sep20v7_MET_2018_rerecoCleaningCuts')),
        'ULv8' : dir_archive(bucoffea_path('submission/merged_2021-04-12_vbfhinv_ULv8_MET_2018_rerecoCleaningCuts')),
    }
    for acc in acc_dict.values():
        acc.load('sumw')
        acc.load('sumw2')

    outtag = '12Apr21'

    distributions = ['mjj', 'ak4_eta0', 'ak4_eta1']

    for distribution in distributions:
        plot_data_to_data_comparison(acc_dict, outtag, distribution=distribution)

if __name__ == '__main__':
    main()