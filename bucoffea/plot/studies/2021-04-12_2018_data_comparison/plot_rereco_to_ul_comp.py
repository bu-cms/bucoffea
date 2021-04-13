#!/usr/bin/env python

import os
import sys
import re
import warnings
import argparse
import numpy as np

from coffea import hist
from bucoffea.plot.util import merge_datasets, merge_extensions, scale_xs_lumi, fig_ratio, lumi
from bucoffea.plot.style import matplotlib_rc
from bucoffea.helpers.paths import bucoffea_path
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
from klepto.archives import dir_archive
from pprint import pprint

pjoin = os.path.join

warnings.filterwarnings('ignore')

matplotlib_rc()

AX_YLIMS = {
    'mjj' : (1e-1,1e5),
    'ak4_eta0' : (1e-2,1e4),
    'ak4_eta1' : (1e-2,1e4),
}

RAX_YLIMS = {
    'mjj' : (0,2),
    'ak4_eta0' : (0,2),
    'ak4_eta1' : (0,2),
}

def parse_cli():
    parser = argparse.ArgumentParser()
    parser.add_argument('outtag', help='Outtag for the job.')
    parser.add_argument('--distribution', default='.*', help='Regex to match the distributions to plot.')
    parser.add_argument('--years', type=int, nargs='*', default=[2017,2018], help='The years to run the comparison on.')
    parser.add_argument('--proc', default='QCD_W', help='The process to run the comparison for.')
    args = parser.parse_args()
    return args

def preprocess(h, acc, dataset, distribution, region, _merge_datasets=True):
    h = merge_extensions(h, acc, reweight_pu=False)
    scale_xs_lumi(h)
    if _merge_datasets:
        h = merge_datasets(h)

    h = h.integrate('region', region).integrate('dataset', dataset)

    if distribution == 'mjj':
        mjj_ax = hist.Bin('mjj', r'$M_{jj} \ (GeV)$', [200., 400., 600., 900., 1200., 1500., 2000., 2750., 3500., 5000.])
        h = h.rebin('mjj', mjj_ax)

    return h

def plot_rereco_to_ul_comparison(acc_dict, outtag, dataset='MET_2018', year=2018, dataset_tag='MET_2018', distribution='mjj', region='sr_vbf'):
    '''For the given dataset, plot the ReReco v7 vs. UL v8 comparison.'''
    histos = {}
    for key, acc in acc_dict.items():
        acc.load(distribution)
        histos[key] = preprocess(acc[distribution], acc, dataset, distribution, region, _merge_datasets='MET' in dataset_tag)

    fig, ax, rax = fig_ratio()
    hist.plot1d(histos['03Sep20v7'], ax=ax)
    hist.plot1d(histos['ULv8'], ax=ax, clear=False)

    ax.set_yscale('log')
    ax.set_ylim(AX_YLIMS[distribution])

    ax.legend(title='Version', labels=[
        'ReReco v7',
        'UL v8',
    ])

    ax.text(0.,1.,f'{dataset_tag.replace("_", " ")} {year}',
        fontsize=14,
        ha='left',
        va='bottom',
        transform=ax.transAxes
    )

    ax.text(0.05,0.9,r'$400 < H_T^{MG} < 1200 \ GeV$',
        fontsize=14,
        ha='left',
        va='bottom',
        transform=ax.transAxes
    )

    ax.text(1.,1.,'VBF SR (No Cleaning Cuts)',
        fontsize=14,
        ha='right',
        va='bottom',
        transform=ax.transAxes
    )

    ax.yaxis.set_ticks_position('both')

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
    rax.set_ylim(RAX_YLIMS[distribution])

    if RAX_YLIMS[distribution] == (0,2):
        loc = MultipleLocator(0.5)
        rax.yaxis.set_major_locator(loc)

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
    outpath = pjoin(outdir, f'{dataset_tag}_rereco_to_ul_comp_{region}_{distribution}_{year}.pdf')
    fig.savefig(outpath)
    plt.close(fig)

    print(f'File saved: {outpath}')

def main():
    args = parse_cli()
    outtag = args.outtag
    if not outtag:
        raise RuntimeError('Please provide an outtag for this job!')
    
    acc_dict = {
        '03Sep20v7' : dir_archive(bucoffea_path('submission/merged_2021-04-13_vbfhinv_03Sep20v7_QCD_W_HT_btw_400_1200_noCleaningCuts')),
        'ULv8' : dir_archive(bucoffea_path('submission/merged_2021-04-13_vbfhinv_ULv8_QCD_W_HT_btw_400_1200_noCleaningCuts')),
        # '03Sep20v7' : dir_archive(bucoffea_path('submission/merged_2021-04-12_vbfhinv_03Sep20v7_MET_2018_rerecoCleaningCuts')),
        # 'ULv8' : dir_archive(bucoffea_path('submission/merged_2021-04-12_vbfhinv_ULv8_MET_2018_rerecoCleaningCuts')),
    }
    for acc in acc_dict.values():
        acc.load('sumw')
        acc.load('sumw2')

    distributions = ['mjj', 'ak4_eta0', 'ak4_eta1']

    for year in args.years:
        dataset_regex = {
            'QCD_W' : re.compile(f'WJetsToLNu_HT.*{year}'),
            'Data' : re.compile(f'MET.*{year}'),
        }
        dataset = dataset_regex[args.proc]

        for distribution in distributions:
            if not re.match(args.distribution, distribution):
                continue
            plot_rereco_to_ul_comparison(acc_dict, outtag, dataset=dataset, year=year, dataset_tag=args.proc, distribution=distribution)

if __name__ == '__main__':
    main()