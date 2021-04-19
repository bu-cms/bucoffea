#!/usr/bin/env python

import os
import re
import argparse
import numpy as np

from coffea import hist
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator

from bucoffea.plot.util import merge_extensions, merge_datasets, scale_xs_lumi, fig_ratio
from bucoffea.plot.style import matplotlib_rc
from bucoffea.helpers.paths import bucoffea_path
from klepto.archives import dir_archive
from pprint import pprint
from distributions import binnings

matplotlib_rc()

pjoin = os.path.join
# Suppress true_divide warnings
np.seterr(divide='ignore', invalid='ignore')

def parse_cli():
    parser = argparse.ArgumentParser()
    parser.add_argument('--regex', default='.*', help='The regex specifying distributions to run.')
    parser.add_argument('--years', type=int, nargs='*', default=[2017,2018], help='Years to run the comparison for.')
    args = parser.parse_args()
    return args

def preprocess(h, acc, region, distribution, dataset_regex):
    h = merge_extensions(h, acc, reweight_pu=False)
    scale_xs_lumi(h)
    h = merge_datasets(h)

    if distribution in binnings.keys():
        new_ax = binnings[distribution]
        h = h.rebin(new_ax.name, new_ax)

    h = h.integrate('region', region).integrate('dataset', dataset_regex)

    return h

def compare_mc_to_mc(acc_dict, outtag, dataset_regex, dataset_tag, year, distribution='mjj', region='sr_vbf_no_veto_all', fformat='pdf'):
    '''MC to MC comparison between EOY and UL MC.'''
    histos = {}
    for key, acc in acc_dict.items():
        acc.load(distribution)
        histos[key] = preprocess(acc[distribution], acc, region, distribution, dataset_regex)

    if distribution  == 'mjj':
        overflow = 'over'
    else:
        overflow = 'none'

    fig, ax, rax = fig_ratio()
    hist.plot1d(histos['ReRecov7'], ax=ax, overflow=overflow)
    hist.plot1d(histos['ULv8'], ax=ax, clear=False, overflow=overflow)

    ax.legend(title='Version', labels=['ReRecov7', 'ULv8'])

    ax.set_yscale('log')
    ax.set_ylim(1e0,1e6)

    ax.text(0.,1.,f'{dataset_tag.replace("_"," ")}',
        fontsize=14,
        ha='left',
        va='bottom',
        transform=ax.transAxes
    )

    ax.text(1.,1.,f'VBF SR {year}',
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
        histos['ReRecov7'],
        ax=rax,
        unc='num',
        error_opts=data_err_opts,
        overflow=overflow
    )

    rax.grid(True)
    if 'ak4_eta' in distribution:
        rax.set_ylim(0,2)
        loc = MultipleLocator(0.5)
        rax.yaxis.set_major_locator(loc)
    else:
        rax.set_ylim(0.5,1.5)
        loc = MultipleLocator(0.25)
        rax.yaxis.set_major_locator(loc)
    rax.set_ylabel('UL / ReReco')

    new_xlabels = {
        'ak4_eta0' : r'Leading Jet $\eta$',
        'ak4_eta1' : r'Trailing Jet $\eta$',
    }

    if distribution in new_xlabels.keys():
        ax.set_xlabel(new_xlabels[distribution])
        rax.set_xlabel(new_xlabels[distribution])

    outdir = f'./output/{outtag}'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outpath = pjoin(outdir, f'{dataset_tag}_{distribution}_{year}.{fformat}')
    fig.savefig(outpath)
    plt.close(fig)

    print(f'File saved: {outpath}')

def main():
    args = parse_cli()
    acc_dict = {
        'ReRecov7' : dir_archive(bucoffea_path('submission/merged_2021-04-08_vbfhinv_03Sep20v7_MC_withHFReweight_ULv8_DATA_one_fifth_unblind_jetPt80')),
        'ULv8' : dir_archive(bucoffea_path('submission/merged_2021-04-18_vbfhinv_ULv8_05Feb21_withMuonRegions_v4')),
    }

    for acc in  acc_dict.values():
        acc.load('sumw')
        acc.load('sumw2')

    outtag = 'mc_to_mc_19Apr21'

    distributions = [
        'mjj',
        'met',
        'ak4_eta0',
        'ak4_eta1',
        'ak4_pt0',
        'ak4_pt1',
        'ak4_nef0',
        'ak4_nef1',
        'ak4_nhf0',
        'ak4_nhf1',
        'ak4_chf0',
        'ak4_chf1',
    ]

    for year in args.years:
        datasets = [
            (re.compile(f'ZJetsToNuNu.*{year}'), 'QCD_Zvv'),
            (re.compile(f'WJetsToLNu.*{year}'), 'QCD_Wlv'),
            (re.compile(f'EWKZ2Jets.*ZToNuNu.*{year}'), 'EWK_Zvv'),
            (re.compile(f'EWKW2Jets.*{year}'), 'EWK_Wlv'),
        ]
    
        for dataset_regex, dataset_tag in datasets:
            for distribution in distributions:
                if not re.match(args.regex, distribution):
                    continue
                for fformat in ['pdf', 'png']:
                    compare_mc_to_mc(acc_dict, outtag,
                        dataset_regex=dataset_regex,
                        dataset_tag=dataset_tag,
                        year=year,
                        distribution=distribution,
                        fformat=fformat,
                    )
    
if __name__ == '__main__':
    main()