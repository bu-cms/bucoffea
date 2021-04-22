#!/usr/bin/env python
import os
import re
import sys
import numpy as np
import mplhep as hep

from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
from coffea import hist
from bucoffea.plot.util import merge_datasets, merge_extensions, scale_xs_lumi, fig_ratio
from bucoffea.helpers.paths import bucoffea_path
from klepto.archives import dir_archive
from pprint import pprint

pjoin = os.path.join

def preprocess(h, acc, distribution, region, dataset_regex):
    h = merge_extensions(h, acc, reweight_pu=False)
    scale_xs_lumi(h)
    h = merge_datasets(h)

    if distribution == 'mjj':
        new_ax = hist.Bin('mjj', r'$M_{jj} \ (GeV)$', [200., 400., 600., 900., 1200., 1500., 2000., 2750., 3500., 5000.])
        h = h.rebin('mjj', new_ax)

    h = h.integrate('region', region).integrate('dataset', dataset_regex)
    return h

def compare_muon_sf(acc_dict, outtag, dataset_tag, dataset_regex, distribution='mjj', year=2017, region='cr_2m_vbf'):
    '''Compare the distribution of the given dataset with EOY and UL muon SF.'''
    histos = {}
    for key, acc in acc_dict.items():
        acc.load(distribution)
        histos[key] = preprocess(acc[distribution], acc, distribution, region, dataset_regex)

    fig, ax, rax = fig_ratio()
    hist.plot1d(histos['ReReco'], ax=ax)
    hist.plot1d(histos['UL'], ax=ax, clear=False)

    ax.set_yscale('log')
    ax.set_ylim(1e-1,1e5)

    ax.text(0.,1.,dataset_tag,
        fontsize=14,
        ha='left',
        va='bottom',
        transform=ax.transAxes
    )

    ax.legend(title='Muon SF', labels=['ReReco', 'UL'])

    data_err_opts = {
        'linestyle':'none',
        'marker': '.',
        'markersize': 10.,
        'color':'k',
    }

    hist.plotratio(
        histos['UL'],
        histos['ReReco'],
        ax=rax,
        unc='num',
        error_opts=data_err_opts
    )

    rax.grid(True)
    rax.set_ylim(0.95,1.05)
    rax.set_ylabel('UL / ReReco')

    outdir = f'./output/{outtag}'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outpath = pjoin(outdir, f'{dataset_tag}_muon_sf_{distribution}_{year}.pdf')
    fig.savefig(outpath)
    plt.close(fig)

    print(f'File saved: {outpath}')

def main():
    # Accumulators with EOY and UL muon SF
    acc_dict = {
        'ReReco' : dir_archive( bucoffea_path('submission/merged_2021-04-21_vbfhinv_ULv8_05Feb21_one_fifth_unblind_prefireUL') ),
        'UL' : dir_archive( bucoffea_path('submission/merged_2021-04-22_vbfhinv_ULv8_05Feb21_Zmm_muonSF_UL') ),
    }

    for acc in acc_dict.values():
        acc.load('sumw')
        acc.load('sumw2')

    outtag = '22Apr21_EOY_UL_muonSF'

    for year in [2017, 2018]:
        dataset_info = [
            ('qcd_zmm', re.compile(f'DYJetsToLL.*{year}')),
            ('ewk_zmm', re.compile(f'EWKZ2Jets.*ZToLL.*{year}')),
        ]

        for dataset_tag, dataset_regex in dataset_info:
            compare_muon_sf(acc_dict, outtag,
                dataset_tag=dataset_tag,
                dataset_regex=dataset_regex,
                year=year,
            )

if __name__ == '__main__':
    main()