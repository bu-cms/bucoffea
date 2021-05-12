#!/usr/bin/env python
import os
import re
import sys
import warnings
import numpy as np

from matplotlib import pyplot as plt
from coffea import hist
from scipy.stats import distributions
from bucoffea.plot.util import merge_datasets, merge_extensions, scale_xs_lumi, fig_ratio
from bucoffea.helpers.paths import bucoffea_path
from klepto.archives import dir_archive
from pprint import pprint

pjoin = os.path.join

warnings.filterwarnings('ignore')

BINNINGS = {
    'mjj' : hist.Bin('mjj', r'$M_{jj} \ (GeV)$', [200., 400., 600., 900., 1200., 1500., 2000., 2750., 3500.]),
}

def preprocess(h, acc, region, distribution='mjj'):
    h = merge_extensions(h, acc)
    scale_xs_lumi(h)
    # h = merge_datasets(h)
    
    # Rebin if necessary
    if distribution in BINNINGS.keys():
        new_ax = BINNINGS[distribution]
        h = h.rebin(new_ax.name, new_ax)

    return h.integrate('region', region)

def compare_per_dataset(acc_dict, outtag, region='sr_vbf_no_veto_all', distribution='mjj'):
    histos = {}
    for key, acc in acc_dict.items():
        acc.load(distribution)
        histos[key] = preprocess(acc[distribution], acc, region, distribution)

    outdir = f'./output/{outtag}'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    for year in [2017, 2018]:
        h_eoy = histos['EOY'][re.compile(f'(ST|TTJets).*{year}')]
        h_ul = histos['UL'][re.compile(f'(ST|TTJets).*{year}')]

        dataset_identifiers_eoy = list(map(str, h_eoy.identifiers('dataset')))
        dataset_identifiers_ul = list(map(str, h_ul.identifiers('dataset')))

        # Dump the dataset identifiers for both EOY and UL to compare
        outpath = pjoin(outdir, f'top_datasets_{year}.txt')
        with open(outpath, 'w+') as f:
            f.write(f'Top Samples {year}\n')
            f.write('\nEOY\n')

            f.write('\n'.join(dataset_identifiers_eoy))
            f.write('\nUL\n')
            f.write('\n'.join(dataset_identifiers_ul))

        # One-to-one comparison: Remove the ones that DO NOT exist in UL
        new_dataset_identifiers_eoy = []
        for d in dataset_identifiers_eoy:
            if 'hadronicDecays' in d or 'pow-pythia' in d or 'MLM' in d:
                continue
            new_dataset_identifiers_eoy.append(d)

        # Compare one by one
        for d_eoy, d_ul in zip(new_dataset_identifiers_eoy, dataset_identifiers_ul):
            _h_eoy = h_eoy[d_eoy]
            _h_ul = h_ul[d_ul]

            fig, ax, rax = fig_ratio()
            hist.plot1d(_h_eoy, ax=ax, overlay='dataset')
            hist.plot1d(_h_ul, ax=ax, overlay='dataset', clear=False)

            data_err_opts = {
                'linestyle':'none',
                'marker': '.',
                'markersize': 10.,
                'color':'k',
            }
            
            hist.plotratio(
                _h_eoy.integrate('dataset', d_eoy),
                _h_ul.integrate('dataset', d_ul),
                ax=rax,
                unc='num',
                error_opts=data_err_opts
            )
            
            rax.grid(True)
            rax.set_ylim(0,2)
            rax.set_ylabel('EOY / UL')
            
            outpath = pjoin(outdir, f'{d_eoy}.pdf')
            fig.savefig(outpath)
            plt.close(fig)

            print(f'File saved: {outpath}')

def main():
    acc_dict = {
        'EOY' : dir_archive( bucoffea_path('submission/merged_2021-03-19_vbfhinv_03Sep20v7_leadak4_notinHF') ),
        'UL'  : dir_archive( bucoffea_path('submission/merged_2021-05-03_vbfhinv_ULv8_05Feb21_trigSF_update') ),
    }

    for acc in acc_dict.values():
        acc.load('sumw')
        acc.load('sumw_pileup')
        acc.load('nevents')

    outtag = '11May21_top'

    compare_per_dataset(acc_dict, outtag)

if __name__ == '__main__':
    main()