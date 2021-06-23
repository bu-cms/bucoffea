#!/usr/bin/env python

import os
import sys
import re
import numpy as np

from matplotlib import pyplot as plt
from coffea import hist
from bucoffea.helpers.paths import bucoffea_path
from bucoffea.plot.util import merge_datasets, merge_extensions, scale_xs_lumi, fig_ratio
from klepto.archives import dir_archive
from pprint import pprint

pjoin = os.path.join

pretty_dataset_label = {
    'WJetsToLNu' : r'$W(e\nu)$ {year}',
    'DYJetsToLL' : r'$Z(ee)$ {year}',
}

def ele_trig_sf_uncs(acc, outtag, dataset, region, distribution='mjj_trig_weight'):
    acc.load(distribution)
    h = acc[distribution]

    h = merge_extensions(h, acc)
    scale_xs_lumi(h)
    h = merge_datasets(h)

    if 'mjj' in distribution:
        mjj_bins = [200., 400., 600., 900., 1200., 1500., 2000., 2750., 3500., 5000.]
        mjj_ax = hist.Bin('mjj', r'$M_{jj} \ (GeV)$', mjj_bins)
        h = h.rebin('mjj', mjj_ax)            

    h = h.integrate('region', region)

    outdir = f'./output/{outtag}'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    data_err_opts = {
        'linestyle':'none',
        'marker': '.',
        'markersize': 10.,
    }

    for year in [2017, 2018]:
        _h = h.integrate('dataset', re.compile(f'{dataset}.*{year}'))
        fig, ax, rax = fig_ratio()
        hist.plot1d(_h, ax=ax, overlay='variation')

        ax.text(0., 1., pretty_dataset_label[dataset].format(year=year),
            fontsize=14,
            ha='left',
            va='bottom',
            transform=ax.transAxes
        )

        ax.set_yscale('log')
        ax.set_ylim(1e-2,1e6)

        hist.plotratio(
            _h.integrate('variation', 'trigger_ele_up'),
            _h.integrate('variation', 'nominal'),
            ax=rax,
            unc='num',
            label='Up',
            error_opts=data_err_opts
        )

        hist.plotratio(
            _h.integrate('variation', 'trigger_ele_down'),
            _h.integrate('variation', 'nominal'),
            ax=rax,
            unc='num',
            label='Down',
            error_opts=data_err_opts,
            clear=False
        )

        rax.legend(ncol=2)
        rax.grid(True)
        rax.set_ylim(0.95,1.05)
        rax.set_ylabel('Variation')

        outpath = pjoin(outdir, f'{dataset}_ele_trigsf_uncs_{year}.pdf')        
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

    datasets_regions = [
        ('WJetsToLNu', 'cr_1e_vbf'),
        ('DYJetsToLL', 'cr_2e_vbf'),
    ]

    for dataset, region in datasets_regions:
        ele_trig_sf_uncs(acc, outtag,
            dataset=dataset,
            region=region
        )

if __name__ == '__main__':
    main()