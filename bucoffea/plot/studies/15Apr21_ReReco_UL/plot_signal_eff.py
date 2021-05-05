#!/usr/bin/env python

import os
import re
import sys
import numpy as np

from coffea import hist
from matplotlib import pyplot as plt
from bucoffea.plot.util import merge_datasets, merge_extensions, scale_xs_lumi, fig_ratio
from klepto.archives import dir_archive

pjoin = os.path.join

def plot_signal_eff(acc, outtag, distribution='mjj'):
    '''Plot signal efficiency of the HF shape cuts in terms of the given variable.'''
    acc.load(distribution)
    h = acc[distribution]

    h = merge_extensions(h, acc)
    scale_xs_lumi(h)
    h = merge_datasets(h)

    if distribution == 'mjj':
        mjj_ax = hist.Bin('mjj', r'$M_{jj} \ (GeV)$', [200., 400., 600., 900., 1200., 1500.,2000., 2750., 3500., 5000.])
        h = h.rebin('mjj', mjj_ax)
    
    for year in [2017, 2018]:
        _h = h.integrate('dataset', re.compile(f'VBF_HToInv.*{year}'))
        fig, ax, rax = fig_ratio()
        hist.plot1d(_h, ax=ax, overlay='region')

        ax.set_yscale('log')
        ax.set_ylim(1e-3,1e5)

        ax.text(0.,1.,f'VBF H(inv) {year}',
            fontsize=14,
            ha='left',
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
            _h.integrate('region', 'sr_vbf'),
            _h.integrate('region', 'sr_vbf_no_hfcut'),
            ax=rax,
            unc='num',
            error_opts=data_err_opts
        )

        rax.grid(True)
        rax.set_ylim(0,1.1)
        rax.set_ylabel('Efficiency')
        rax.axhline(1, xmin=0, xmax=1, color='red')

        outdir = f'./output/{outtag}'
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        outpath = pjoin(outdir, f'signal_eff_{distribution}_{year}.pdf')
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

    plot_signal_eff(acc, outtag)

if __name__ == '__main__':
    main()