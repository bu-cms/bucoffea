#!/usr/bin/env python

import os
import sys
import re
import numpy as np
import matplotlib.colors as colors

from bucoffea.plot.util import merge_datasets, merge_extensions, scale_xs_lumi, fig_ratio
from coffea import hist
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
from klepto.archives import dir_archive
from pprint import pprint

pjoin = os.path.join

REBIN = {
    'met' : hist.Bin('met',r'GEN $p_{T}^{miss}$ (GeV)',list(range(0,500,50)) + list(range(500,1000,100)) + list(range(1000,2000,250))),
    'mjj' : hist.Bin('mjj',r'$M_{jj}$ (GeV)', list(range(0,3000,100)) + list(range(3000,5500,500)))
}

def plot_gen_met_vs_mjj(acc, outtag, distribution='gen_met_mjj'):
    '''Plot GenMET vs. mjj for signal'''
    acc.load(distribution)
    h = acc[distribution]

    h = merge_extensions(h, acc)
    scale_xs_lumi(h)
    h = merge_datasets(h)

    for axis in ['met', 'mjj']:
        h = h.rebin(axis, REBIN[axis])

    h = h.integrate('region', 'sr_vbf_no_veto_all')

    for year in [2017, 2018]:
        _h = h.integrate('dataset', re.compile(f'VBF_HToInvisible.*{year}'))
        fig, ax = plt.subplots()
        hist.plot2d(_h, ax=ax, xaxis='mjj', patch_opts={'norm' : colors.LogNorm(1e-2,1e2)})

        ax.text(0.,1.,f'VBF H(inv) {year}',
            fontsize=14,
            ha='left',
            va='bottom',
            transform=ax.transAxes
        )

        outdir = f'./output/{outtag}'
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        outpath = pjoin(outdir, f'gen_met_mjj_vbf_{year}.pdf')
        fig.savefig(outpath)
        plt.close(fig)

        print(f'File saved: {outpath}')

def main():
    inpath = sys.argv[1]
    acc = dir_archive(inpath)
    acc.load('sumw')
    acc.load('sumw_pileup')
    acc.load('nevents')

    outtag = re.findall('merged_.*', inpath)[0].replace('/', '')

    plot_gen_met_vs_mjj(acc, outtag)

if __name__ == '__main__':
    main()