#!/usr/bin/env python

import os
import sys
import re
import numpy as np

from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
from coffea import hist
from coffea.hist import poisson_interval
from bucoffea.plot.util import merge_datasets, merge_extensions, scale_xs_lumi, fig_ratio
from klepto.archives import dir_archive
from pprint import pprint

pjoin = os.path.join

def prefire_comparison(acc, outtag, mc, distribution='mjj'):
    '''Plot comparison of given MC with and without prefiring weights applied.'''
    acc.load(distribution)
    h = acc[distribution]

    h = merge_extensions(h, acc)
    scale_xs_lumi(h)
    h = merge_datasets(h)

    if distribution == 'mjj':
        mjj_bins = [200., 400., 600., 900., 1200., 1500., 2000., 2750., 3500., 5000.]
        mjj_ax = hist.Bin('mjj', r'$M_{jj} \ (GeV)$', mjj_bins)
        h = h.rebin('mjj', mjj_ax)

    h = h.integrate('dataset', mc)[re.compile('sr_vbf(_no_pref)?')]

    fig, ax, rax = fig_ratio()
    hist.plot1d(h, ax=ax, overlay='region', binwnorm=1)

    ax.set_yscale('log')
    ax.set_ylim(1e-4,1e4)
    ax.set_ylabel('Events / GeV')
    
    data_err_opts = {
        'linestyle':'none',
        'marker': '.',
        'markersize': 10.,
        'color':'k',
    }

    signal_err_opts = {
        'color':'crimson',
        'elinewidth': 1,
    }

    hist.plotratio(
        h.integrate('region', 'sr_vbf'),
        h.integrate('region', 'sr_vbf_no_pref'),
        ax=rax,
        unc='num',
        error_opts=data_err_opts
    )

    rax.grid(True)
    rax.set_ylim(0.8,1.2)
    rax.set_ylabel('With / Without Pref.')

    outdir = f'./output/{outtag}'
    try:
        os.makedirs(outdir)
    except FileExistsError:
        pass

    outpath = pjoin(outdir, f'prefire_comp.pdf')
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

    mc = re.compile(f'(ZJetsToNuNu.*|EW.*|Top_FXFX.*|Diboson.*|.*DYJetsToLL_M-50_HT_MLM.*|.*WJetsToLNu.*HT.*).*2017')

    prefire_comparison(acc, outtag, mc=mc)

if __name__ == '__main__':
    main()