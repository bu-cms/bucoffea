#!/usr/bin/env python

import os
import sys
import re
import numpy as np
import mplhep as hep

from matplotlib import pyplot as plt
from bucoffea.plot.util import merge_datasets, merge_extensions, scale_xs_lumi
from coffea import hist
from klepto.archives import dir_archive

pjoin = os.path.join

def plot_mc_in_sr(acc, outtag, region='sr_vbf', distribution='mjj'):
    '''Plot total MC background in SR'''
    acc.load(distribution)
    h = acc[distribution]

    h = merge_extensions(h, acc, reweight_pu=False)
    scale_xs_lumi(h)
    h = merge_datasets(h)

    if distribution == 'mjj':
        mjj_ax = hist.Bin('mjj', r'$M_{jj} \ (GeV)$', [200., 400., 600., 900., 1200., 1500., 2000., 2750., 3500., 5000.])
        h = h.rebin('mjj', mjj_ax)

    dataset = re.compile('(ZJetsToNuNu.*|EW.*|Top_FXFX.*|Diboson.*|.*DYJetsToLL_M-50_HT_MLM.*|.*WJetsToLNu.*HT.*).*2017')
    h = h.integrate('region', region).integrate('dataset', dataset)

    fig, ax = plt.subplots()
    hist.plot1d(h, ax=ax)

    ax.set_yscale('log')
    ax.set_ylim(1e0,1e6)

    ax.text(0.,1.,'Total MC Bkg',
        fontsize=14,
        ha='left',
        va='bottom',
        transform=ax.transAxes
    )

    ax.text(1.,1.,'VBF SR',
        fontsize=14,
        ha='right',
        va='bottom',
        transform=ax.transAxes
    )

    ax.get_legend().remove()

    outdir = f'./output/{outtag}'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outpath = pjoin(outdir, f'mc_in_sr_{distribution}.pdf')
    fig.savefig(outpath)
    plt.close(fig)

def main():
    inpath = sys.argv[1]
    acc = dir_archive(inpath)
    acc.load('sumw')
    acc.load('sumw2')

    outtag = re.findall('merged_.*', inpath)[0].replace('/','')

    plot_mc_in_sr(acc, outtag)

if __name__ == '__main__':
    main()
