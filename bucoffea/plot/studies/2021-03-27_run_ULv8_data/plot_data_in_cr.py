#!/usr/bin/env python

import os
from sre_compile import dis
import sys
import re
import numpy as np
import mplhep as hep

from matplotlib import pyplot as plt
from bucoffea.plot.util import merge_datasets, merge_extensions, scale_xs_lumi
from coffea import hist
from klepto.archives import dir_archive

pjoin = os.path.join

def plot_data_in_qcd_cr(acc, outtag, dataset='MET_2017', distribution='mjj', region='sr_vbf_fail_hf_cuts'):
    acc.load(distribution)
    h = acc[distribution]

    h = merge_extensions(h, acc, reweight_pu=False)
    scale_xs_lumi(h)
    h = merge_datasets(h)

    if distribution == 'mjj':
        mjj_ax = hist.Bin('mjj', r'$M_{jj} \ (GeV)$', [200., 400., 600., 900., 1200., 1500., 2000., 2750., 3500., 5000.])
        h = h.rebin('mjj', mjj_ax)

    h = h.integrate('region', region).integrate('dataset', dataset)

    fig, ax = plt.subplots()
    hist.plot1d(h, ax=ax)
    
    ax.set_yscale('log')
    ax.set_ylim(1e0,1e4)

    ax.text(0.,1.,'UL v8 MET 2017',
        fontsize=14,
        ha='left',
        va='bottom',
        transform=ax.transAxes
    )

    ax.text(1.,1.,'QCD CR',
        fontsize=14,
        ha='right',
        va='bottom',
        transform=ax.transAxes
    )

    ax.get_legend().remove()
    ax.yaxis.set_ticks_position('both')

    outdir = f'./output/{outtag}'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outpath = pjoin(outdir, f'data_in_qcd_cr_{distribution}.pdf')
    fig.savefig(outpath)
    plt.close(fig)

def main():
    inpath = sys.argv[1]
    acc = dir_archive(inpath)
    acc.load('sumw')
    acc.load('sumw2')

    outtag = re.findall('merged_.*', inpath)[0].replace('/','')

    plot_data_in_qcd_cr(acc, outtag)

if __name__ == '__main__':
    main()
