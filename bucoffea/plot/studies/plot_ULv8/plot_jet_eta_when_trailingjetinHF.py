#!/usr/bin/env python

import os
import sys
import re
import numpy as np

from coffea import hist
from bucoffea.plot.util import merge_datasets, merge_extensions, scale_xs_lumi
from matplotlib import pyplot as plt
from klepto.archives import dir_archive
from pprint import pprint

pjoin = os.path.join

def plot_jet_eta(acc, outtag, distribution, dataset='MET_2017', region='cr_vbf_qcd'):
    acc.load(distribution)
    h = acc[distribution]

    h = merge_extensions(h, acc, reweight_pu=False)
    scale_xs_lumi(h)
    h = merge_datasets(h)

    h = h.integrate('dataset', dataset).integrate('region', region)

    fig, ax = plt.subplots()
    hist.plot1d(h, ax=ax)
    
    ax.set_yscale('log')
    ax.set_ylim(1e-1,1e5)

    xlabels = {
        'ak4_eta0_trailjetHF' : r'Leading Jet $\eta$',
        'ak4_eta1_trailjetHF' : r'Trailing Jet $\eta$',
    }
    ax.set_xlabel(xlabels[distribution])

    ax.text(0., 1., r'Trailing Jet: $|\eta|>2.9$',
        fontsize=14,
        ha='left',
        va='bottom',
        transform=ax.transAxes
    )

    outdir = f'./output/{outtag}'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    outpath = pjoin(outdir, f'MET_2017_{region}_{distribution}.pdf')
    fig.savefig(outpath)
    plt.close(fig)
    print(f'File saved: {outpath}')

def main():
    inpath = sys.argv[1]
    acc = dir_archive(inpath)
    acc.load('sumw')
    acc.load('sumw2')

    outtag = re.findall('merged_.*', inpath)[0].replace('/','')

    distributions = [
        'ak4_eta0_trailjetHF',
        'ak4_eta1_trailjetHF',
    ]

    for distribution in distributions:
        plot_jet_eta(acc, outtag, distribution=distribution)

if __name__ == '__main__':
    main()