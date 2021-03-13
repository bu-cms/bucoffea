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

def plot_met_for_hf(acc, outtag, distribution, dataset='MET_2017', region='cr_2m_vbf_relaxed_sel'):
    acc.load(distribution)
    h = acc[distribution]

    h = merge_extensions(h, acc, reweight_pu=False)
    scale_xs_lumi(h)
    h = merge_datasets(h)

    h = h.integrate('dataset', dataset).integrate('region', region)

    new_met_ax = hist.Bin('met', r'$p_T^{miss}$ (GeV)', 50, 0, 500)
    h = h.rebin('met', new_met_ax)

    fig, ax = plt.subplots()
    hist.plot1d(h, ax=ax)

    ax.set_yscale('log')
    ax.set_ylim(1e-1, 1e5)

    ax.axvline(50, ymin=0, ymax=1, color='k')

    if distribution == 'met_pt_ak40_hf':
        _text = 'Leading jet in HF'
    elif distribution == 'met_pt_ak41_hf':
        _text = 'Trailing jet in HF'
    else:
        raise RuntimeError(f'Check distribution: {distribution}')

    ax.text(0., 1., _text,
        fontsize=14,
        ha='left',
        va='bottom',
        transform=ax.transAxes
    )
    
    ax.text(1., 1., r'$Z(\mu\mu)$',
        fontsize=14,
        ha='right',
        va='bottom',
        transform=ax.transAxes
    )

    outdir = f'./output/{outtag}/met_check'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outpath = pjoin(outdir, f'{distribution}.pdf')

    fig.savefig(outpath)
    plt.close(fig)
    print(f'File saved: {outpath}')

def main():
    inpath = sys.argv[1]
    acc = dir_archive(inpath)
    acc.load('sumw')
    acc.load('sumw2')

    outtag = re.findall('merged_.*', inpath)[0].replace('/','')

    distributions = ['met_pt_ak40_hf', 'met_pt_ak41_hf']
    for distribution in distributions:
        plot_met_for_hf(acc, outtag, distribution=distribution)

if __name__ == '__main__':
    main()