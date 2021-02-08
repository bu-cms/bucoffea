#!/usr/bin/env python

import os
import sys
import re
import numpy as np

from matplotlib import pyplot as plt
from bucoffea.plot.util import merge_datasets, merge_extensions, lumi
from coffea import hist
from klepto.archives import dir_archive

pjoin = os.path.join

def get_region_tag(region):
    mapping = {
        'sr_vbf' : 'VBF SR'
    }
    return mapping[region]

def plot_hf_distributions(acc, outtag, distribution, year, region='sr_vbf'):
    acc.load(distribution)
    h = acc[distribution]

    h = merge_extensions(h, acc, reweight_pu=False)
    h = merge_datasets(h)

    h = h.integrate('dataset', f'MET_{year}').integrate('region', region)

    fig, ax = plt.subplots()
    hist.plot1d(h, ax=ax)

    ax.set_yscale('log')
    ax.set_ylim(1e-2,1e6)

    ax.get_legend().remove()

    ax.text(0., 1., f'MET {year}',
        fontsize=14,
        ha='left',
        va='bottom',
        transform=ax.transAxes
    )

    ax.text(1., 1., f'{get_region_tag(region)}, {lumi(year):.1f} fb$^{{-1}}$',
        fontsize=14,
        ha='right',
        va='bottom',
        transform=ax.transAxes
    )

    outdir = f'./output/{outtag}'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    outpath = pjoin(outdir, f'met_{year}_{distribution}.pdf')
    fig.savefig(outpath)
    plt.close(fig)
    print(f'File saved: {outpath}')

def main():
    inpath = sys.argv[1]
    acc = dir_archive(inpath)
    acc.load('sumw')
    acc.load('sumw2')

    outtag = re.findall('merged_.*', inpath)[0].replace('/', '')

    distributions = [
        'ak4_sigma_eta_eta0',
        'ak4_sigma_eta_eta1',
        'ak4_sigma_phi_phi0',
        'ak4_sigma_phi_phi1',
    ]

    for distribution in distributions:
        plot_hf_distributions(acc, outtag, distribution, year=2017)

if __name__ == '__main__':
    main()