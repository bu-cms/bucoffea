#!/usr/bin/env python

import os
import sys
import re
import numpy as np
import matplotlib.colors as colors

from coffea import hist
from bucoffea.plot.util import merge_datasets, merge_extensions, scale_xs_lumi
from matplotlib import pyplot as plt
from klepto.archives import dir_archive
from pprint import pprint

pjoin = os.path.join

def plot_leading_jet_eta_vs_pt(acc, outtag, region='cr_2m_vbf_relaxed_sel', dataset_tag='MET_2017'):
    distribution = 'ak4_pt0_eta0'
    acc.load(distribution)
    h = acc[distribution]

    h = merge_extensions(h, acc, reweight_pu=False)
    scale_xs_lumi(h)
    h = merge_datasets(h)

    datasetregex = {
        'MET_2017' : 'MET_2017',
        'DY_2017' : re.compile('DYJetsToLL.*2017'),
    }

    h = h.integrate('region', region).integrate('dataset', datasetregex[dataset_tag])

    pt_ax = hist.Bin('jetpt', r'Leading Jet $p_T \ (GeV)$', list(range(80,780,50)))
    h = h.rebin('jetpt', pt_ax)

    fig, ax = plt.subplots()
    hist.plot2d(h, ax=ax, xaxis='jeteta',patch_opts={'norm' : colors.LogNorm(1e-1,1e3)})
    ax.set_xlabel(r'Leading Jet $\eta$')

    ax.text(0.,1.,dataset_tag.replace('_', ' '),
        fontsize=14,
        ha='left',
        va='bottom',
        transform=ax.transAxes
    )

    ax.text(1.,1.,r'$Z(\mu\mu)$',
        fontsize=14,
        ha='right',
        va='bottom',
        transform=ax.transAxes
    )

    outdir = f'./output/{outtag}'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    outpath = pjoin(outdir, f'{dataset_tag}_{distribution}.pdf')
    fig.savefig(outpath)
    plt.close(fig)

    print(f'File saved: {outpath}')

def main():
    inpath = sys.argv[1]
    acc = dir_archive(inpath)
    acc.load('sumw')
    acc.load('sumw2')
    
    outtag = re.findall('merged_.*', inpath)[0].replace('/','')

    dataset_tag = 'DY_2017' if 'DY' in outtag else 'MET_2017'

    plot_leading_jet_eta_vs_pt(acc, outtag, dataset_tag=dataset_tag)

if __name__ == '__main__':
    main()