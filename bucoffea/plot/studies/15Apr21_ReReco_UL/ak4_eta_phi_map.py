#!/usr/bin/env python

import os
import sys
import re
import numpy as np

from matplotlib import pyplot as plt
from coffea import hist
from bucoffea.plot.util import merge_datasets, merge_extensions, scale_xs_lumi
from klepto.archives import dir_archive
from pprint import pprint

pjoin = os.path.join

def plot_eta_phi_map(acc, outtag, distribution='ak4_eta_phi', region='sr_vbf'):
    '''Plot eta/phi map for data in SR'''
    acc.load(distribution)
    h = acc[distribution]

    h = merge_extensions(h, acc)
    scale_xs_lumi(h)
    h = merge_datasets(h)

    h = h.integrate('region',  region)
    outdir = f'./output/{outtag}'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    for year in [2017, 2018]:
        fig, ax = plt.subplots()
        hist.plot2d(
            h[f'MET_{year}'].integrate('dataset'), 
            ax=ax,
            xaxis='jeteta',
        )

        ax.text(0,1,f'MET {year}',
            fontsize=14,
            ha='left',
            va='bottom',
            transform=ax.transAxes
        )
        
        outpath = pjoin(outdir, f'ak4_eta_phi_{year}.pdf')
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

    plot_eta_phi_map(acc, outtag)

if __name__ == '__main__':
    main()