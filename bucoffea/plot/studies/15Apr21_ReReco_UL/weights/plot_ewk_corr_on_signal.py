#!/usr/bin/env python

import os
import sys
import re
import numpy as np

from bucoffea.plot.util import merge_datasets, merge_extensions, scale_xs_lumi, fig_ratio
from coffea import hist
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
from klepto.archives import dir_archive
from pprint import pprint

pjoin = os.path.join

def plot_ewk_corr_on_signal(acc, outtag, year, dataset='VBF_HToInvisible.*', distribution='weights', region='sr_vbf', weight_type='ewk_corr_signal'):
    acc.load(distribution)
    h = acc[distribution]
    
    h = merge_extensions(h, acc)
    scale_xs_lumi(h)
    h = merge_datasets(h)

    h = h.integrate('region', region).integrate('dataset', re.compile(f'{dataset}.*{year}')).integrate('weight_type', weight_type)

    fig, ax = plt.subplots()
    hist.plot1d(h, ax=ax)
    
    ax.legend(labels=[r'1 + $\alpha_{EWK}^{VBF}$'])

    outdir = f'./output/{outtag}/ewk_correction'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outpath = pjoin(outdir, f'vbf_weights_{year}.pdf')
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

    for year in [2017, 2018]:
        plot_ewk_corr_on_signal(acc, outtag, year=year)

if __name__ == '__main__':
    main()