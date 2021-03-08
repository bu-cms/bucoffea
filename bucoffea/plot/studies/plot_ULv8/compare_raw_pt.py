#!/usr/bin/env python

import os
import sys
import re
import numpy as np

from matplotlib import pyplot as plt
from coffea import hist
from bucoffea.plot.util import merge_datasets, merge_extensions, scale_xs_lumi, fig_ratio
from bucoffea.helpers.paths import bucoffea_path
from bucoffea.plot.style import matplotlib_rc
from klepto.archives import dir_archive
from pprint import pprint

pjoin = os.path.join

matplotlib_rc()

def preprocess(h, acc, region='sr_vbf'):
    h = merge_extensions(h, acc, reweight_pu=False)
    scale_xs_lumi(h)
    h = merge_datasets(h)

    h = h.integrate('region', region).integrate('dataset', 'MET_2017')

    new_ax = hist.Bin('jetpt', r'Jet $p_T$ (GeV)', 25, 0, 500)
    h = h.rebin('jetpt', new_ax)
    
    return h

def compare_raw_pt(acc_ULv7, acc_ULv8, outtag, distribution='ak4_ptraw0', region='sr_vbf'):
    acc_ULv7.load(distribution)
    acc_ULv8.load(distribution)

    h_v7 = preprocess(acc_ULv7[distribution], acc_ULv7, region)
    h_v8 = preprocess(acc_ULv8[distribution], acc_ULv8, region)

    fig, ax, rax = fig_ratio()
    hist.plot1d(h_v7, ax=ax)
    hist.plot1d(h_v8, ax=ax, clear=False)

    ax.set_yscale('log')
    ax.set_ylim(1e-1,1e7)

    ax.legend(title='UL version', labels=['02Dec2019 UL', 'Nanov8 UL'])

    data_err_opts = {
        'linestyle':'none',
        'marker': '.',
        'markersize': 10.,
        'color':'k',
    }

    hist.plotratio(h_v8, h_v7,
        ax=rax,
        unc='num',
        error_opts=data_err_opts
    )

    rax.set_ylabel('v8 / v7')
    rax.set_ylim(0.8,1.2)
    rax.grid(True)

    if distribution == 'ak4_ptraw0':
        xlabel = r'Leading jet raw $p_T$ (GeV)'
    elif distribution == 'ak4_ptraw1':
        xlabel = r'Trailing jet raw $p_T$ (GeV)'

    ax.set_xlabel(xlabel)
    rax.set_xlabel(xlabel)

    # Save figure
    outdir = f'./output/{outtag}'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outpath = pjoin(outdir, f'{distribution}.pdf')
    fig.savefig(outpath)
    plt.close(fig)

    print(f'File saved: {outpath}')

def main():
    acc_ULv7 = dir_archive( bucoffea_path('submission/merged_2021-03-02_vbfhinv_MET2017_oldUL') )
    acc_ULv8 = dir_archive( bucoffea_path('submission/merged_2021-03-02_vbfhinv_MET2017_newUL') )

    for acc in [acc_ULv7, acc_ULv8]:
        acc.load('sumw')
        acc.load('sumw2')

    outtag = '08Mar21_rawpt_comparison'

    for distribution in ['ak4_ptraw0', 'ak4_ptraw1']:
        compare_raw_pt(acc_ULv7, acc_ULv8, outtag, distribution)

if __name__ == '__main__':
    main()