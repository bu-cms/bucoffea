#!/usr/bin/env python

import os
import sys
import re
import numpy as np

from matplotlib import pyplot as plt
from coffea import hist
from bucoffea.plot.util import fig_ratio
from bucoffea.helpers.paths import bucoffea_path
from bucoffea.plot.style import matplotlib_rc
from klepto.archives import dir_archive
from pprint import pprint

pjoin = os.path.join

matplotlib_rc()

def plot_ul_comparison(acc_old, acc_new, jobtag, distribution='ak4_eta0', region='sr_vbf', dataset_regex='MET.*2017C'):
    '''Plot the comparison of two UL samples: Old and new (most recent v8)'''
    acc_old.load(distribution)
    acc_new.load(distribution)

    # Take the data from the desired region and dataset
    h_old = acc_old[distribution].integrate('region', region).integrate('dataset', re.compile(dataset_regex))
    h_new = acc_new[distribution].integrate('region', region).integrate('dataset', re.compile(dataset_regex))

    # Rebin mjj
    if distribution == 'mjj':
        mjj_ax = hist.Bin('mjj', r'$M_{jj} \ (GeV)$', [200., 400., 600., 900., 1200., 1500., 2000., 2750., 3500., 5000.])
        h_old = h_old.rebin('mjj', mjj_ax)
        h_new = h_new.rebin('mjj', mjj_ax)

    fig, ax, rax = fig_ratio()
    hist.plot1d(h_old, ax=ax)
    hist.plot1d(h_new, ax=ax, clear=False)

    ax.set_yscale('log')
    ax.set_ylim(1e-1,1e5)

    ax.legend(labels=['Old UL', 'New v8 UL'])

    if re.match('.*2017[A-F]', dataset_regex):
        dataset_tag = dataset_regex.replace('.*', '_') 
        ax.text(1., 1., dataset_tag,
            fontsize=14,
            ha='right',
            va='bottom',
            transform=ax.transAxes
        )

    ax.text(0., 1., 'VBF Signal Region',
        fontsize=14,
        ha='left',
        va='bottom',
        transform=ax.transAxes
    )

    # Plot ratio
    data_err_opts = {
        'linestyle':'none',
        'marker': '.',
        'markersize': 10.,
        'color':'k',
    }
    
    hist.plotratio(
        h_new, h_old, 
        ax=rax, 
        unc='num',
        error_opts=data_err_opts
        )

    rax.set_ylim(0.5,1.5)
    rax.set_ylabel('New UL / Old UL')

    rax.grid(True)

    new_xlabels = {
        'ak4_eta0' : r'Leading Jet $\eta$',
        'ak4_eta1' : r'Trailing Jet $\eta$',
    }

    if distribution in new_xlabels.keys():
        rax.set_xlabel(new_xlabels[distribution])
        ax.set_xlabel(new_xlabels[distribution])

    # Save figure
    outdir = f'./output/{jobtag}'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    outpath = pjoin(outdir, f'ul_comparison_{region}_{distribution}.pdf')
    fig.savefig(outpath)
    plt.close(fig)
    
    print(f'File saved: {outpath}')

def main():
    acc_old_UL = dir_archive( bucoffea_path('./submission/merged_2021-02-23_vbfhinv_MET2017_RunC_oldUL') )
    acc_new_UL = dir_archive( bucoffea_path('./submission/merged_2021-02-18_vbfhinv_ULv8_MET_2017') )

    acc_old_UL.load('sumw')
    acc_old_UL.load('sumw2')
    acc_new_UL.load('sumw')
    acc_new_UL.load('sumw2')

    jobtag = '23Feb21'

    distributions = ['ak4_eta0', 'ak4_eta1', 'mjj', 'detajj']
    for distribution in distributions:
        plot_ul_comparison(acc_old_UL, acc_new_UL, jobtag, distribution=distribution)

if __name__ == '__main__':
    main()