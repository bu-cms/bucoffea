#!/usr/bin/env python
import os
import sys
import re
import numpy as np

from matplotlib import pyplot as plt
from coffea import hist
from bucoffea.plot.util import merge_extensions, merge_datasets, scale_xs_lumi, fig_ratio
from klepto.archives import dir_archive
from pprint import pprint

pjoin = os.path.join

def preprocess(h, acc, region, dataset):
    h = merge_extensions(h,acc)
    scale_xs_lumi(h)
    h = merge_datasets(h)

    mjj_bins = [200., 400., 600., 900., 1200., 1500., 2000., 2750., 3500., 5000.]
    mjj_ax = hist.Bin('mjj', r'$M_{jj} \ (GeV)$', mjj_bins)
    h = h.rebin('mjj', mjj_ax)

    h = h.integrate('region', region).integrate('dataset', re.compile(f'{dataset}.*2017'))

    return h

def get_process_label(region):
    mapping = {
        'cr_2e_vbf' : r'$Z(ee)$',
        'cr_2m_vbf' : r'$Z(\mu\mu)$',
    }
    return mapping[region]

def plot_mjj_pref_nopref(acc, outtag, region, dataset):
    '''Plot mjj distribution with and without prefire weights.'''
    acc.load('mjj')
    acc.load('mjj_nopref')

    h_withpref = preprocess(acc['mjj'], acc, region=region, dataset=dataset)
    h_nopref = preprocess(acc['mjj_nopref'], acc, region=region, dataset=dataset)

    fig, ax, rax = fig_ratio()
    data_err_opts = {
        'linestyle':'none',
        'marker': '.',
        'markersize': 10.,
        'color':'k',
        'elinewidth': 1,
    }
    hist.plot1d(h_withpref, ax=ax, error_opts=data_err_opts)
    hist.plot1d(h_nopref, ax=ax, clear=False)

    ax.legend(labels=['Applied', 'Not Applied'], title='Prefire Weight')

    ax.set_yscale('log')
    ax.set_ylim(1e-4,1e4)
    ax.set_ylabel('Events')

    ax.text(0,1,get_process_label(region),
        fontsize=14,
        ha='left',
        va='bottom',
        transform=ax.transAxes
    )

    hist.plotratio(
        h_withpref,
        h_nopref,
        ax=rax,
        unc='num',
        error_opts=data_err_opts
    )

    rax.grid(True)
    rax.set_ylim(0.5,1.5)
    rax.set_ylabel('With pref / No pref')

    outdir = f'./output/{outtag}/prefire'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outpath = pjoin(outdir, f'{dataset}_{region}.pdf')
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

    regions = ['cr_2e_vbf', 'cr_2m_vbf']
    for region in regions:
        plot_mjj_pref_nopref(acc, outtag, region=region, dataset='DYJetsToLL')

if __name__ == '__main__':
    main()