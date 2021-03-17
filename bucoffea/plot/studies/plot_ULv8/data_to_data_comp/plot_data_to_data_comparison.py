#!/usr/bin/env python

import os
import sys
import re
import numpy as np

from coffea import hist
from bucoffea.plot.util import merge_datasets, merge_extensions, scale_xs_lumi, fig_ratio
from matplotlib import pyplot as plt
from klepto.archives import dir_archive
from pprint import pprint

pjoin = os.path.join

def plot_data_to_data_comparison(acc, outtag, distribution='mjj', datatag='data', regionbase='sr_vbf'):
    '''Plot data to data comparison before/after HF shape cuts are applied.'''
    acc.load(distribution)
    h = acc[distribution]

    h = merge_extensions(h, acc, reweight_pu=False)
    scale_xs_lumi(h)
    h = merge_datasets(h)

    if datatag == 'data':
        h = h.integrate('dataset', 'MET_2017')
    else:
        h = h.integrate('dataset', re.compile('DYJets.*2017'))

    if distribution == 'mjj':
        mjj_ax = hist.Bin('mjj', r'$M_{jj} \ (GeV)$', [200., 400., 600., 900., 1200., 1500., 2000., 2750., 3500., 5000.])
        h = h.rebin('mjj', mjj_ax)
            
    # Get the two regions: With and without the HF shape cuts applied
    h = h[re.compile(f'^{regionbase}$|^{regionbase}_with_hfcuts$')]

    fig, ax, rax = fig_ratio()
    hist.plot1d(h, ax=ax, overlay='region')

    ax.set_yscale('log')
    if datatag == 'data':
        ax.set_ylim(1e0,1e6)
    else:
        ax.set_ylim(1e-1,1e5)

    if 'nocleaningcuts' in regionbase:
        plottag = 'Usual Cleaning Cuts Removed'
    else:
        plottag = 'With Usual Cleaning Cuts'

    ax.text(0., 1., plottag,
        fontsize=12,
        ha='left',
        va='bottom',
        transform=ax.transAxes
    )

    if datatag == 'data':
        _text = 'MET 2017 UL Nanov8'
    else:
        _text = r'$Z(\mu\mu)$ 2017 UL Nanov8'

    ax.text(0.02, 0.02, _text,
        fontsize=12,
        ha='left',
        va='bottom',
        transform=ax.transAxes
    )
    
    ax.text(1., 1., r'$\sigma_{\eta\eta} - \sigma_{\phi\phi} < 0.03$ & $CSSize < 3$',
        fontsize=12,
        ha='right',
        va='bottom',
        transform=ax.transAxes
    )

    handles, labels = ax.get_legend_handles_labels()
    for handle, label in zip(handles, labels):
        if '_with_hfcuts' in label:
            handle.set_label('Applied')
        else:
            handle.set_label('Not applied')

    ax.legend(title='HF Shape Cuts', handles=handles)

    # Plot the ratio of the two
    data_err_opts = {
        'linestyle':'none',
        'marker': '.',
        'markersize': 10.,
        'color':'k',
    }

    hist.plotratio(
        h.integrate('region', f'{regionbase}_with_hfcuts'),
        h.integrate('region', f'{regionbase}'),
        ax=rax,
        unc='num',
        error_opts=data_err_opts
    )

    rax.grid(True)
    rax.set_ylabel('With HF cut / Without')
    if 'nocleaningcuts' in regionbase:
        rax.set_ylim(0,2)
    else:
        rax.set_ylim(0.5,1.5)

    # Save figure
    outdir = f'./output/{outtag}'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outpath = pjoin(outdir, f'data_to_data_comp_{regionbase}.pdf')
    fig.savefig(outpath)
    plt.close(fig)

    print(f'File saved: {outpath}')

def main():
    inpath = sys.argv[1]
    acc = dir_archive(inpath)
    acc.load('sumw')
    acc.load('sumw2')

    outtag = re.findall('merged_.*', inpath)[0].replace('/','')

    # Two regions: Regular SR + SR without the usual cleaning cuts
    regionbases = ['sr_vbf', 'sr_vbf_nocleaningcuts']

    datatag = 'mc' if 'DY' in inpath else 'data'

    regionbases_for_dtype = {
        'data' : ['sr_vbf', 'sr_vbf_nocleaningcuts'],
        'mc' : ['cr_2m_vbf_relaxed_sel', 'cr_2m_vbf_relaxed_sel_nocleaningcuts'],
    }

    regionbases = regionbases_for_dtype[datatag]

    for regionbase in regionbases:
        plot_data_to_data_comparison(acc, outtag, regionbase=regionbase, datatag=datatag)

if __name__ == '__main__':
    main()