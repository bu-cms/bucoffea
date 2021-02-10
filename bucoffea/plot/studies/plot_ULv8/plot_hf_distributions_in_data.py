#!/usr/bin/env python

import os
import sys
import re
import numpy as np
import matplotlib.colors as colors

from matplotlib import pyplot as plt
from bucoffea.plot.util import merge_datasets, merge_extensions, lumi
from coffea import hist
from klepto.archives import dir_archive

pjoin = os.path.join

def get_new_xlabel(distribution):
    mapping = {
        'ak4_sigma_eta_eta0' : r'Leading Jet $\sigma_{\eta\eta}$',
        'ak4_sigma_eta_eta1' : r'Trailing Jet $\sigma_{\eta\eta}$',
        'ak4_sigma_phi_phi0' : r'Leading Jet $\sigma_{\phi\phi}$',
        'ak4_sigma_phi_phi1' : r'Trailing Jet $\sigma_{\phi\phi}$',
    }

    return mapping[distribution]

def get_new_legend_label(label):
    mapping = {
        'sr_vbf' : 'VBF SR',
        'sr_vbf_hfhf' : 'VBF SR (HF-HF)',
        'sr_vbf_at_least_one_jet_in_hf' : 'VBF SR ($\\geq 1$ jet in HF)',
    }
    return mapping[label]

def plot_hf_distributions(acc, outtag, distribution, year=2017, region_regex='sr_vbf.*'):
    acc.load(distribution)
    h = acc[distribution]

    h = merge_extensions(h, acc, reweight_pu=False)
    h = merge_datasets(h)

    # h = h.integrate('dataset', f'MET_{year}').integrate('region', region)
    h = h.integrate('dataset', f'MET_{year}')[re.compile(region_regex)]

    fig, ax = plt.subplots()
    hist.plot1d(h, ax=ax, overlay='region')

    ax.set_yscale('log')
    ax.set_ylim(1e-2,1e6)
    ax.set_xlabel(get_new_xlabel(distribution))

    ax.text(0., 1., f'MET UL v8 {year}',
        fontsize=14,
        ha='left',
        va='bottom',
        transform=ax.transAxes
    )

    ax.text(1., 1., f'VBF SR',
        fontsize=14,
        ha='right',
        va='bottom',
        transform=ax.transAxes
    )

    # Update legend labels
    handles, labels = ax.get_legend_handles_labels()
    for handle, label in zip(handles, labels):
        new_label = get_new_legend_label(label)
        handle.set_label(new_label)

    ax.legend(title='Region', handles=handles)

    outdir = f'./output/{outtag}'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    outpath = pjoin(outdir, f'met_{year}_{distribution}.pdf')
    fig.savefig(outpath)
    plt.close(fig)
    print(f'File saved: {outpath}')

def plot_sigma_eta_phi(acc, outtag, distribution, year=2017, region='sr_vbf_hfhf'):
    '''Plot 2D sigma eta vs. phi distribution in the given region, for the leading or trailing jet.'''
    acc.load(distribution)
    h = acc[distribution]

    h = merge_extensions(h, acc, reweight_pu=False)
    h = merge_datasets(h)

    h = h.integrate('dataset', f'MET_{year}').integrate('region', region)

    fig, ax = plt.subplots()
    hist.plot2d(h, ax=ax, xaxis='sigmaetaeta', patch_opts={'norm': colors.LogNorm(1e-3,1e1)})
    
    ax.text(0., 1., f'MET UL v8 {year}',
        fontsize=14,
        ha='left',
        va='bottom',
        transform=ax.transAxes
    )

    regiontag = {
        'sr_vbf' : 'VBF SR',
        'sr_vbf_hfhf' : 'VBF SR (HF-HF)',
    }

    ax.text(1., 1., regiontag[region],
        fontsize=14,
        ha='right',
        va='bottom',
        transform=ax.transAxes
    )

    jettag = 'Leading Jet' if distribution == 'ak4_sigma_eta_phi0' else 'Trailing Jet'

    ax.text(0.05, 0.9, jettag,
        fontsize=14,
        ha='left',
        va='bottom',
        transform=ax.transAxes
    )

    # Save figure
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

    regions = [
        'sr_vbf',
        'sr_vbf_hfhf'
    ]
    for distribution in distributions:
        plot_hf_distributions(acc, outtag, distribution=distribution)

    # 2D eta vs. phi histograms
    distributions_2d = [
        'ak4_sigma_eta_phi0',
        'ak4_sigma_eta_phi1'
    ]

    for distribution in distributions_2d:
        plot_sigma_eta_phi(acc, outtag, distribution=distribution)

if __name__ == '__main__':
    main()