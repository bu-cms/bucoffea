#!/usr/bin/env python

import os
import sys
import re
import numpy as np
import matplotlib.colors as colors

from coffea import hist
from scipy.stats import distributions
from bucoffea.plot.util import merge_datasets, merge_extensions, scale_xs_lumi
from matplotlib import pyplot as plt
from klepto.archives import dir_archive
from pprint import pprint

pjoin = os.path.join

def plot_eta_vs_phi(acc, outtag, distribution, region, dataset='MET_2017', plot_diag=True, xmax=0.3, ymax=0.3):
    '''2D sigma eta vs. phi plot for the given region.'''
    acc.load(distribution)
    h = acc[distribution]

    h = merge_extensions(h, acc, reweight_pu=False)
    scale_xs_lumi(h)
    h = merge_datasets(h)

    h = h.integrate('dataset', dataset).integrate('region', region)

    fig, ax = plt.subplots()
    hist.plot2d(h, ax=ax, xaxis='sigmaetaeta',  patch_opts={'norm': colors.LogNorm(1e-3,1e1)})

    ax.set_xlim(0,xmax)
    ax.set_ylim(0,ymax)

    if re.match('cr_2m_vbf.*', region):
        regiontag = r'$Z(\mu\mu)$'
    elif region == 'cr_vbf_qcd':
        regiontag = 'QCD CR'
    else:
        raise RuntimeError(f'Check region: {region}')

    ax.text(0., 1., f'{regiontag}, Data',
        fontsize=14,
        ha='left',
        va='bottom',
        transform=ax.transAxes
    )

    ax.text(1., 1., r'HF Jets, $p_T > 80 \ GeV$',
        fontsize=14,
        ha='right',
        va='bottom',
        transform=ax.transAxes
    )

    if distribution == 'ak4_sigma_eta_phi0':
        ax.set_xlabel(r'Leading Jet $\sigma_{\eta\eta}$')
        ax.set_ylabel(r'Leading Jet $\sigma_{\phi\phi}$')

    elif distribution == 'ak4_sigma_eta_phi1':
        ax.set_xlabel(r'Trailing Jet $\sigma_{\eta\eta}$')
        ax.set_ylabel(r'Trailing Jet $\sigma_{\phi\phi}$')

    if plot_diag:
        x = np.linspace(0, xmax)
        y = np.linspace(0, ymax)
        ax.plot(x,y,color='k',lw=2, label='Diagonal')

        slope = 0.5
        ycut = slope * x
        ax.plot(x,ycut,color='red',lw=2,label=f'$\\sigma_{{\phi\phi}} = {slope:.1f} \\sigma_{{\\eta\\eta}}$')

    ax.legend()

    # Save figure
    outdir = f'./output/{outtag}/2d'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    outpath = pjoin(outdir, f'MET_2017_{region}_{distribution}.pdf')
    fig.savefig(outpath)
    plt.close(fig)
    print(f'File saved: {outpath}')

def main():
    inpath = sys.argv[1]
    acc = dir_archive(inpath)
    acc.load('sumw')
    acc.load('sumw2')

    outtag = re.findall('merged_.*', inpath)[0].replace('/','')

    # 2D sigma eta/phi distributions for data and MC
    distributions = [
        'ak4_sigma_eta_phi',
        'ak4_sigma_eta_phi0',
        'ak4_sigma_eta_phi1',
    ]

    # Two regions: Z(mumu) physics-enriched and QCD CR
    regions = [
        'cr_2m_vbf_no_noisecuts',
        'cr_2m_vbf_relaxed_sel',
        'cr_vbf_qcd',
    ]

    for region in regions:
        try:
            for distribution in distributions:
                plot_eta_vs_phi(acc, outtag, 
                        distribution=distribution,
                        region=region
                        )
        except KeyError:
            print(f'Region not found in this input: {region}, moving on.')
            continue

if __name__ == '__main__':
    main()