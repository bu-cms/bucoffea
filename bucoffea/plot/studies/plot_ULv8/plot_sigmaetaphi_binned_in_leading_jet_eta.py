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

def plot_sigma_eta_phi_binned_in_leading_jet_eta(acc, outtag, distribution, region, etaslice, dataset='MET_2017', xmax=0.3, ymax=0.3):
    acc.load(distribution)
    h = acc[distribution]

    h = merge_extensions(h, acc, reweight_pu=False)
    scale_xs_lumi(h)
    h = merge_datasets(h)

    h = h.integrate('dataset', dataset).integrate('region', region)
    h = h.integrate('jeta', etaslice)

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
        fontsize=12,
        ha='left',
        va='bottom',
        transform=ax.transAxes
    )

    ax.text(1., 1., f'Leading Jet: ${etaslice.start:.1f} < |\eta| < {etaslice.stop:.1f}$',
        fontsize=12,
        ha='right',
        va='bottom',
        transform=ax.transAxes
    )

    if distribution == 'ak4_sigma_eta_phi1_binned_in_ak4_eta0':
        ax.set_xlabel(r'Trailing Jet $\sigma_{\eta\eta}$')
        ax.set_ylabel(r'Trailing Jet $\sigma_{\phi\phi}$')
    elif distribution == 'ak4_sigma_eta_phi_binned_in_ak4_eta0':
        ax.set_xlabel(r'All Jet $\sigma_{\eta\eta}$')
        ax.set_ylabel(r'All Jet $\sigma_{\phi\phi}$')
    else:
        raise RuntimeError(f'Check distribution: {distribution}')

    x = np.linspace(0, xmax)
    y = np.linspace(0, ymax)
    ax.plot(x,y,color='k',lw=2)

    # Save figure
    outdir = f'./output/{outtag}/2d/binned_in_ak4_eta0'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    etaslicetag = f'_eta_{str(etaslice.start).replace(".", "p")}_{str(etaslice.stop).replace(".", "p")}'

    outpath = pjoin(outdir, f'MET_2017_{region}_{distribution}{etaslicetag}.pdf')
    fig.savefig(outpath)
    plt.close(fig)
    print(f'File saved: {outpath}')

def main():
    inpath = sys.argv[1]
    acc = dir_archive(inpath)
    acc.load('sumw')
    acc.load('sumw2')

    outtag = re.findall('merged_.*', inpath)[0].replace('/','')

    # Eta slices for the leading jet: Central or forward (or both combined)
    etaslices = [
        slice(0, 2.5),
        slice(2.5, 5),
        slice(0., 5),
    ]

    # Sigma eta vs. phi of all jets vs the trailing jets
    distributions = [
        'ak4_sigma_eta_phi_binned_in_ak4_eta0',
        'ak4_sigma_eta_phi1_binned_in_ak4_eta0',
    ]

    for distribution in distributions:
        try:
            for etaslice in etaslices:
                plot_sigma_eta_phi_binned_in_leading_jet_eta(acc, outtag, 
                    distribution=distribution,
                    region='cr_vbf_qcd',
                    etaslice=etaslice    
                )
        except KeyError:
            print(f'Distribution not found: {distribution}, skipping')
            continue

if __name__ == '__main__':
    main()