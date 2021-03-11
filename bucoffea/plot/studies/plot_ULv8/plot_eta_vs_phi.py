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

def plot_eta_vs_phi(acc, outtag, distribution, region, dataset='MET_2017', plot_diag=True, plot_cutslope=False, mindphislice=None, etaslice=None, xmax=0.3, ymax=0.3):
    '''2D sigma eta vs. phi plot for the given region.'''
    acc.load(distribution)
    h = acc[distribution]

    h = merge_extensions(h, acc, reweight_pu=False)
    scale_xs_lumi(h)
    h = merge_datasets(h)

    h = h.integrate('dataset', dataset).integrate('region', region)

    # If the histogram is also split by jet eta, take the slice we're interested in
    # Typically: 2.9 < |eta| < 3.25 and 3.25 < |eta| < 5.0
    if etaslice is not None:
        h = h.integrate('jeta', etaslice)

    if mindphislice is not None:
        h = h.integrate('mindphi', mindphislice)
    
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

    if etaslice in [slice(2.9, 3.25), slice(3.25, 5.0)]:
        tagforslice = {
            (2.9, 3.25): r'$2.9 < |\eta| < 3.25$',
            (3.25, 5.0): r'$3.25 < |\eta| < 5.0$',
        }
        
        for (lo, hi), _tag in tagforslice.items():
            if (etaslice.start, etaslice.stop) == (lo, hi):
                _tagforslice = _tag
                break

        ax.text(1., 0., _tagforslice,
            fontsize=14,
            ha='right',
            va='bottom',
            transform=ax.transAxes
        )

    elif mindphislice is not None:
        if mindphislice == slice(0, 0.1):
            tagforslice = f'$\\Delta\\phi(j,MET) < 0.1$'
        else:
            tagforslice = f'${mindphislice.start:.1f} < \\Delta\\phi(j,MET) < {mindphislice.stop:.1f}$'

        ax.text(1., 0., tagforslice,
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

        if plot_cutslope:
            slope = 0.5
            ycut = slope * x
            ax.plot(x,ycut,color='red',lw=2,label=f'$\\sigma_{{\phi\phi}} = {slope:.1f} \\sigma_{{\\eta\\eta}}$')

    ax.legend()

    # Save figure
    outdir = f'./output/{outtag}/2d'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    if etaslice in [slice(2.9, 3.25), slice(3.25, 5.0)]:
        etaslicetag = f'_eta_{str(etaslice.start).replace(".", "p")}_{str(etaslice.stop).replace(".", "p")}'
    elif etaslice == slice(2.9, 5.0):
        etaslicetag= '_eta_combined'
    else:
        etaslicetag = ''

    if mindphislice in [slice(0,0.1),slice(0.1,0.3),slice(0.3,0.5)]:
        mindphislicetag = f'_mindphi_{str(mindphislice.start).replace(".","p")}_{str(mindphislice.stop).replace(".","p")}'
    else:
        mindphislicetag = ''

    outpath = pjoin(outdir, f'MET_2017_{region}_{distribution}{etaslicetag}{mindphislicetag}.pdf')
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
        # 'cr_2m_vbf_relaxed_sel',
        'cr_vbf_qcd',
    ]

    etaslices = [
        slice(2.9, 3.25),
        slice(3.25, 5.0),
        slice(2.9, 5.0),
    ]

    # By default, no mindphi slices
    # List these if you need it
    mindphislices = {
        'ak4_sigma_eta_phi' : [
            None
        ],    
        'ak4_sigma_eta_phi0' : [
            None
        ],    
        'ak4_sigma_eta_phi1' : [
            None
        ],    
    }

    for region in regions:
        try:
            for distribution in distributions:
                for mindphislice in mindphislices[distribution]:
                    for etaslice in etaslices:
                        plot_eta_vs_phi(acc, outtag, 
                                distribution=distribution,
                                region=region,
                                etaslice=etaslice,
                                mindphislice=mindphislice,
                                plot_cutslope=True
                                )
        except KeyError:
            print(f'Region not found in this input: {region}, moving on.')
            continue

if __name__ == '__main__':
    main()