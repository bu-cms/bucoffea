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

def plot_sigma_eta_phi_thirdjet(acc, outtag, region='cr_vbf_qcd', dataset='MET_2017'):
    '''Plot 2D sigma eta vs. phi for third jets in HF (if any).'''
    distribution = 'ak4_sigma_eta_phi_third_jet'
    acc.load(distribution)
    h = acc[distribution]

    h = merge_extensions(h, acc, reweight_pu=False)
    scale_xs_lumi(h)
    h = merge_datasets(h)

    h = h.integrate('dataset', dataset).integrate('region', region)

    fig, ax = plt.subplots()
    hist.plot2d(h, ax=ax, xaxis='sigmaetaeta',  patch_opts={'norm': colors.LogNorm(1e-3,1e1)})

    ax.set_xlim(0,0.3)
    ax.set_ylim(0,0.3)

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
    
    ax.text(1., 1., r'Third Jet in HF ($p_T > 80 \ GeV$)',
        fontsize=12,
        ha='right',
        va='bottom',
        transform=ax.transAxes
    )
    
    x = np.linspace(0, 0.3)
    y = np.linspace(0, 0.3)
    ax.plot(x,y,color='k',lw=2)

    ax.set_xlabel(r'Third Jet $\sigma_{\eta\eta}$')
    ax.set_ylabel(r'Third Jet $\sigma_{\phi\phi}$')

    # Save figure
    outdir = f'./output/{outtag}/2d/third_jet'
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

    plot_sigma_eta_phi_thirdjet(acc, outtag)

if __name__ == '__main__':
    main()