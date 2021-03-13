#!/usr/bin/env python

import os
import sys
import re
import uproot
import numpy as np

from coffea import hist
from bucoffea.plot.util import merge_datasets, merge_extensions, scale_xs_lumi
from matplotlib import pyplot as plt
from klepto.archives import dir_archive
from pprint import pprint

pjoin = os.path.join

def plot_efficiency(acc, outtag, distribution='ak4_pt0_eta0_hf', region='noise_enriched', dataset='MET_2017', etaslice=slice(2.9,3.25), rebin_pt=True):
    acc.load(distribution)
    h = acc[distribution]

    h = merge_extensions(h, acc, reweight_pu=False)
    scale_xs_lumi(h)
    h = merge_datasets(h)

    # Get the dataset and the regions we're interested in
    if region == 'noise_enriched':
        regionbase = 'cr_vbf_qcd'
    elif region == 'physics_enriched':
        regionbase = 'cr_2m_vbf_relaxed_sel'
    else:
        raise RuntimeError(f'Check region: {region}')

    regionregex = re.compile(f'{regionbase}.*')

    h = h.integrate('dataset', dataset).integrate('jeta', etaslice)[regionregex]

    # Rebin jet pt
    if rebin_pt:
        if 'ak4_pt' in distribution:
            new_ax = hist.Bin('jetpt', r'Leading Jet $p_T$ (GeV)', list(range(80,280,50)) + list(range(280,480,100)) + [600])
            h = h.rebin('jetpt', new_ax)

    data_err_opts = {
        'linestyle':'none',
        'marker': '.',
        'markersize': 10.,
    }

    fig, ax = plt.subplots()
    # Plot the efficiencies
    try:
        hist.plotratio(
            h.integrate('region', f'{regionbase}_sphietacut'),
            h.integrate('region', f'{regionbase}'),
            ax=ax,
            error_opts=data_err_opts,
            label=r'$\sigma_{\phi\phi} / \sigma_{\eta\eta} > 0.5$ ($2.9 < |\eta| < 3.25$)'
        )
    except KeyError:
        pass

    try:
        hist.plotratio(
            h.integrate('region', f'{regionbase}_cssizecut'),
            h.integrate('region', f'{regionbase}'),
            ax=ax,
            error_opts=data_err_opts,
            label=r'CStripSize < 3',
            clear=False
        )
    except KeyError:
        pass
    
    hist.plotratio(
        h.integrate('region', f'{regionbase}_both_cuts'),
        h.integrate('region', f'{regionbase}'),
        ax=ax,
        error_opts=data_err_opts,
        label=r'Both cuts applied',
        clear=False
    )

    ax.legend()

    if 'ak4_pt' in distribution:
        ax.set_ylim(0., 1.1)
    elif 'ak4_eta' in distribution:
        ax.set_ylim(0., 1.1)
    
    ax.set_ylabel('Efficiency')
    if distribution == 'ak4_eta0':
        ax.set_xlabel(r'Leading Jet $\eta$')

    ax.axhline(1, xmin=0, xmax=1, color='k')

    regiontags = {
        'noise_enriched' : r'QCD CR',
        'physics_enriched' : r'$Z(\mu\mu)$',
    }

    ax.text(0., 1., regiontags[region],
        fontsize=14,
        ha='left',
        va='bottom',
        transform=ax.transAxes
    )

    plottag = f'${etaslice.start:.2f} < |\\eta| < {etaslice.stop:.2f}$'
    ax.text(1., 1., plottag,
        fontsize=14,
        ha='right',
        va='bottom',
        transform=ax.transAxes
    )
    

    if etaslice in [slice(2.9, 3.25), slice(3.25, 5.0)]:
        etaslicetag = f'_eta_{str(etaslice.start).replace(".", "p")}_{str(etaslice.stop).replace(".", "p")}'
    elif etaslice == slice(2.9,5.0):
        etaslicetag = ''

    # Save figure
    outdir = f'./output/{outtag}/eff'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outpath = pjoin(outdir, f'{region}_{distribution}{etaslicetag}.pdf')
    fig.savefig(outpath)
    plt.close(fig)

    print(f'File saved: {outpath}')

def main():
    inpath = sys.argv[1]
    acc = dir_archive(inpath)
    acc.load('sumw')
    acc.load('sumw2')

    outtag = re.findall('merged_.*', inpath)[0].replace('/','')

    regions = [
        'physics_enriched',
        'noise_enriched',
    ]

    etaslices = [
        slice(2.9, 3.25),
        slice(3.25, 5.0),
        slice(2.9,5.0)
    ]

    for region in regions:
        for etaslice in etaslices:
            plot_efficiency(acc, outtag, 
                    etaslice=etaslice, 
                    region=region
                    )

if __name__ == '__main__':
    main()
