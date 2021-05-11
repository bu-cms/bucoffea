#!/usr/bin/env python

import os
import sys
import re
import uproot
import numpy as np
import pandas as pd
import mplhep as hep

from bucoffea.plot.util import fig_ratio
from coffea.hist import poisson_interval
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
from pprint import pprint

pjoin = os.path.join

LABELS = {
    'ewk_zvv' : r'EWK $Z(\nu\nu)$',
    'ewk_zee' : r'EWK $Z(ee)$',
    'ewk_zmm' : r'EWK $Z(\mu\mu)$',
    'ewk_wminusmv' : r'EWK $W^{-}(\mu\nu)$',
    'ewk_wminusev' : r'EWK $W^{-}(e\nu)$',
    'ewk_wplusmv' : r'EWK $W^{+}(\mu\nu)$',
    'ewk_wplusev' : r'EWK $W^{+}(e\nu)$',
}

BINNINGS = {
    'mjj' : [200., 400., 600., 900., 1200., 1500., 2000., 2750., 3500., 5000.],
    'leadak4_ptnano' : list(range(80,1000,20)),
    'leadak4_eta' : np.linspace(-5,5),
    'trailak4_ptnano' : list(range(40,800,20)),
    'trailak4_eta' : np.linspace(-5,5),
    'recoil_pt_nano' : [ 250,  280,  310,  340,  370,  400,  430,  470,  510, 550,  590,  640,  690,  740,  790,  840,  900,  960, 1020, 1090, 1160, 1250, 1400],
    'recoil_phi_nano' : np.linspace(-np.pi, np.pi)
}

XLABELS = {
    'mjj' : r'$M_{jj} \ (GeV)$',
    'leadak4_ptnano' : r'Leading Jet $p_T \ (GeV)$',
    'leadak4_eta' : r'Leading Jet $\eta$',
    'trailak4_ptnano' : r'Trailing Jet $p_T \ (GeV)$',
    'trailak4_eta' : r'Trailing Jet $\eta$',
    'recoil_pt_nano' : r'Recoil $p_T \ (GeV)$',
    'recoil_phi_nano' : r'Recoil $\phi$',
}

XSECTIONS = {
    'ewk_wminusmv' : 32.27,
    'ewk_wminusev' : 32.27,
    'ewk_wplusmv' : 39.33,
    'ewk_wplusev' : 39.33,
    'ewk_zee' : 6.22,
    'ewk_zmm' : 6.22,
    'ewk_zvv' : 10.72,
}

# Total sumw and sumw2
SUMW = {
    'ewk_wminusmv' : 2571000.0,
    'ewk_wminusev' : 2571000.0,
    'ewk_wplusmv' : 1941000.0,
    'ewk_wplusev' : 1941000.0,
    'ewk_zee' : 997000.0,
    'ewk_zmm' : 997000.0,
    'ewk_zvv' : 2985000.0,
}

SUMW2 = {
    'ewk_wminusmv' : 2571000.0,
    'ewk_wminusev' : 2571000.0,
    'ewk_wplusmv' : 1941000.0,
    'ewk_wplusev' : 1941000.0,
    'ewk_zee' : 997000.0,
    'ewk_zmm' : 997000.0,
    'ewk_zvv' : 2985000.0,
}

def event_by_event_comparison(path_central, path_private, region, tag, distribution='mjj', year=2017):
    df_central = uproot.open(path_central)[region].pandas.df()
    df_private = uproot.open(path_private)[region].pandas.df()

    # Determine sample year
    year = re.findall('201\d', path_central)[0]

    # Merge dataframes based on run,event,lumi
    merged_df = df_central.merge(df_private, on=['run', 'event', 'lumi'], suffixes=['_central', '_private'])

    central_vals = merged_df[f'{distribution}_central']
    private_vals = merged_df[f'{distribution}_private']

    diff = (central_vals - private_vals) / (central_vals + private_vals)
    fig, ax = plt.subplots()
    bins= np.linspace(-0.05,0.05)
    ax.hist(diff, bins=bins)
    
    ax.text(0.,1.,LABELS[tag],
        fontsize=14,
        ha='left',
        va='bottom',
        transform=ax.transAxes
    )

    ax.text(1.,1.,year,
        fontsize=14,
        ha='right',
        va='bottom',
        transform=ax.transAxes
    )

    outdir = './output/11May21/with_event_match/deltas'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outpath = pjoin(outdir, f'{tag}_{distribution}_{year}.pdf')
    fig.savefig(outpath)
    plt.close(fig)
    print(f'File saved: {outpath}')


def make_comparison(path_central, path_private, region, tag, distribution='mjj'):
    df_central = uproot.open(path_central)[region].pandas.df()
    df_private = uproot.open(path_private)[region].pandas.df()

    # Merge dataframes based on run,event,lumi
    merged_df = df_central.merge(df_private, on=['run', 'event', 'lumi'], suffixes=['_central', '_private'])

    # Determine sample year
    year = re.findall('201\d', path_central)[0]

    h_central, _ = np.histogram(merged_df[f'{distribution}_central'], bins=BINNINGS[distribution])
    h_private, _ = np.histogram(merged_df[f'{distribution}_private'], bins=BINNINGS[distribution])
    # Scale with xs & lumi
    def scale_xs_lumi(h, tag):
        return h * XSECTIONS[tag] * 41.5 * 1e3 / SUMW[tag] 

    h_central = scale_xs_lumi(h_central, tag)
    h_private = scale_xs_lumi(h_private, tag)

    # Compute stat errors on templates
    def w2err(sumw, sumw2):
        return np.abs(poisson_interval(sumw, sumw2) - sumw)

    # Compute errors per mjj bin
    yerr_central = w2err(h_central, h_central)
    yerr_private = w2err(h_private, h_private)

    fig, ax, rax = fig_ratio()
    hep.histplot(h_central, bins=BINNINGS[distribution], ax=ax, yerr=yerr_central, histtype='step', label='Central')
    hep.histplot(h_private, bins=BINNINGS[distribution], ax=ax, yerr=yerr_private, histtype='step', label='Private')

    ax.legend(title='Production')
    ax.set_yscale('log')
    ax.set_ylim(1e0,1e4)
    ax.set_ylabel('Events')
    ax.set_xlabel(XLABELS[distribution])
    
    ax.text(0.,1.,LABELS[tag],
        fontsize=14,
        ha='left',
        va='bottom',
        transform=ax.transAxes
    )

    ax.text(1.,1.,year,
        fontsize=14,
        ha='right',
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

    xedges = BINNINGS[distribution]

    r = h_central / h_private
    mjj_centers = 0.5 * (xedges + np.roll(xedges, -1))[:-1]

    rax.plot(mjj_centers, r, **data_err_opts)

    rax.grid(True)
    rax.set_ylim(0.9,1.1)
    rax.set_ylabel('Central / Private')
    rax.set_xlabel(XLABELS[distribution])

    outdir = './output/11May21/with_event_match'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outpath = pjoin(outdir, f'{tag}_{distribution}_{year}.pdf')
    fig.savefig(outpath)
    plt.close(fig)
    print(f'File saved: {outpath}')

def main():
    treetag = 'trees_UL_11May21'
    for year in [2017,2018]:
        treepaths = {
            'ewk_zvv' : {
                'path_central': f'./input/{treetag}/central/tree_EWKZ2Jets_ZToNuNu_M-50_withDipoleRecoil-mg_{year}.root',
                'path_private': f'./input/{treetag}/private/tree_EWKZ2Jets_ZToNuNu_M-50_withDipoleRecoil-mg_{year}.root',
                'region': 'sr_vbf',
            },
            'ewk_zmm' : {
                'path_central': f'./input/{treetag}/central/tree_EWKZ2Jets_ZToLL_M-50_withDipoleRecoil-mg_{year}.root',
                'path_private': f'./input/{treetag}/private/tree_EWKZ2Jets_ZToLL_M-50_withDipoleRecoil-mg_{year}.root',
                'region' : 'cr_2m_vbf',
            },
            'ewk_zee' : {
                'path_central': f'./input/{treetag}/central/tree_EWKZ2Jets_ZToLL_M-50_withDipoleRecoil-mg_{year}.root',
                'path_private': f'./input/{treetag}/private/tree_EWKZ2Jets_ZToLL_M-50_withDipoleRecoil-mg_{year}.root',
                'region' : 'cr_2e_vbf',
            },
            'ewk_wminusmv' : {
                'path_central': f'./input/{treetag}/central/tree_EWKWMinus2Jets_WToLNu_M-50_withDipoleRecoil-mg_{year}.root',
                'path_private': f'./input/{treetag}/private/tree_EWKWMinus2Jets_WToLNu_M-50_withDipoleRecoil-mg_{year}.root',
                'region' : 'cr_1m_vbf',
            },
            'ewk_wplusmv' : {
                'path_central': f'./input/{treetag}/central/tree_EWKWPlus2Jets_WToLNu_M-50_withDipoleRecoil-mg_{year}.root',
                'path_private': f'./input/{treetag}/private/tree_EWKWPlus2Jets_WToLNu_M-50_withDipoleRecoil-mg_{year}.root',
                'region' : 'cr_1m_vbf',
            },
            'ewk_wminusev' : {
                'path_central': f'./input/{treetag}/central/tree_EWKWMinus2Jets_WToLNu_M-50_withDipoleRecoil-mg_{year}.root',
                'path_private': f'./input/{treetag}/private/tree_EWKWMinus2Jets_WToLNu_M-50_withDipoleRecoil-mg_{year}.root',
                'region' : 'cr_1e_vbf',
            },
            'ewk_wplusev' : {
                'path_central': f'./input/{treetag}/central/tree_EWKWPlus2Jets_WToLNu_M-50_withDipoleRecoil-mg_{year}.root',
                'path_private': f'./input/{treetag}/private/tree_EWKWPlus2Jets_WToLNu_M-50_withDipoleRecoil-mg_{year}.root',
                'region' : 'cr_1e_vbf',
            },
        }
    
        distributions = [
            'leadak4_ptnano',
            'leadak4_eta',
            'trailak4_ptnano',
            'trailak4_eta',
            'recoil_pt_nano',
            'recoil_phi_nano',
        ]
    
        tags  = list(treepaths.keys())
        for tag in tags:
            if year == 2018 and tag == 'ewk_zvv':
                continue
            for distribution in distributions:
                # Compare the event-matched distributions
                make_comparison(
                    treepaths[tag]['path_central'],
                    treepaths[tag]['path_private'],
                    treepaths[tag]['region'],
                    tag=tag,
                    distribution=distribution,
                )
                # Plot the difference between the two quantities
                event_by_event_comparison(
                    treepaths[tag]['path_central'],
                    treepaths[tag]['path_private'],
                    treepaths[tag]['region'],
                    tag=tag,
                    distribution=distribution
                )

if __name__ == '__main__':
    main()