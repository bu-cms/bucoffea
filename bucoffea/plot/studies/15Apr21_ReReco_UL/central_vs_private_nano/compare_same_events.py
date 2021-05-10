#!/usr/bin/env python

from bucoffea.plot.util import fig_ratio
import os
import sys
import re
import uproot
import numpy as np
import pandas as pd
import mplhep as hep

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

def make_comparison(path_central, path_private, region, tag, distribution='mjj'):
    df_central = uproot.open(path_central)[region].pandas.df()
    df_private = uproot.open(path_private)[region].pandas.df()

    # Merge dataframes based on run,event,lumi
    merged_df = df_central.merge(df_private, on=['run', 'event', 'lumi'], suffixes=['_central', '_private'])

    mjj_edges = [200., 400., 600., 900., 1200., 1500., 2000., 2750., 3500., 5000.]

    central_weights = merged_df['weight_total_central']
    private_weights = merged_df['weight_total_private']

    h_central, _ = np.histogram(merged_df[f'{distribution}_central'], bins=mjj_edges, weights=central_weights)
    h_private, _ = np.histogram(merged_df[f'{distribution}_private'], bins=mjj_edges, weights=private_weights)
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
    hep.histplot(h_central, ax=ax, bins=mjj_edges, yerr=yerr_central, histtype='step', label='Central')
    hep.histplot(h_private, ax=ax, bins=mjj_edges, yerr=yerr_private, histtype='step', label='Private')

    ax.legend(title='Production')
    ax.set_yscale('log')
    ax.set_ylim(1e0,1e4)
    ax.set_ylabel('Events')
    ax.set_xlabel(r'$M_{jj} \ (GeV)$')
    
    ax.text(0.,1.,LABELS[tag],
        fontsize=14,
        ha='left',
        va='bottom',
        transform=ax.transAxes
    )

    ax.text(1.,1.,'2017',
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

    r = h_central / h_private
    mjj_centers = 0.5 * (mjj_edges + np.roll(mjj_edges, -1))[:-1]

    rax.plot(mjj_centers, r, **data_err_opts)

    rax.grid(True)
    rax.set_ylim(0.9,1.1)
    rax.set_ylabel('Central / Private')
    rax.set_xlabel(r'$M_{jj} \ (GeV)$')

    outdir = './output/10May21/with_event_match'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outpath = pjoin(outdir, f'{tag}_{distribution}_2017.pdf')
    fig.savefig(outpath)
    plt.close(fig)
    print(f'File saved: {outpath}')

def main():
    treepaths = {
        'ewk_zvv' : {
            'path_central': './input/trees_10May21_EWK_V/central/tree_EWKZ2Jets_ZToNuNu_M-50_withDipoleRecoil-mg_2017.root',
            'path_private': './input/trees_10May21_EWK_V/private/tree_EWKZ2Jets_ZToNuNu_M-50_withDipoleRecoil-mg_2017.root',
            'region': 'sr_vbf',
        },
        'ewk_zmm' : {
            'path_central': './input/trees_10May21_EWK_V/central/tree_EWKZ2Jets_ZToLL_M-50_withDipoleRecoil-mg_2017.root',
            'path_private': './input/trees_10May21_EWK_V/private/tree_EWKZ2Jets_ZToLL_M-50_withDipoleRecoil-mg_2017.root',
            'region' : 'cr_2m_vbf',
        },
        'ewk_zee' : {
            'path_central': './input/trees_10May21_EWK_V/central/tree_EWKZ2Jets_ZToLL_M-50_withDipoleRecoil-mg_2017.root',
            'path_private': './input/trees_10May21_EWK_V/private/tree_EWKZ2Jets_ZToLL_M-50_withDipoleRecoil-mg_2017.root',
            'region' : 'cr_2e_vbf',
        },
        'ewk_wminusmv' : {
            'path_central': './input/trees_10May21_EWK_V/central/tree_EWKWMinus2Jets_WToLNu_M-50_withDipoleRecoil-mg_2017.root',
            'path_private': './input/trees_10May21_EWK_V/private/tree_EWKWMinus2Jets_WToLNu_M-50_withDipoleRecoil-mg_2017.root',
            'region' : 'cr_1m_vbf',
        },
        'ewk_wplusmv' : {
            'path_central': './input/trees_10May21_EWK_V/central/tree_EWKWPlus2Jets_WToLNu_M-50_withDipoleRecoil-mg_2017.root',
            'path_private': './input/trees_10May21_EWK_V/private/tree_EWKWPlus2Jets_WToLNu_M-50_withDipoleRecoil-mg_2017.root',
            'region' : 'cr_1m_vbf',
        },
        'ewk_wminusev' : {
            'path_central': './input/trees_10May21_EWK_V/central/tree_EWKWMinus2Jets_WToLNu_M-50_withDipoleRecoil-mg_2017.root',
            'path_private': './input/trees_10May21_EWK_V/private/tree_EWKWMinus2Jets_WToLNu_M-50_withDipoleRecoil-mg_2017.root',
            'region' : 'cr_1e_vbf',
        },
        'ewk_wplusev' : {
            'path_central': './input/trees_10May21_EWK_V/central/tree_EWKWPlus2Jets_WToLNu_M-50_withDipoleRecoil-mg_2017.root',
            'path_private': './input/trees_10May21_EWK_V/private/tree_EWKWPlus2Jets_WToLNu_M-50_withDipoleRecoil-mg_2017.root',
            'region' : 'cr_1e_vbf',
        },
    }

    tags  = list(treepaths.keys())
    for tag in tags:
        make_comparison(
            treepaths[tag]['path_central'],
            treepaths[tag]['path_private'],
            treepaths[tag]['region'],
            tag=tag
        )

if __name__ == '__main__':
    main()