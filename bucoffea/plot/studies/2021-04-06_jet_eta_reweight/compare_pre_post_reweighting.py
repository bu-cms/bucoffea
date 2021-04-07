#!/usr/bin/env python
import os
import re
import sys
import numpy as np

from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
from coffea import hist
from bucoffea.plot.util import merge_datasets, merge_extensions, scale_xs_lumi, fig_ratio
from bucoffea.plot.style import matplotlib_rc
from bucoffea.helpers.paths import bucoffea_path
from klepto.archives import dir_archive

pjoin = os.path.join

matplotlib_rc()

def preprocess(h, acc, region, distribution, year=2017):
    h = merge_extensions(h, acc, reweight_pu=False)
    scale_xs_lumi(h)
    h = merge_datasets(h)

    h = h.integrate('region', region)

    if distribution == 'mjj':
        mjj_ax = hist.Bin('mjj', r'$M_{jj} \ (GeV)$', [200., 400., 600., 900., 1200., 1500., 2000., 2750., 3500., 5000.])
        h = h.rebin('mjj',mjj_ax)

    # Get total MC bkg in SR
    mc = re.compile(f'(ZJetsToNuNu.*|EW.*|Top_FXFX.*|Diboson.*|.*DYJetsToLL_M-50_HT_MLM.*|.*WJetsToLNu.*HT.*).*{year}')
    h = h.integrate('dataset', mc)

    return h

def compare_pre_post_reweighting(acc_dict, outtag, distribution='mjj', region='sr_vbf'):
    '''Compare pre and post reweighting distributions of the total bkg.'''
    histos = {}
    for key, acc in acc_dict.items():
        acc.load(distribution)
        histos[key] = preprocess(acc[distribution], acc, region, distribution)

    fig, ax, rax = fig_ratio()
    hist.plot1d(histos['withReweighting'], ax=ax)
    hist.plot1d(histos['noReweighting'], ax=ax, clear=False)

    ax.set_yscale('log')
    ax.set_ylim(1e0,1e6)

    ax.yaxis.set_ticks_position('both')

    ax.legend(labels=[
        'With Reweighting',
        'No Reweighting',
    ])

    ax.text(0.,1.,'ReReco MC Bkg 2017',
        fontsize=14,
        ha='left',
        va='bottom',
        transform=ax.transAxes
    )

    ax.text(1.,1.,'VBF SR',
        fontsize=14,
        ha='right',
        va='bottom',
        transform=ax.transAxes
    )

    data_err_opts = {
        'linestyle':'none',
        'marker': '.',
        'markersize': 10.,
        'color':'k',
    }

    hist.plotratio(
        histos['withReweighting'],
        histos['noReweighting'],
        ax=rax,
        unc='num',
        error_opts=data_err_opts
    )

    rax.grid(True)
    rax.set_ylim(0.6,1.4)
    rax.set_ylabel('With Rw / Nominal')

    loc = MultipleLocator(0.2)
    rax.yaxis.set_major_locator(loc)

    new_xlabels = {
        'ak4_eta0': r'Leading Jet $\eta$',
        'ak4_eta1': r'Trailing Jet $\eta$',
    }

    if distribution in new_xlabels.keys():
        ax.set_xlabel(new_xlabels[distribution])
        rax.set_xlabel(new_xlabels[distribution])

    outdir = f'./output/{outtag}'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outpath = pjoin(outdir, f'{region}_{distribution}.pdf')
    fig.savefig(outpath)
    plt.close(fig)

    print(f'File saved: {outpath}')

def main():
    acc_dict = {
        'withReweighting' : dir_archive( bucoffea_path('submission/merged_2021-04-06_vbfhinv_03Sep20v7_MC_reweightWithJetEta') ),
        'noReweighting' : dir_archive( bucoffea_path('submission/merged_2021-04-06_vbfhinv_03Sep20v7_MC_noReweighting_JetPtCut80') ),
    }

    for acc in acc_dict.values():
        acc.load('sumw')
        acc.load('sumw2')

    outtag = '06Apr21'

    for distribution in ['mjj', 'ak4_eta0', 'ak4_eta1']:
        compare_pre_post_reweighting(acc_dict, outtag, distribution=distribution)

if __name__ == '__main__':
    main()