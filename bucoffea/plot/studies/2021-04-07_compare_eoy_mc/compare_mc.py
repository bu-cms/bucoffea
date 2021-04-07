#!/usr/bin/env python
import os
import re
import sys
import numpy as np

from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
from coffea import hist
from bucoffea.plot.util import merge_datasets, merge_extensions, scale_xs_lumi, fig_ratio
from bucoffea.helpers.paths import bucoffea_path
from klepto.archives import dir_archive

pjoin = os.path.join

def preprocess(h, acc, region, distribution, year=2017, plot='total_bkg'):
    h = merge_extensions(h, acc, reweight_pu=False)
    scale_xs_lumi(h)
    h = merge_datasets(h)

    h = h.integrate('region', region)

    if distribution == 'mjj':
        mjj_ax = hist.Bin('mjj', r'$M_{jj} \ (GeV)$', [200., 400., 600., 900., 1200., 1500., 2000., 2750., 3500., 5000.])
        h = h.rebin('mjj',mjj_ax)

    # Get total MC bkg in SR
    if plot == 'total_bkg':
        mc = re.compile(f'(ZJetsToNuNu.*|EW.*|Top_FXFX.*|Diboson.*|DYJetsToLL_M-50_HT_MLM.*|WJetsToLNu.*HT.*).*{year}')
    elif plot == 'signal':
        mc = re.compile(f'VBF_HToInv.*M125.*{year}')
    else:
        raise RuntimeError(f'Invalid parameter name for plot: {plot}')
    h = h.integrate('dataset', mc)

    return h

def compare_mc(acc_dict, outtag, distribution='mjj', region='sr_vbf_no_veto_all', year=2017, plot='total_bkg'):
    '''
    Compare total ReReco background or signal for three cases:
    1. No HF or endcap reweighting
    2. With HF reweighting
    3. With HF + endcap reweighting
    '''
    histos = {}
    for key, acc in acc_dict.items():
        acc.load(distribution)
        histos[key] = preprocess(acc[distribution], acc, region, distribution, year, plot)

    fig, ax, rax = fig_ratio()
    hist.plot1d(histos['noReweighting'], ax=ax)
    hist.plot1d(histos['withHFReweighting'], ax=ax, clear=False)
    hist.plot1d(histos['withEndcapAndHFReweighting'], ax=ax, clear=False)

    ax.set_yscale('log')
    ax.set_ylim(1e0,1e6)
    ax.yaxis.set_ticks_position('both')

    ax.legend(title='Additional Reweighting', labels=[
        'None',
        'HF',
        'Endcap + HF'  
    ])

    plot_to_text = {
        'total_bkg' : 'ReReco MC Bkg 2017',
        'signal' : 'VBF H(inv) 2017',
    }

    ax.text(0.,1., plot_to_text[plot],
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
    }

    h_nom = histos['noReweighting']

    legend_labels = {
        'withHFReweighting' : 'HF rw',
        'withEndcapAndHFReweighting' : 'Endcap + HF rw',
    }

    for key in ['withHFReweighting', 'withEndcapAndHFReweighting']:
        hist.plotratio(
            histos[key],
            h_nom,
            ax=rax,
            unc='num',
            error_opts=data_err_opts,
            clear=False,
            label=legend_labels[key]
        )

    rax.grid(True)
    if distribution == 'mjj':
        rax.set_ylim(0.5,1.5)
    else:
        rax.set_ylim(0,2)
    rax.set_ylabel('Ratio to Nom.')

    rax.legend(ncol=2)

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
    outpath = pjoin(outdir, f'{plot}_{region}_{distribution}.pdf')
    fig.savefig(outpath)
    plt.close(fig)

    print(f'File saved: {outpath}')

def main():
    acc_dict = {
        'noReweighting': dir_archive(bucoffea_path('submission/merged_2021-04-07_vbfhinv_03Sep20v7_MC_noEndcapReweight_noHFReweight_jetPt80')),
        'withHFReweighting': dir_archive(bucoffea_path('submission/merged_2021-04-07_vbfhinv_03Sep20v7_MC_noEndcapReweight_withHFReweight_jetPt80')),
        'withEndcapAndHFReweighting': dir_archive(bucoffea_path('submission/merged_2021-04-07_vbfhinv_03Sep20v7_MC_withEndcapReweight_withHFReweight_jetPt80')),
    }

    for acc in acc_dict.values():
        acc.load('sumw')
        acc.load('sumw2')

    outtag = '07Apr21'

    # Make the comparison for VBF signal and total background
    for plot in ['signal', 'total_bkg']:
        for distribution in ['mjj', 'ak4_eta0', 'ak4_eta1']:
            compare_mc(acc_dict, outtag, distribution=distribution, plot=plot)

if __name__ == '__main__':
    main()