#!/usr/bin/env python

from io import UnsupportedOperation
import os
import sys
import re
import uproot
import numpy as np
import mplhep as hep

from coffea import hist
from bucoffea.plot.util import merge_datasets, merge_extensions, scale_xs_lumi, fig_ratio
from matplotlib import pyplot as plt
from klepto.archives import dir_archive
from pprint import pprint

pjoin = os.path.join

def get_qcd_estimate(acc, outtag, outputrootfile, distribution='mjj'):
    '''Get HF-based QCD estimate from ReReco data.'''
    acc.load(distribution)
    h = acc[distribution]

    h = merge_extensions(h, acc, reweight_pu=False)
    scale_xs_lumi(h)
    h = merge_datasets(h)

    if distribution == 'mjj':
        mjj_ax = hist.Bin('mjj', r'$M_{jj} \ (GeV)$', [200., 400., 600., 900., 1200., 1500., 2000., 2750., 3500., 5000.])
        h = h.rebin('mjj', mjj_ax)

    outdir = f'./output/{outtag}'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    for year in [2017, 2018]:
        # Get data in CR, scaled by the transfer factor
        try:
            h_data = h.integrate('dataset', f'MET_{year}').integrate('region', 'sr_vbf_qcd_data')
        except KeyError:
            print(f'WARNING: Data input for year {year} not found, skipping')
            continue

        # Total MC in SR, scaled by mistag rate * TF
        mc_regex = re.compile(f'(ZJetsToNuNu.*|EW.*|Top_FXFX.*|Diboson.*|DYJetsToLL_M-50_HT_MLM.*|WJetsToLNu.*HT.*).*{year}')
        h_mc = h.integrate('dataset', mc_regex).integrate('region', 'sr_vbf_qcd_mc')

        # Let's plot these separately and then plot the diff
        fig, ax, rax = fig_ratio()
        hist.plot1d(h_data, ax=ax)
        hist.plot1d(h_mc, ax=ax, clear=False)

        ax.text(0.,1.,year,
            fontsize=14,
            ha='left',
            va='bottom',
            transform=ax.transAxes
        )

        ax.text(1.,1.,r'$TF = P(pass | noise) \ / \ P(fail | noise)$',
            fontsize=14,
            ha='right',
            va='bottom',
            transform=ax.transAxes
        )

        ax.legend(labels=[
            '1: (Data in CR x TF)',
            '2: (MC in SR x TF x Mistag Rate)',
        ])

        ax.set_yscale('log')
        ax.set_ylim(1e-2,1e3)

        # Our estimate is: h_data - h_mc
        h_mc.scale(-1)
        h_data.add(h_mc)
        hist.plot1d(h_data, ax=rax)

        rax.legend(labels=['QCD Estimate: 1 - 2'])
        rax.set_yscale('log')
        rax.set_ylim(1e-2,1e3)

        ax.yaxis.set_ticks_position('both')
        rax.yaxis.set_ticks_position('both')

        outpath = pjoin(outdir,f'qcd_estimate_{distribution}_{year}.pdf')
        fig.savefig(outpath)
        plt.close(fig)

        print(f'File saved: {outpath}')

        # Save to output ROOT file
        outputrootfile[f'sr_template_{year}_{distribution}'] = (h_data.values()[()], h_data.axes()[0].edges())

def main():
    inpath = sys.argv[1]
    acc = dir_archive(inpath)
    acc.load('sumw')
    acc.load('sumw2')

    outtag = re.findall('merged_.*', inpath)[0].replace('/','')

    outdir = f'./output/{outtag}'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    outputrootpath = pjoin(outdir, 'qcd_estimate_sr.root')
    outputrootfile = uproot.recreate(outputrootpath)

    for distribution in ['mjj', 'ak4_eta0', 'ak4_eta1']:
        get_qcd_estimate(acc, outtag, outputrootfile=outputrootfile, distribution=distribution)

if __name__ == '__main__':
    main()