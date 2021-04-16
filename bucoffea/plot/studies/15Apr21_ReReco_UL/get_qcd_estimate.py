#!/usr/bin/env python

import os
import sys
import re
import uproot
import numpy as np

from matplotlib import pyplot as plt
from bucoffea.plot.util import merge_extensions, merge_datasets, scale_xs_lumi
from coffea import hist
from klepto.archives import dir_archive
from pprint import pprint
from distributions import distributions, binnings

pjoin = os.path.join


def get_qcd_estimate(acc, outtag, outrootfile, distribution):
    '''Calculate the QCD template in SR'''
    acc.load(distribution)
    h = acc[distribution]

    h = merge_extensions(h, acc, reweight_pu=False)
    scale_xs_lumi(h)
    h = merge_datasets(h)

    overflow = 'none'
    if distribution == 'mjj':
        overflow = 'over'

    if distribution in binnings.keys():
        new_ax = binnings[distribution]
        h = h.rebin(new_ax.name, new_ax)
    
    outdir = f'./output/{outtag}/qcd_estimate'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    new_xlabels = {
        'ak4_eta0': r'Leading Jet $\eta$',
        'ak4_eta1': r'Trailing Jet $\eta$',
    }

    # Get data and MC yields in the QCD CR
    h = h.integrate('region', 'cr_vbf_qcd')
    for year in [2017, 2018]:
        data = f'MET_{year}'
        # mc = re.compile(f'(ZJetsToNuNu.*|EW.*|Top_FXFX.*|Diboson.*|DYJetsToLL_M-50_HT_MLM.*|WJetsToLNu.*HT.*).*{year}')
        mc = re.compile(f'(ZJetsToNuNu.*|EW.*|DYJetsToLL_M-50_HT_MLM.*|WJetsToLNu.*HT.*).*{year}')

        fig, ax = plt.subplots()
        hist.plot1d(h.integrate('dataset', data), ax=ax, overflow=overflow)
        hist.plot1d(h.integrate('dataset', mc), ax=ax, clear=False, overflow=overflow)

        ax.set_yscale('log')
        ax.set_ylim(1e-2,1e6)
        
        ax.text(0.,1.,r'QCD CR $\times$ $CR \rightarrow SR$ TF',
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

        ax.legend(labels=['Data', 'MC']) 

        if distribution in new_xlabels.keys():
            ax.set_xlabel(new_xlabels[distribution])

        outpath = pjoin(outdir, f'qcd_cr_{distribution}_{year}.pdf')
        fig.savefig(outpath)
        plt.close(fig)
        print(f'File saved: {outpath}')

        fig, ax = plt.subplots()
        # Calculate the QCD template in SR: Data - MC in CR (already weighted by TF)
        h_qcd = h.integrate('dataset', data)        
        h_mc = h.integrate('dataset', mc)        
        h_mc.scale(-1)
        h_qcd.add(h_mc)

        hist.plot1d(h_qcd, ax=ax, overflow=overflow)
        ax.set_yscale('log')
        ax.set_ylim(1e-2,1e6)
        ax.get_legend().remove()
        
        ax.text(0.,1.,'QCD Estimate in SR',
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

        if distribution in new_xlabels.keys():
            ax.set_xlabel(new_xlabels[distribution])

        outpath = pjoin(outdir, f'qcd_estimation_{distribution}_{year}.pdf')
        fig.savefig(outpath)
        plt.close(fig)
        print(f'File saved: {outpath}')

        # Write the estimate into the output root file
        sumw = h_qcd.values(overflow=overflow)[()]
        xedges = h_qcd.axes()[0].edges(overflow=overflow)
        outrootfile[f'qcd_estimate_{distribution}_{year}'] = (sumw, xedges)

def main():
    inpath = sys.argv[1]
    acc = dir_archive(inpath)
    acc.load('sumw')
    acc.load('sumw2')

    outtag = re.findall('merged_.*', inpath)[0].replace('/','')
    outdir = f'./output/{outtag}/qcd_estimate'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    outrootpath = pjoin(outdir, 'hf_qcd_estimate.root')
    outrootfile = uproot.recreate(outrootpath)
    print(f'ROOT file initiated: {outrootpath}')

    for distribution in distributions:
        get_qcd_estimate(acc, outtag, outrootfile, distribution=distribution)

if __name__ == '__main__':
    main()