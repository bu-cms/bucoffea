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

pjoin = os.path.join

def get_qcd_estimate(acc, outtag, outrootfile, distribution):
    '''Calculate the QCD template in SR'''
    acc.load(distribution)
    h = acc[distribution]

    h = merge_extensions(h, acc, reweight_pu=False)
    scale_xs_lumi(h)
    h = merge_datasets(h)

    if distribution == 'mjj':
        mjj_ax = hist.Bin('mjj', r'$M_{jj} \ (GeV)$', [200., 400., 600., 900., 1200., 1500., 2000., 2750., 3500.])
        h = h.rebin('mjj', mjj_ax)
    
    outdir = f'./output/{outtag}/qcd_estimate'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # Get data and MC yields in the QCD CR
    h = h.integrate('region', 'cr_vbf_qcd')
    for year in [2017, 2018]:
        data = f'MET_{year}'
        mc = re.compile(f'(ZJetsToNuNu.*|EW.*|Top_FXFX.*|Diboson.*|DYJetsToLL_M-50_HT_MLM.*|WJetsToLNu.*HT.*).*{year}')

        fig, ax = plt.subplots()
        hist.plot1d(h.integrate('dataset', data), ax=ax, overflow='over')
        hist.plot1d(h[mc], ax=ax, clear=False, overlay='dataset', stack=True, overflow='over')

        ax.set_yscale('log')
        ax.set_ylim(1e-4,1e4)
        
        ax.text(0.,1.,'QCD CR',
            fontsize=14,
            ha='left',
            va='bottom',
            transform=ax.transAxes
        )

        handles, labels = ax.get_legend_handles_labels()
        for handle, label in zip(handles, labels):
            if label == 'None':
                handle.set_label('Data')

        ax.legend(handles=handles)

        outpath = pjoin(outdir, f'qcd_cr_{year}.pdf')
        fig.savefig(outpath)
        plt.close(fig)
        print(f'File saved: {outpath}')

        fig, ax = plt.subplots()
        # Calculate the QCD template in SR: Data - MC in CR (already weighted by TF)
        h_qcd = h.integrate('dataset', data)        
        h_mc = h.integrate('dataset', mc)        
        h_mc.scale(-1)
        h_qcd.add(h_mc)

        hist.plot1d(h_qcd, ax=ax, overflow='over')
        ax.set_yscale('log')
        ax.set_ylim(1e-4,1e2)
        ax.get_legend().remove()
        
        ax.text(0.,1.,'QCD Estimate in SR',
            fontsize=14,
            ha='left',
            va='bottom',
            transform=ax.transAxes
        )

        outpath = pjoin(outdir, f'qcd_estimation_{year}.pdf')
        fig.savefig(outpath)
        plt.close(fig)
        print(f'File saved: {outpath}')

        # Write the estimate into the output root file
        sumw = h_qcd.values(overflow='over')[()]
        xedges = h_qcd.axes()[0].edges(overflow='over')
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

    for distribution in ['mjj','ak4_eta0','ak4_eta1']:
        get_qcd_estimate(acc, outtag, outrootfile, distribution=distribution)

if __name__ == '__main__':
    main()