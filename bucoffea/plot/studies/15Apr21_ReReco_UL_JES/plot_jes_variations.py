#!/usr/bin/env python
import os
import re
import sys
import uproot
import numpy as np
import mplhep as hep

from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
from coffea import hist
from bucoffea.plot.util import merge_datasets, merge_extensions, scale_xs_lumi, fig_ratio
from klepto.archives import dir_archive
from pprint import pprint

pjoin = os.path.join

def plot_jes_variations(acc, outtag, rootfile, distribution='mjj', year=2017):
    '''Plot JES/JER variations of the total MC in VBF SR.'''
    acc.load(distribution)
    h = acc[distribution]

    h = merge_extensions(h, acc, reweight_pu=False)
    scale_xs_lumi(h)
    h = merge_datasets(h)
    
    if distribution == 'mjj':
        new_ax = hist.Bin('mjj', r'$M_{jj} \ (GeV)$', [200., 400., 600., 900., 1200., 1500., 2000., 2750., 3500., 5000.])
        h = h.rebin('mjj', new_ax)

    mc = re.compile(f'(ZJetsToNuNu.*|EW.*|Top_FXFX.*|Diboson.*|DYJetsToLL_M-50_HT_MLM.*|WJetsToLNu.*HT.*).*{year}')
    h = h.integrate('dataset', mc)

    regions = [
        'sr_vbf_jerUp',
        'sr_vbf_jerDown',
        'sr_vbf_jesTotalUp',
        'sr_vbf_jesTotalDown',
    ]

    fig, ax, rax = fig_ratio()
    hist.plot1d(h, ax=ax, overlay='region')

    ax.set_yscale('log')
    ax.set_ylim(1e0,1e6)

    handles, labels = ax.get_legend_handles_labels()
    new_legend_labels = {
        'sr_vbf': 'Nominal',
        'sr_vbf_jerUp': 'JER up',
        'sr_vbf_jerDown': 'JER down',
        'sr_vbf_jesTotalUp': 'JES up',
        'sr_vbf_jesTotalDown': 'JES down',
    }
    for handle, label in zip(handles, labels):
        handle.set_label(new_legend_labels[label])

    ax.legend(handles=handles)

    ax.text(0.,1.,f'{year} Total MC Bkg',
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

    ax.yaxis.set_ticks_position('both')

    data_err_opts = {
        'linestyle':'none',
        'marker': '.',
        'markersize': 10.,
    }

    h_nom = h.integrate('region', 'sr_vbf')
    for region in regions:
        h_var = h.integrate('region', region)
        hist.plotratio(
            h_var, 
            h_nom,
            ax=rax,
            unc='num',
            error_opts=data_err_opts,
            label=region.split('_')[-1],
            clear=False
        )

    rax.axhline(1.,xmin=0,xmax=1,color='k')

    rax.set_ylim(0.5,1.5)
    rax.set_ylabel('Var / Nom')
    # rax.legend(ncol=2)

    loc = MultipleLocator(0.25)
    rax.yaxis.set_minor_locator(loc)
    rax.grid(True, which='both')

    outdir = f'./output/{outtag}'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outpath = pjoin(outdir, f'mc_JES_JER_{distribution}_{year}.pdf')
    fig.savefig(outpath)
    plt.close(fig)

    print(f'File saved: {outpath}')

    # Save the individual variations into the given ROOT file
    xedges = h_nom.axis('mjj').edges()
    for region in regions:
        tag = region.split('_')[-1]
        ratio = h.integrate('region', region).values()[()] / h_nom.values()[()]
        rootfile[f'MTR_{year}_{tag}'] = (ratio, xedges)

def main():
    inpath = sys.argv[1]
    acc = dir_archive(inpath)
    acc.load('sumw')
    acc.load('sumw2')

    outtag = re.findall('merged_.*', inpath)[0].replace('/','')
    
    outdir = f'./output/{outtag}'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    rootpath = pjoin(outdir, 'jes_uncs.root')
    rootfile = uproot.recreate(rootpath)

    for year in [2017, 2018]:
        plot_jes_variations(acc, outtag, rootfile=rootfile, year=year)

if __name__ == '__main__':
    main()