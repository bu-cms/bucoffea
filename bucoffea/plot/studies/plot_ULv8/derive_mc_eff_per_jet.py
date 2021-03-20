#!/usr/bin/env python

from bucoffea.plot.studies.plot_ULv8.derive_mc_eff_with_eta import derive_mc_efficiency
import os
import sys
import re
import uproot
import numpy as np
import mplhep as hep

from coffea import hist
from bucoffea.plot.util import merge_datasets, merge_extensions, scale_xs_lumi
from matplotlib import pyplot as plt
from klepto.archives import dir_archive
from pprint import pprint

pjoin = os.path.join

def preprocess(h, acc, region, dataset=re.compile('DYJets.*2017')):
    h = merge_extensions(h, acc, reweight_pu=False)
    scale_xs_lumi(h)
    h = merge_datasets(h)

    h = h.integrate('dataset', dataset).integrate('region', region)
    return h

def derive_efficiency(acc, outtag, region='cr_2m_vbf_relaxed_sel'):
    '''Derive per-jet efficiency of the HF shape cuts as a function of jet eta.'''
    distributions = ['ak4_eta', 'ak4_eta_hf_filtered']
    histos = {}
    for distribution in distributions:
        acc.load(distribution)
        histos[distribution] = preprocess(acc[distribution], acc, region)
    
    data_err_opts = {
        'linestyle':'none',
        'marker': '.',
        'markersize': 10.,
        'color' : 'k'
    }

    fig, ax = plt.subplots()
    hist.plotratio(
        histos['ak4_eta_hf_filtered'],
        histos['ak4_eta'],
        ax=ax,
        error_opts=data_err_opts
    )

    ax.set_ylim(0.7,1.1)
    ax.set_xlabel(r'Jet $\eta$', fontsize=14)
    ax.set_ylabel('Efficiency', fontsize=14)
    
    ax.axhline(1., xmin=0, xmax=1, color='red')
    
    if region == 'cr_2m_vbf_relaxed_sel':
        ax.text(0.,1.,r'$Z(\mu\mu)$',
            fontsize=14,
            ha='left',
            va='bottom',
            transform=ax.transAxes
        )

    outdir = f'./output/{outtag}'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outpath = pjoin(outdir, f'mc_eff_per_jet.pdf')
    fig.savefig(outpath)
    plt.close(fig)

    print(f'File saved: {outpath}')

    # Save the per-jet efficiency to an output ROOT file
    outfile = uproot.recreate(pjoin(outdir, 'mc_eff_per_jet.root'))
    eff = histos['ak4_eta_hf_filtered'].values()[()] / histos['ak4_eta'].values()[()]
    edges = histos['ak4_eta'].axis('jeteta').edges()

    outfile['mc_eff_per_jet'] = (eff, edges)

def main():
    inpath = sys.argv[1]
    acc = dir_archive(inpath)
    acc.load('sumw')
    acc.load('sumw2')

    outtag = re.findall('merged_.*', inpath)[0].replace('/','')

    derive_efficiency(acc, outtag)

if __name__ == '__main__':
    main()