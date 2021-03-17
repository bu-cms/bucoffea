#!/usr/bin/env python

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

def derive_mc_efficiency(acc, outtag, distribution='ak4_eta0', regionbase='cr_2m_vbf_relaxed_sel'):
    '''Derive and save the efficiency of the HF cuts on MC, as a function of jet eta.'''
    acc.load(distribution)
    h = acc[distribution]
    
    h = merge_extensions(h, acc, reweight_pu=False)
    scale_xs_lumi(h)
    h = merge_datasets(h)

    h = h.integrate('dataset', re.compile('DYJets.*2017'))

    data_err_opts = {
        'linestyle':'none',
        'marker': '.',
        'markersize': 10.,
    }

    fig, ax = plt.subplots()
    hist.plotratio(
        h.integrate('region', f'{regionbase}_with_hfcuts'),
        h.integrate('region', f'{regionbase}'),
        ax=ax,
        error_opts=data_err_opts,
    )

    ax.set_ylim(0.6,1.1)
    ax.set_ylabel('Efficiency')
    ax.set_xlabel(r'Leading Jet $\eta$')

    ax.axhline(1., xmin=0, xmax=1, color='k')

    ax.text(0.,1.,'DY ULv8 MC 2017',
        fontsize=14,
        ha='left',
        va='bottom',
        transform=ax.transAxes
    )

    # Save figure
    outdir = f'./output/{outtag}/eff'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outpath = pjoin(outdir, f'mc_eff_{distribution}.pdf')
    fig.savefig(outpath)
    plt.close(fig)

    print(f'File saved: {outpath}')

    outfile = uproot.recreate(pjoin(outdir, 'mc_eff.root'))
    eff = h.integrate('region', f'{regionbase}_with_hfcuts').values()[()] / h.integrate('region', f'{regionbase}').values()[()]
    eff[np.isnan(eff) | np.isinf(eff)] =1.
    edges = h.axis('jeteta').edges()
    outfile['eff_with_leading_jet_eta'] = (eff, edges)

def main():
    inpath = sys.argv[1]
    acc = dir_archive(inpath)
    acc.load('sumw')
    acc.load('sumw2')
    
    outtag = re.findall('merged_.*', inpath)[0].replace('/','')

    derive_mc_efficiency(acc, outtag)

if __name__ == '__main__':
    main()