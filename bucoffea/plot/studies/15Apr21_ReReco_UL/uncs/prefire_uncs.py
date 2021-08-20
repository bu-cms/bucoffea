#!/usr/bin/env python

import os
import sys
import re
import uproot
import numpy as np

from coffea import hist
from bucoffea.plot.util import merge_datasets, merge_extensions, scale_xs_lumi, fig_ratio
from matplotlib import pyplot as plt
from klepto.archives import dir_archive
from pprint import pprint

pjoin = os.path.join

def plot_prefire_uncs(acc, outtag, dataset, year, region='sr_vbf_no_veto_all', outrootfile=None):
    distribution = 'mjj_pref_weights'
    acc.load(distribution)
    h = acc[distribution]

    h = merge_extensions(h, acc)
    scale_xs_lumi(h)
    h = merge_datasets(h)

    mjj_bins = [200., 400., 600., 900., 1200., 1500., 2000., 2750., 3500., 5000.]
    mjj_ax = hist.Bin('mjj', r'$M_{jj} \ (GeV)$', mjj_bins)    
    h = h.rebin('mjj', mjj_ax)

    h = h.integrate('region', region).integrate('dataset', re.compile(f'{dataset}.*{year}'))

    fig, ax, rax = fig_ratio()
    hist.plot1d(h, ax=ax, overlay='variation')

    ax.set_yscale('log')
    ax.set_ylim(1e-2,1e6)

    data_err_opts = {
        'linestyle':'none',
        'marker': '.',
        'markersize': 10.,
    }

    h_nom = h.integrate('variation', 'nom')
    for variation in ['up', 'down']:
        hist.plotratio(
            h.integrate('variation', variation),
            h_nom,
            ax=rax,
            unc='num',
            error_opts=data_err_opts,
            label=f'Prefire {variation}',
            clear=False
        )

    rax.set_ylabel('Variation / Nom')
    rax.set_ylim(0.8,1.2)
    rax.legend(ncol=2)
    rax.grid(True)

    outdir = f'./output/{outtag}'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    outpath = pjoin(outdir, f'{dataset}_prefire_uncs_{year}.pdf')
    fig.savefig(outpath)
    plt.close(fig)

    print(f'File saved: {outpath}')

    r_up = h.integrate('variation', 'up').values()[()] / h_nom.values()[()]
    r_down = h.integrate('variation', 'down').values()[()] / h_nom.values()[()]
    outrootfile[f'{dataset}_{year}_CMS_L1prefire_2017_up'] = (r_up, h_nom.axis('mjj').edges())
    outrootfile[f'{dataset}_{year}_CMS_L1prefire_2017_down'] = (r_down, h_nom.axis('mjj').edges())

def main():
    inpath = sys.argv[1]
    acc = dir_archive(inpath)
    acc.load('sumw')
    acc.load('sumw_pileup')
    acc.load('nevents')

    outtag = re.findall('merged_.*', inpath)[0].replace('/','')

    outdir = f'./output/{outtag}'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outrootpath = pjoin(outdir, 'vbf_prefire_uncs.root')
    outrootfile = uproot.recreate(outrootpath)

    for year in [2017]:
        plot_prefire_uncs(acc, outtag, 
            dataset='VBF_HToInvisible', 
            year=year,
            outrootfile=outrootfile
            )

if __name__ == '__main__':
    main()