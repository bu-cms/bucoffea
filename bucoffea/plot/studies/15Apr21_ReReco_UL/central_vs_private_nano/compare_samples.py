#!/usr/bin/env python
import os
import re
import sys
import warnings
import numpy as np
import mplhep as hep

from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
from coffea import hist
from scipy.stats import distributions
from bucoffea.plot.util import merge_datasets, merge_extensions, scale_xs_lumi, fig_ratio
from bucoffea.helpers.paths import bucoffea_path
from klepto.archives import dir_archive
from pprint import pprint

pjoin = os.path.join

warnings.filterwarnings('ignore')

BINNINGS = {
    'mjj' : hist.Bin('mjj', r'$M_{jj} \ (GeV)$', [200., 400., 600., 900., 1200., 1500., 2000., 2750., 3500.]),
    'ak4_pt0' : hist.Bin('jetpt',r'Leading AK4 jet $p_{T}$ (GeV)',list(range(80,600,20)) + list(range(600,1000,20)) ),
    'ak4_pt1' : hist.Bin('jetpt',r'Trailing AK4 jet $p_{T}$ (GeV)',list(range(40,600,20)) + list(range(600,1000,20)) ),
}

TAGS = {
    'ewk_wlv' : r'EWK $W(\mu\nu)$',
    'ewk_zvv' : r'EWK $Z(\nu\nu)$',
}

RAX_YLIMS = {
    'mjj' : (0.8,1.2)
}

def preprocess(h, acc, region, distribution='mjj'):
    h = merge_extensions(h, acc)
    scale_xs_lumi(h)
    h = merge_datasets(h)
    
    # Rebin if necessary
    if distribution in BINNINGS.keys():
        new_ax = BINNINGS[distribution]
        h = h.rebin(new_ax.name, new_ax)

    return h

def compare_samples(acc_dict, region, distribution='mjj', years=[2017]):
    '''Compare samples from private and central Nano production.'''
    histos = {}
    for key, acc in acc_dict.items():
        acc.load(distribution)
        histos[key] = preprocess(acc[distribution], acc, region, distribution)
    
    data_err_opts = {
        'linestyle':'none',
        'marker': '.',
        'markersize': 10.,
        'color':'k',
    }

    for year in years:
        datasets = [
            ('ewk_wlv', re.compile(f'EWKW2Jets.*{year}'), 'cr_1m_vbf'),
            ('ewk_zvv', re.compile(f'EWKZ2Jets.*ZToNuNu.*{year}'), 'sr_vbf_no_veto_all'),
        ]
        for tag, regex, region in datasets:
            h_private = histos['private'].integrate('dataset', regex).integrate('region', region)        
            h_central = histos['central'].integrate('dataset', regex).integrate('region', region)

            fig, ax, rax = fig_ratio()
            hist.plot1d(h_private, ax=ax, overflow='over')
            hist.plot1d(h_central, ax=ax, overflow='over', clear=False)

            ax.set_yscale('log')
            ax.set_ylim(1e-2,1e4)
            ax.legend(title='Nano', labels=['Private', 'Central'])

            ax.text(0.,1.,TAGS[tag],
                fontsize=14,
                ha='left',
                va='bottom',
                transform=ax.transAxes
            )

            # Plot ratio of the two versions
            hist.plotratio(
                h_private,
                h_central,
                ax=rax,
                unc='num',
                overflow='over',
                error_opts=data_err_opts
            )

            new_xlabels = {
                'ak4_eta0' : r'Leading Jet $\eta$',
                'ak4_eta1' : r'Trailing Jet $\eta$',
            }

            if distribution in new_xlabels.keys():
                ax.set_xlabel(new_xlabels[distribution])
                rax.set_xlabel(new_xlabels[distribution])

            rax.grid(True)
            if distribution in RAX_YLIMS.keys():
                rax.set_ylim(*RAX_YLIMS[distribution])
            else:
                rax.set_ylim(0.5,1.5)
            rax.set_ylabel('Private / Central')

            outdir = './output'
            if not os.path.exists(outdir):
                os.makedirs(outdir)
            outpath = pjoin(outdir, f'private_vs_central_{tag}_{distribution}_{year}.pdf')
            fig.savefig(outpath)
            plt.close(fig)
            print(f'File saved: {outpath}')

def main():
    acc_dict = {
        'private' : dir_archive( bucoffea_path('submission/merged_2021-05-03_vbfhinv_ULv8_05Feb21_trigSF_update') ),
        'central' : dir_archive( bucoffea_path('submission/merged_2021-05-07_vbfhinv_ULv8_05Feb21_central_EWK_V') ),
    }
    
    for acc in acc_dict.values():
        acc.load('sumw')
        acc.load('sumw_pileup')
        acc.load('nevents')

    distributions = [
        'mjj', 
        'ak4_eta0', 
        'ak4_eta1',
        'ak4_pt0',
        'ak4_pt1',
    ]

    for distribution in distributions:
        compare_samples(acc_dict, region='sr_vbf_no_veto_all', distribution=distribution)

if __name__ == '__main__':
    main()