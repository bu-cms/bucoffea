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
    'mjj_nano' : hist.Bin('mjj', r'$M_{jj} \ (GeV)$', [200., 400., 600., 900., 1200., 1500., 2000., 2750., 3500.]),
    'ak4_pt0_nano' : hist.Bin('jetpt',r'Leading AK4 jet $p_{T}$ (GeV)',list(range(80,600,20)) + list(range(600,1000,20)) ),
    'ak4_pt1_nano' : hist.Bin('jetpt',r'Trailing AK4 jet $p_{T}$ (GeV)',list(range(40,600,20)) + list(range(600,1000,20)) ),
    'recoil_pt_nano' : hist.Bin('recoil', r'Recoil (GeV)', [ 250,  280,  310,  340,  370,  400,  430,  470,  510, 550,  590,  640,  690,  740,  790,  840,  900,  960, 1020, 1090, 1160, 1250, 1400])
}

TAGS = {
    'ewk_wmv' : r'EWK $W(\mu\nu)$',
    'ewk_wev' : r'EWK $W(e\nu)$',
    'ewk_zvv' : r'EWK $Z(\nu\nu)$',
    'ewk_zee' : r'EWK $Z(ee)$',
    'ewk_zmm' : r'EWK $Z(\mu\mu)$',
}

RAX_YLIMS = {
    'mjj' : (0.8,1.2)
}

def preprocess(h, acc, distribution='mjj'):
    h = merge_extensions(h, acc)
    scale_xs_lumi(h)
    h = merge_datasets(h)
    
    # Rebin if necessary
    if distribution in BINNINGS.keys():
        new_ax = BINNINGS[distribution]
        h = h.rebin(new_ax.name, new_ax)

    return h

def compare_samples(acc_dict, distribution='mjj', years=[2017,2018]):
    '''Compare samples from private and central Nano production.'''
    histos = {}
    for key, acc in acc_dict.items():
        acc.load(distribution)
        histos[key] = preprocess(acc[distribution], acc, distribution)
    
    data_err_opts = {
        'linestyle':'none',
        'marker': '.',
        'markersize': 10.,
        'color':'k',
    }

    for year in years:
        datasets = [
            ('ewk_wev', re.compile(f'EWKW2Jets.*{year}'), 'cr_1e_vbf'),
            ('ewk_wmv', re.compile(f'EWKW2Jets.*{year}'), 'cr_1m_vbf'),
            ('ewk_zee', re.compile(f'EWKZ2Jets.*ZToLL.*{year}'), 'cr_2e_vbf'),
            ('ewk_zmm', re.compile(f'EWKZ2Jets.*ZToLL.*{year}'), 'cr_2m_vbf'),
            ('ewk_zvv', re.compile(f'EWKZ2Jets.*ZToNuNu.*{year}'), 'sr_vbf_no_veto_all'),
        ]
        for tag, regex, region in datasets:
            if year == 2018 and tag == 'ewk_zvv':
                continue
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

            outdir = './output/11May21_nano'
            if not os.path.exists(outdir):
                os.makedirs(outdir)
            outpath = pjoin(outdir, f'private_vs_central_{tag}_{distribution}_{year}.pdf')
            fig.savefig(outpath)
            plt.close(fig)
            print(f'File saved: {outpath}')

def main():
    acc_dict = {
        'private' : dir_archive( bucoffea_path('submission/merged_2021-05-11_vbfhinv_ULv8_05Feb21_EWK_V_pureNANO') ),
        'central' : dir_archive( bucoffea_path('submission/merged_2021-05-11_vbfhinv_ULv8_05Feb21_central_EWK_V_pureNANO') ),
    }
    
    for acc in acc_dict.values():
        acc.load('sumw')
        acc.load('sumw_pileup')
        acc.load('nevents')

    distributions = [
        'mjj_nano', 
        'ak4_eta0', 
        'ak4_eta1',
        'ak4_pt0_nano',
        'ak4_pt1_nano',
        'recoil_pt_nano',
        'recoil_phi_nano',
    ]

    for distribution in distributions:
        compare_samples(acc_dict, distribution=distribution)

if __name__ == '__main__':
    main()