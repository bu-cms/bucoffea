#!/usr/bin/env python
import os
import re
import sys
import numpy as np

from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
from coffea import hist
from bucoffea.plot.util import merge_datasets, merge_extensions, scale_xs_lumi, fig_ratio
from klepto.archives import dir_archive
from pprint import pprint

pjoin = os.path.join

def compare_nlo_eoy_ul(acc, outtag, region='cr_2m_vbf', distribution='mjj'):
    '''Compare NLO samples between UL and EOY versions.'''
    acc.load(distribution)
    h = acc[distribution]
    
    h = merge_extensions(h, acc)
    scale_xs_lumi(h)
    h = merge_datasets(h)
    
    if distribution == 'mjj':
        mjj_ax = hist.Bin('mjj', r'$M_{jj} \ (GeV)$', [200., 400., 600., 900., 1200., 1500.,2000., 2750., 3500., 5000.])
        h = h.rebin('mjj', mjj_ax)

    pretty_region_tag = {
        'sr_vbf_no_veto_all' : 'VBF SR',       
        'cr_2m_vbf' : r'$Z(\mu\mu)$ CR',       
        'cr_2e_vbf' : r'$Z(ee)$ CR',       
    }

    h = h.integrate('region', region)

    for year in [2017, 2018]:

        fig, ax, rax = fig_ratio()
        hist.plot1d(h[re.compile(f'.*{year}')], ax=ax, overlay='dataset')

        ax.set_yscale('log')
        ax.set_ylim(1e-2,1e6)

        handles, labels = ax.get_legend_handles_labels()
        for handle, label in zip(handles, labels):
            if label.startswith('DYJetsToLL_Pt'):
                handle.set_label('UL')
            else:
                handle.set_label('EOY')

        ax.legend(title='NLO Version', handles=handles)

        ax.text(0.,1.,pretty_region_tag[region],
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

        data_err_opts = {
            'linestyle':'none',
            'marker': '.',
            'markersize': 10.,
            'color':'k',
        }

        hist.plotratio(
            h.integrate('dataset', f'DYJetsToLL_Pt_FXFX_{year}'),
            h.integrate('dataset', f'DYNJetsToLL_M-50_LHEZpT-FXFX_{year}'),
            ax=rax,
            unc='num',
            error_opts=data_err_opts
        )

        rax.grid(True)
        rax.set_ylim(0.5,1.5)
        rax.set_ylabel('UL / EOY')

        outdir = f'./output/{outtag}/nlo_vs_nlo'
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        outpath = pjoin(outdir, f'nlo_eoy_ul_comp_{region}_{year}.pdf')
        fig.savefig(outpath)
        plt.close(fig)
    
        print(f'File saved: {outpath}')

def main():
    inpath = sys.argv[1]
    acc = dir_archive(inpath)
    acc.load('sumw')
    acc.load('sumw_pileup')
    acc.load('nevents')

    outtag = re.findall('merged_.*', inpath)[0].replace('/','')

    for region in ['cr_2m_vbf', 'cr_2e_vbf']:
        compare_nlo_eoy_ul(acc, 
            outtag=outtag,
            region=region,
            )

if  __name__ == '__main__':
    main()