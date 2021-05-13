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

def compare_lo_vs_nlo(acc, outtag, tag, lo_regex, nlo_regex, distribution='mjj', region='sr_vbf_no_veto_all'):
    '''Compare MC templates with LO+SF and NLO.'''
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

    for year in [2017, 2018]:
        _h = h.integrate('region', region)[re.compile(f'(({lo_regex})|({nlo_regex})).*{year}')]
        
        fig, ax, rax = fig_ratio()
    
        hist.plot1d(_h, ax=ax, overlay='dataset')

        ax.set_yscale('log')
        ax.set_ylim(1e-2,1e6)
        ax.yaxis.set_ticks_position('both')

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

        handles, labels = ax.get_legend_handles_labels()
        for handle, label in zip(handles, labels):
            if re.match('DYJets.*HT.*', label):
                handle.set_label('LO+SF')
            elif re.match('DYJets.*Pt.*FXFX.*', label):
                handle.set_label('NLO')

        ax.legend(title='DY Sample', handles=handles)

        data_err_opts = {
            'linestyle':'none',
            'marker': '.',
            'markersize': 10.,
            'color':'k',
        }
    
        hist.plotratio(
            _h.integrate('dataset', re.compile(f'{nlo_regex}.*{year}')),
            _h.integrate('dataset', re.compile(f'{lo_regex}.*{year}')),
            ax=rax,
            unc='num',
            error_opts=data_err_opts
        )
    
        rax.grid(True)
        rax.set_ylim(0.5,1.5)
        rax.set_ylabel('NLO / LO+SF')
    
        outdir = f'./output/{outtag}/lo_vs_nlo'
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        outpath = pjoin(outdir, f'{tag}_{region}_{year}.pdf')
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

    for region in ['sr_vbf_no_veto_all', 'cr_2m_vbf', 'cr_2e_vbf']:
        compare_lo_vs_nlo(acc, outtag,
            tag='dy',
            lo_regex='DYJetsToLL.*HT.*',
            nlo_regex='DYJets.*Pt.*',
            region=region,
        )

if  __name__ == '__main__':
    main()