#!/usr/bin/env python

from bucoffea.plot.studies.plot_ULv8.plot_v7_vs_v8_data import plot_v7_vs_v8_data
import os
import sys
import re
import matplotlib
import numpy as np

from coffea import hist
from bucoffea.plot.util import merge_datasets, merge_extensions, fig_ratio
from bucoffea.plot.style import matplotlib_rc
from matplotlib import pyplot as plt
from klepto.archives import dir_archive
from pprint import pprint

pjoin = os.path.join

matplotlib_rc()
np.seterr(all='ignore')

def preprocess(h, acc, distribution, region='sr_vbf', dataset='MET_2017'):
    h = merge_extensions(h, acc, reweight_pu=False)
    h = merge_datasets(h)

    h = h.integrate('region', region).integrate('dataset', dataset)

    if distribution == 'mjj':
        mjj_ax = hist.Bin('mjj', r'$M_{jj} \ (GeV)$', [200., 400., 600., 900., 1200., 1500., 2000., 2750., 3500., 5000.])
        h = h.rebin('mjj', mjj_ax)

    return h

def plot_cleaning_variable(acc_v7, acc_v8, distribution, region):
    '''Plot comparison of the cleaning variable in the given region.'''
    acc_v7.load(distribution)
    acc_v8.load(distribution)

    h_v7 = preprocess(acc_v7[distribution], acc_v7, distribution, region)
    h_v8 = preprocess(acc_v8[distribution], acc_v8, distribution, region)

    fig, ax, rax = fig_ratio()
    hist.plot1d(h_v7, ax=ax)
    hist.plot1d(h_v8, ax=ax, clear=False)

    ax.set_yscale('log')
    ax.set_ylim(1e-1,1e5)

    ax.legend(title='MET 2017',labels=['ReReco Nano v7', 'UL Nano v8'])

    # Plot v8 / v7 ratio
    data_err_opts = {
        'linestyle':'none',
        'marker': '.',
        'markersize': 10.,
        'color':'k',
    }

    hist.plotratio(h_v8, h_v7, 
            ax=rax,
            unc='num',
            error_opts=data_err_opts
            )

    rax.set_ylabel('v8 / v7')
    rax.grid(True)
    rax.set_ylim(0,2)

    new_xlabels = {
        'ak4_nef0_eeonly': r'Leading Jet EM Fraction',
        'ak4_nef1_eeonly': r'Trailing Jet EM Fraction',
        'dphitkpf': r'$\Delta\phi$(TK,PF)',
    }

    if distribution in new_xlabels.keys():
        ax.set_xlabel(new_xlabels[distribution])
        rax.set_xlabel(new_xlabels[distribution])

    if 'eeonly' in distribution:
        ltext = 'EE Jets Only'
        rtext = 'VBF SR (No EM frac veto)'
    elif distribution == 'dPFTkMET':
        ltext = r'Jets: $2.8<|\eta|<3.2$'
        rtext = 'VBF SR (No horn veto)'
        ax.set_xlim(0,1)
    elif region == 'sr_vbf_no_eemitigation':
        ltext = ''
        rtext = 'VBF SR (No VecB/VecDPhi cut)'
    
    ax.text(0.,1., ltext,
        fontsize=14,
        ha='left',
        va='bottom',
        transform=ax.transAxes
    )
    ax.text(1.,1., rtext,
        fontsize=14,
        ha='right',
        va='bottom',
        transform=ax.transAxes
    )

    # Save figure
    outdir = './output/v7_vs_v8/03Mar21/cleaning_variables'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    outpath = pjoin(outdir, f'met_2017_v7_v8_{distribution}.pdf')
    fig.savefig(outpath)
    plt.close(fig)
    print(f'File saved: {outpath}')

def main():
    inpath_v7 = '../../../submission/merged_2021-03-04_vbfhinv_03Sep20v7_MET_2017'
    inpath_v8 = '../../../submission/merged_2021-03-04_vbfhinv_ULv8_MET_2017'
    acc_v7 = dir_archive(inpath_v7) 
    acc_v8 = dir_archive(inpath_v8) 

    for acc in [acc_v7, acc_v8]:
        acc.load('sumw')
        acc.load('sumw2')

    distributions_regions = [
        ('ak4_nef0_eeonly','sr_vbf_no_emfraccut'),
        ('ak4_nef1_eeonly','sr_vbf_no_emfraccut'),
        ('dPFTkMET','sr_vbf_nohornveto'),
        ('vecb','sr_vbf_no_eemitigation'),
        ('dphitkpf','sr_vbf_no_eemitigation'),
    ]

    for distribution, region in distributions_regions:
        plot_cleaning_variable(acc_v7, acc_v8, distribution, region)

if __name__ == '__main__':
    main()