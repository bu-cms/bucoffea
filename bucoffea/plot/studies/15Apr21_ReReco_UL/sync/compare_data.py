#!/usr/bin/env python

import os
import sys
import re
import uproot
import mplhep as hep
import numpy as np

from matplotlib import pyplot as plt
from coffea import hist
from bucoffea.plot.util import merge_datasets, merge_extensions, scale_xs_lumi, fig_ratio
from klepto.archives import dir_archive
from pprint import pprint

pjoin = os.path.join

def dump_to_root(acc, dataset, region, outputrootfile, distribution='mjj'):
    '''Given the coffea accumulator, prepare a ROOT file with the data templates from BU.'''
    acc.load(distribution)
    h = acc[distribution]

    h = merge_extensions(h, acc)
    scale_xs_lumi(h)
    h = merge_datasets(h)

    if distribution == 'mjj':
        mjj_bins = [200., 400., 600., 900., 1200., 1500., 2000., 2750., 3500., 5000.]
        mjj_ax = hist.Bin('mjj', r'$M_{jj} \ (GeV)$', mjj_bins)
        h = h.rebin('mjj', mjj_ax)

    h = h.integrate('region', region).integrate('dataset', dataset)

    # Save to output ROOT file
    outputrootfile[f'{region}_mjj'] = hist.export1d(h)

def compare_data(bu_file, ic_file, year):
    '''Compare data templates between BU and IC.'''
    region_name_mapping = {
        # BU region name: IC name
        'sr_vbf'    : f'MTR_{year}_SR',
        'cr_1m_vbf' : f'MTR_{year}_WMUNU',
        'cr_1e_vbf' : f'MTR_{year}_WENU',
        'cr_2m_vbf' : f'MTR_{year}_ZMUMU',
        'cr_2e_vbf' : f'MTR_{year}_ZEE',
    }

    bu_regions = region_name_mapping.keys()
    ic_regions = region_name_mapping.values()

    pretty_labels = {
        'sr_vbf' : f'VBF SR {year}',        
        'cr_1m_vbf' : f'VBF $W(\\mu\\nu)$ {year}',        
        'cr_1e_vbf' : f'VBF $W(e\\nu)$ {year}',        
        'cr_2m_vbf' : f'VBF $Z(\\mu\\mu)$ {year}',
        'cr_2e_vbf' : f'VBF $Z(ee)$ {year}',
    }

    f_ic = ic_file['shapes_prefit']
    for bu_region, ic_region in zip(bu_regions, ic_regions):
        h_bu = bu_file[f'{bu_region}_mjj']
        h_ic = f_ic[ic_region]['data']

        xedges = h_bu.edges

        fig, ax, rax = fig_ratio()
        hep.histplot(h_bu.values, xedges, ax=ax, binwnorm=1, label='BU')
        hep.histplot(h_ic.yvalues, xedges, ax=ax, label='IC')

        ax.set_yscale('log')
        ax.set_ylim(1e-4,1e6)
        ax.legend()

        ax.set_xlim(0,5000.)
        ax.set_xlabel(r'$M_{jj} \ (GeV)$')
        ax.yaxis.set_ticks_position('both')

        ax.text(0.,1., pretty_labels[bu_region],
            fontsize=14,
            ha='left',
            va='bottom',
            transform=ax.transAxes
        )

        data_err_opts = {
            'linestyle':'none',
            'marker': '.',
            'markersize': 10.,
            'color':'k',
        }

        binw = np.diff(xedges)
        ratio = (h_bu.values / binw) / h_ic.yvalues

        xcenters = 0.5 * (xedges + np.roll(xedges, -1))[:-1]
        rax.plot(xcenters, ratio, **data_err_opts)

        rax.grid(True)
        rax.set_ylim(0.8,1.2)
        rax.set_ylabel('BU / IC')

        outdir = './output/data_comparison'
        try:
            os.makedirs(outdir)
        except FileExistsError:
            pass

        outpath = pjoin(outdir, f'{bu_region}_data_comp_{year}.pdf')
        fig.savefig(outpath)
        plt.close(fig)

        print(f'File saved: {outpath}')

def main():
    inpath = sys.argv[1]
    acc = dir_archive(inpath)

    outtag = re.findall('merged.*', inpath)[0].replace('/','')

    outdir = f'./output/{outtag}'
    try:
        os.makedirs(outdir)

    except FileExistsError:
        pass

    for year in [2017, 2018]:
        # Fit diag file of IC, to get data templates to compare with BU
        ic_fitdiag_file = uproot.open(pjoin(f'input/07Jun21_sync/IC_fitDiagnosticsMTR_{year}.root'))
        bu_root_path = pjoin(outdir, f'BU_data_templates_{year}.root')
        if not os.path.exists(bu_root_path):
            # Convert BU inputs to coffea histograms (if the file doesn't exist already)
            bu_root_file = uproot.open(bu_root_path)
            datasets_regions = [
                (f'MET_{year}', 'sr_vbf'),
                (f'MET_{year}', 'cr_1m_vbf'),
                (f'MET_{year}', 'cr_2m_vbf'),
                (f'EGamma_{year}', 'cr_1e_vbf'),
                (f'EGamma_{year}', 'cr_2e_vbf'),
            ]
    
            for dataset, region in datasets_regions:
                dump_to_root(acc, dataset=dataset, region=region, outputrootfile=bu_root_file)

        compare_data(bu_file=uproot.open(bu_root_path), ic_file=ic_fitdiag_file, year=year)

if __name__ == '__main__':
    main()