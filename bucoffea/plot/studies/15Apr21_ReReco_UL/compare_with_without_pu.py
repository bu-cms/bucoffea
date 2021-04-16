#!/usr/bin/env python

import os
import re
import sys
import numpy as np
import mplhep as hep

from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
from coffea import hist
from bucoffea.plot.util import merge_datasets, merge_extensions, scale_xs_lumi, fig_ratio
from bucoffea.plot.style import matplotlib_rc
from klepto.archives import dir_archive
from pprint import pprint
from distributions import distributions, binnings

pjoin = os.path.join

matplotlib_rc()

def compare_with_without_pu(acc, outtag, dataset_tag, year, distribution='mjj', fformat='pdf'):
    acc.load(distribution)
    h = acc[distribution]
    
    if distribution == 'mjj':
        overflow = 'over'
    else:
        overflow = 'none'

    h = merge_extensions(h, acc, reweight_pu=False)
    scale_xs_lumi(h)
    h = merge_datasets(h)

    dataset_regex_dict = {
        'EWK_zll' : re.compile(f'EWKZ2Jets.*ZToLL.*{year}'),
        'EWK_zvv' : re.compile(f'EWKZ2Jets.*ZToNuNu.*{year}'),
    }
    
    h = h.integrate('dataset', dataset_regex_dict[dataset_tag])[re.compile('sr_vbf(_no_pu)?$')]

    if distribution in binnings.keys():
        new_ax = binnings[distribution]
        h = h.rebin(new_ax.name, new_ax)

    fig, ax, rax = fig_ratio()
    hist.plot1d(h, ax=ax, overlay='region', binwnorm=1)

    ax.set_yscale('log')
    ax.set_ylim(1e-2,1e4)
    ax.set_ylabel('Events / GeV')

    handles, labels = ax.get_legend_handles_labels()
    legend_labels = {
        'sr_vbf' : 'With PU Weights',
        'sr_vbf_no_pu' : 'Without PU Weights',
    }
    for handle, label in zip(handles, labels):
        handle.set_label(legend_labels[label])

    ax.legend(handles=handles)

    data_err_opts = {
        'linestyle':'none',
        'marker': '.',
        'markersize': 10.,
        'color':'k',
    }

    hist.plotratio(
        h.integrate('region', 'sr_vbf_no_pu'),
        h.integrate('region', 'sr_vbf'),
        ax=rax,
        unc='num',
        error_opts=data_err_opts
    )

    rax.grid(True)
    rax.set_ylim(0.95,1.05)
    rax.set_ylabel('Without PU / With')

    outdir = f'./output/{outtag}/pu_comparison'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outpath = pjoin(outdir, f'{dataset_tag}_{distribution}_{year}.{fformat}')
    fig.savefig(outpath)
    plt.close(fig)

    print(f'File saved: {outpath}')

def main():
    inpath = sys.argv[1]
    acc = dir_archive(inpath)
    acc.load('sumw')
    acc.load('sumw2')

    outtag = re.findall('merged_.*', inpath)[0].replace('/','')

    dataset_tag = 'EWK_zvv'
    for year in [2017, 2018]:
        for fformat in ['pdf', 'png']:
            compare_with_without_pu(acc, outtag, 
                dataset_tag=dataset_tag,
                year=year,
                fformat=fformat
            )

if __name__ == '__main__':
    main()