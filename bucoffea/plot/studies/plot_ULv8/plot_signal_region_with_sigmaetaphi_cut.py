#!/usr/bin/env python

import os
import sys
import re
import numpy as np
import matplotlib.colors as colors

from coffea import hist
from numpy.core.numeric import extend_all
from scipy.stats import distributions
from scipy.stats.stats import _validate_distribution
from bucoffea.plot.util import fig_ratio, merge_datasets, merge_extensions, scale_xs_lumi
from matplotlib import pyplot as plt
from klepto.archives import dir_archive
from pprint import pprint

pjoin = os.path.join

def plot_signal_region(acc, outtag, dataset='MET_2017', distribution='mjj'):
    acc.load(distribution)
    h = acc[distribution]

    h = merge_extensions(h, acc, reweight_pu=False)
    scale_xs_lumi(h)
    h = merge_datasets(h)

    h = h.integrate('dataset', dataset)

    # Rebin mjj
    if distribution == 'mjj':
        mjj_bins = [200., 400., 600., 900., 1200., 1500., 2000., 2750., 3500., 5000.]
        new_ax = hist.Bin('mjj', r'$M_{jj} \ (GeV)$', mjj_bins)
        h = h.rebin('mjj', new_ax)
    
    fig, ax, rax = fig_ratio()
    hist.plot1d(h, ax=ax, overlay='region') 

    ax.set_yscale('log')
    ax.set_ylim(1e-2,1e6)

    ax.text(0., 1., 'MET 2017 Nanov8',
        fontsize=14,
        ha='left',
        va='bottom',
        transform=ax.transAxes
    )
    newlabels = {
        'sr_vbf_with_hfcut' : r'With $\sigma_{\phi\phi} / \sigma_{\eta\eta}$ cut',
        'sr_vbf_no_hfhfveto' : r'HF-HF events included',
        'sr_vbf' : r'With HF-HF veto',
    }

    handles, labels = ax.get_legend_handles_labels()
    for handle, label in zip(handles, labels):
        handle.set_label(newlabels[label])
    
    ax.legend(title='VBF Signal Region', handles=handles)

    data_err_opts = {
        'linestyle':'none',
        'marker': '.',
        'markersize': 10.,
    }

    h_with_hfhf = h.integrate('region', 'sr_vbf_no_hfhfveto')
    colors = ['blue', 'green']
    regions = ['sr_vbf_with_hfcut', 'sr_vbf']
    for region, color in zip(regions, colors):
        _data_err_opts = data_err_opts
        _data_err_opts['color'] = color
        hist.plotratio(
            h.integrate('region', region),
            h_with_hfhf,
            ax=rax,
            unc='num',
            label=newlabels[region],
            error_opts=_data_err_opts,
            clear=False
        )

    rax.legend(ncol=2)
    rax.grid(True)
    if 'ak4_pt' in distribution:
        rax.set_ylim(0.8,1.2)
    else:
        rax.set_ylim(0.4,1.6)
    rax.set_ylabel('Ratio to HF-HF Included')

    # Save figure
    outdir = f'./output/{outtag}'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outpath = pjoin(outdir, f'{distribution}_signal_region.pdf')
    fig.savefig(outpath)
    plt.close(fig)

    print(f'File saved: {outpath}')

def main():
    inpath = sys.argv[1]
    acc = dir_archive(inpath)
    acc.load('sumw')
    acc.load('sumw2')

    outtag = re.findall('merged_.*', inpath)[0].replace('/','')

    distributions = ['mjj', 'ak4_eta0', 'ak4_eta1', 'ak4_pt0', 'ak4_pt1']

    for distribution in distributions:
        plot_signal_region(acc, outtag, distribution=distribution)

if __name__ == '__main__':
    main()
