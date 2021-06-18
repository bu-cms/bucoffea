#!/usr/bin/env python

import os
import re
import sys
import uproot
import numpy as np
import mplhep as hep

from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
from coffea import hist
from bucoffea.plot.util import merge_datasets, merge_extensions, scale_xs_lumi, fig_ratio, lumi
from bucoffea.helpers.paths import bucoffea_path
from klepto.archives import dir_archive
from pprint import pprint

pjoin = os.path.join

# Suppress true_divide warnings
np.seterr(divide='ignore', invalid='ignore')

def make_plot(acc, outtag, distribution='npv_nopu', region='cr_1m_vbf'):
    acc.load(distribution)
    h = acc[distribution]

    h = merge_extensions(h, acc)
    scale_xs_lumi(h)
    h = merge_datasets(h)

    h = h.integrate('region', region)

    for year in [2017, 2018]:
        data = f'MET_{year}'
        mc = re.compile(f'(EWKW.*|Top_FXFX.*|Diboson.*|DYJetsToLL_M-50_HT_MLM.*|WJetsToLNu.*HT.*).*{year}')

        data_err_opts = {
            'linestyle':'none',
            'marker': '.',
            'markersize': 10.,
            'color':'k',
            'elinewidth': 1,
        }

        fig, ax, rax = fig_ratio()

        hist.plot1d(
            h[data],
            ax=ax,
            overlay='dataset',
            error_opts=data_err_opts
        )

        hist.plot1d(
            h[mc],
            ax=ax,
            stack=True,
            overlay='dataset',
            clear=False
        )

        # Plot ratio
        hist.plotratio(
            h[data].integrate('dataset'),
            h[mc].integrate('dataset'),
            ax=rax,
            unc='num',
            denom_fill_opts={},
            guide_opts={},
            error_opts=data_err_opts
        )
        
        ax.text(0., 1., '$\\bf{CMS}$ internal',
                fontsize=14,
                horizontalalignment='left',
                verticalalignment='bottom',
                transform=ax.transAxes
               )

        ax.text(1., 1., f'VBF, {lumi(year):.1f} fb$^{{-1}}$ ({year})',
                fontsize=14,
                horizontalalignment='right',
                verticalalignment='bottom',
                transform=ax.transAxes
               )

        ax.legend(title='Single muon CR')

        rax.set_ylabel('Data / MC')
        rax.grid(True)
        rax.set_ylim(0.5,1.5)

        outdir = f'./output/{outtag}'
        if not os.path.exists(outdir):
            os.makedirs(outdir)

        outpath = pjoin(outdir, f'cr_1m_vbf_{distribution}_{year}.pdf')
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

    distributions = [
        'npv',
        'npv_nopu',
        'rho_all',
        'rho_all_nopu',
    ]

    for distribution in distributions:
        make_plot(acc, outtag, distribution=distribution)

if __name__ == '__main__':
    main()