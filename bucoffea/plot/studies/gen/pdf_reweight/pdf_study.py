#!/usr/bin/env python

import os
import re
import copy
from collections import defaultdict
from pprint import pprint

import numpy as np
from coffea import hist
from coffea.util import load
from matplotlib import pyplot as plt

from bucoffea.execute.dataset_definitions import short_name
from bucoffea.helpers.dataset import is_data
from bucoffea.helpers.paths import bucoffea_path
from bucoffea.plot.util import (acc_from_dir, lumi, merge_datasets,
                                merge_extensions, scale_xs_lumi,fig_ratio)

data_err_opts = {
    'linestyle':'none',
    'marker': '.',
    'markersize': 10.,
    'color':'k',
    'elinewidth': 0,
    'emarker': ''
}

colors = {
    "303600" : 'darkslategrey',
    "263000" : 'skyblue',
    "262000" : 'darkorange',
    'none' : 'magenta'
        }

pjoin = os.path.join
def pdf_plot(acc):
    outdir = './output/pdfstudy/'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    datasets = [
        'WJetsToLNu_HT_MLM_2017',
        'DYJetsToLL_M-50_HT_MLM_2017' ,
    ]
    for ds in datasets:
        fig, ax, rax = fig_ratio()
        h = acc['gen_vpt']
        h=h.rebin(h.axis('vpt'), hist.Bin("vpt",r"$p_{T}^{V}$ (GeV)", 10, 0, 2000))
        h = merge_extensions(h, acc, reweight_pu=False)
        scale_xs_lumi(h)
        h = merge_datasets(h)

        h = h.project(h.axis('dataset'),ds)
        
        

        for pdf in h.axis('pdf').identifiers():

            if str(pdf)=='none':
                continue
            data_err_opts['color'] = colors[str(pdf)]
            hist.plot1d(
                h.project('pdf',pdf),
                # overlay='pdf',
                error_opts=data_err_opts,
                ax=ax,
                overflow='all',
                clear=False)
            
            hist.plotratio(h.project('pdf',pdf), h.project('pdf','none'),
                ax=rax,
                denom_fill_opts={},
                guide_opts={},
                unc='num',
                overflow='all',
                error_opts=data_err_opts,
                clear=False,
                )
        ax.set_ylim(1e-3,1e8)
        rax.set_ylim(0.9,1.6)
        ax.set_yscale('log')
        leg=ax.legend()
        for i, pdf in enumerate(h.axis('pdf').identifiers()):
            if str(pdf)=='none':
                continue
            leg.get_texts()[i].set_text(str(pdf))
        fig.savefig(pjoin(outdir,f'{ds}.pdf'))
        plt.close(fig)

def main():
    acc = acc_from_dir('input/pdfweight/')
    pdf_plot(acc)

if __name__ == "__main__":
    main()
