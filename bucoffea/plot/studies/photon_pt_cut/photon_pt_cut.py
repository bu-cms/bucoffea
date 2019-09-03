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
    outdir = './output/photon_pt_cut/'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    for year in [2017,2018]:
        fig = plt.gcf()
        fig.clf()
        ax = plt.gca()
        h = copy.deepcopy(acc['photon_pt0_recoil'])
        h=h.rebin(h.axis('pt'), hist.Bin("pt",r"$p_{T}^{\gamma}$ (GeV)", [0,175,215,10000]))
        h=h.rebin(h.axis('recoil'),hist.Bin('recoil','recoil',list(range(200,500,50)) + list(range(500,1000,100)) + list(range(1000,2000,250))))
        h = merge_extensions(h, acc, reweight_pu=False)
        scale_xs_lumi(h)
        h = merge_datasets(h)


        # hlow = h.integrate(h.axis('pt'),)
        pprint(h.axis('dataset').identifiers())
        # h = h.integrate(h.axis('dataset'),f'GJets_HT_MLM_{year}')
        h = h.integrate(h.axis('dataset'),f'GJets_HT_MLM_{year}')
        h = h.integrate(h.axis('region'),'tr_g_notrig_num')
        pprint(h)
        hist.plot1d(
            h,
            overlay='pt',
            # error_opts=data_err_opts,
            ax=ax,
            overflow='all',
            clear=False)
        
        ax.set_ylim(0,2e5)
        ax.set_xlim(200,500)
        ax.set_ylabel('Expected GJets events (a.u.)')
        # rax.set_ylim(0.9,1.6)
        # ax.set_yscale('log')
        leg=ax.legend(['< 175', '175 - 215', '> 215'],title='Photon $p_{T}$')
        # for i, pdf in enumerate(h.axis('pdf').identifiers()):
        #     if str(pdf)=='none':
        #         continue
        #     leg.get_texts()[i].set_text(str(pdf))

        ax.text(0.97, 0.65, 'Photon CR, no trigger applied',
                fontsize=10,
                horizontalalignment='right',
                verticalalignment='bottom',
                transform=ax.transAxes
        )
        ax.plot([250,250],[0,1e8],'--',color='grey')
        
        fig.savefig(pjoin(outdir,f'photon_pt_cut_{year}.pdf'))
        plt.close(fig)

def main():
    acc = acc_from_dir('input/photon_pt_cut/')
    pdf_plot(acc)

if __name__ == "__main__":
    main()
