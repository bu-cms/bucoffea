#!/usr/bin/env python

import os
import re
import sys
import warnings
import uproot
import numpy as np
import mplhep as hep

from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
from bucoffea.plot.util import fig_ratio
from tqdm import tqdm
from pprint import pprint

pjoin = os.path.join

warnings.filterwarnings('ignore')

def compare_limit_inputs(infiles):
    '''Compare the prefit MC shapes between the two versions (MTR 2017).'''
    keys = [key.decode('utf-8').replace(';1', '') for key in infiles['01Apr21'].keys()]

    for key in tqdm(keys):
        h_01Apr = infiles['01Apr21'][key]
        h_07Apr = infiles['07Apr21'][key]

        fig, ax, rax = fig_ratio()
        hep.histplot(h_01Apr.values, h_01Apr.edges, ax=ax, label='01 Apr')
        hep.histplot(h_07Apr.values, h_07Apr.edges, ax=ax, label='07 Apr')

        ax.set_yscale('log')
        ax.set_ylim(1e-2,1e6)
        ax.legend(title='Version')

        ax.text(0.,1.,key,
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

        xcenters = (0.5 * (h_01Apr.edges + np.roll(h_01Apr.edges, -1)))[:-1]

        rax.plot(xcenters, h_01Apr.values / h_07Apr.values, **data_err_opts)

        rax.grid(True)
        rax.set_ylim(0,2)
        rax.set_ylabel('01Apr / 07Apr')

        outdir = f'./output/compare_limit_inputs/afterFix'
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        outpath = pjoin(outdir, f'{key}.png')
        fig.savefig(outpath)
        plt.close(fig)

def main():
    infiles = {
        '01Apr21' : uproot.open('input/merged_2021-04-01_vbfhinv_ULv8_DATA_03Sep20v7_MC_noRegularCleaningCuts_correctEta_4_0_cut/legacy_limit_vbf_2017.root'),
        # '07Apr21' : uproot.open('input/merged_2021-04-07_vbfhinv_ULv8_DATA_one_fifth_unblind_03Sep20v7_MC_withHFWeight_jetPt80/legacy_limit_vbf_2017.root'),
        '07Apr21' : uproot.open('input/merged_2021-04-08_vbfhinv_03Sep20v7_MC_noEndcapReweight_withHFReweight_jetPt80/legacy_limit_vbf_2017.root'),
    }

    compare_limit_inputs(infiles)

if __name__ == '__main__':
    main()