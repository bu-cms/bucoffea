#!/usr/bin/env python

import os
import sys
import re
import uproot
import mplhep as hep
import numpy as np

from bucoffea.plot.util import merge_datasets, merge_extensions, scale_xs_lumi, fig_ratio
from coffea import hist
from matplotlib import pyplot as plt
from pprint import pprint

pjoin = os.path.join

BINNINGS = {
    'leadak4_pt' : list(range(0,800,10)),
    'trailak4_pt' : list(range(0,800,10)),
    'recoil_pt' : list(range(200,800,10)),
}

def make_plot(tree, distribution):
    vals_nom = tree[distribution].array()
    vals_jerUp = tree[f'{distribution}_jerUp'].array()
    vals_jerDown = tree[f'{distribution}_jerDown'].array()

    h_nom, edges = np.histogram(vals_nom, bins=BINNINGS[distribution])
    h_jerUp, edges = np.histogram(vals_jerUp, bins=BINNINGS[distribution])
    h_jerDown, edges = np.histogram(vals_jerDown, bins=BINNINGS[distribution])

    fig, ax = plt.subplots()
    hep.histplot(h_nom, edges, ax=ax, label='Nominal')
    hep.histplot(h_jerUp, edges, ax=ax, label='JER up')
    hep.histplot(h_jerDown, edges, ax=ax, label='JER down')

    ax.legend()

    outdir = './output/zvv_jer'
    try: 
        os.makedirs(outdir)
    except FileExistsError:
        pass
    outpath = pjoin(outdir, f'zvv_jer_{distribution}.pdf')
    fig.savefig(outpath)
    plt.close(fig)

    print(f'File saved: {outpath}')

def main():
    infile = uproot.open('tree_ZJetsToNuNu_HT-MLM_2017.root')
    tree = infile['sr_vbf']

    distributions = [
        'leadak4_pt',
        'trailak4_pt',
        'recoil_pt',
    ]

    for distribution in distributions:
        make_plot(tree, distribution)

if __name__ == '__main__':
    main()