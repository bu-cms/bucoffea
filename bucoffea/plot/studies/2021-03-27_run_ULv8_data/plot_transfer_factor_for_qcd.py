#!/usr/bin/env python

import os
import sys
import re
import uproot
import numpy as np
import mplhep as hep

from matplotlib import pyplot as plt
from bucoffea.helpers.paths import bucoffea_path

pjoin = os.path.join

def main():
    infile = uproot.open(bucoffea_path('data/sf/hfmask/EffcySFandNoisePassRate.root'))

    h = infile['NoisePassingRate_2017']
    # Omitting the first eta bin, 2.95-2.99. Not using it anyway
    tf = h.values[1:,:] / (1 - h.values)[1:,:]

    fig, ax = plt.subplots()
    xbins = h.edges[0][1:]
    ybins = h.edges[1]
    hep.hist2dplot(tf, xbins=xbins, ybins=ybins, ax=ax)

    ax.set_xlabel(r'Jet $|\eta|$', fontsize=14)
    ax.set_ylabel(r'Jet $p_T \ (GeV)$', fontsize=14)

    ax.text(0,1,r'CR $\rightarrow$ SR Transfer Factor',
        fontsize=14,
        ha='left',
        va='bottom',
        transform=ax.transAxes
    )

    _xbin_centers = xbins[1:] - np.diff(xbins) / float(2)
    _ybin_centers = ybins[1:] - np.diff(ybins) / float(2)
    for ix, xc in enumerate(_xbin_centers):
        for iy, yc in enumerate(_ybin_centers):
            color = "black" if tf[ix, iy] > 0.5 else "lightgrey"
            ax.text(xc, yc, f'{tf[ix, iy]:.2f}', ha="center", va="center", color=color)

    outdir =  f'./output/qcd_tf'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outpath = pjoin(outdir, f'qcd_tf.pdf')
    fig.savefig(outpath)
    plt.close(fig)

if __name__ == '__main__':
    main()