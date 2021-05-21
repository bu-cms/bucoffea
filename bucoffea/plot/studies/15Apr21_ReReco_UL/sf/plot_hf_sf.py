#!/usr/bin/env python
import os
import re
import sys
import uproot
import numpy as np

from matplotlib import pyplot as plt
from coffea import hist
from bucoffea.helpers.paths import bucoffea_path
from pprint import pprint

pjoin = os.path.join

def main():
    inpath = bucoffea_path('data/sf/hfmask/EffcySFandNoisePassRate.root')
    infile = uproot.open(inpath)

    for year in [2017, 2018]:
        h = infile[f'EfficiencySF_{year}']
        xedges, yedges = h.edges

        # Remove the 2.95 < |eta| < 2.99 bin        
        xedges = xedges[1:]
        vals = h.values[1:, :]
        err = np.sqrt(h.variances[1:, :])

        fig, ax = plt.subplots()
        opts = {"cmap": "viridis"}
        pc = ax.pcolormesh(xedges, yedges, vals.T, **opts)
        fig.colorbar(pc, ax=ax, label='HF Scale Factor')

        xcenters = 0.5 * (xedges + np.roll(xedges, -1))[:-1]
        ycenters = 0.5 * (yedges + np.roll(yedges, -1))[:-1]

        for ix, xcenter in enumerate(xcenters):
            for iy, ycenter in enumerate(ycenters):
                opts = {
                    "horizontalalignment": "center",
                    "verticalalignment": "center",
                    "fontsize" : 8.
                }
                opts['color'] = 'white' if vals[ix, iy] < 0.5 else 'black'
                ax.text(xcenter, ycenter, f'${vals[ix, iy]:.3f}$ \n $\\pm$ \n ${err[ix, iy]:.3f}$', **opts)

        ax.set_xlabel(r'Jet $|\eta|$', fontsize=14)
        ax.set_ylabel(r'Jet $p_T \ (GeV)$', fontsize=14)

        ax.text(0.,1.,'HF Scale Factors',
            fontsize=14,
            ha='left',
            va='bottom',
            transform=ax.transAxes
        )

        ax.text(1.,1.,year,
            fontsize=14,
            ha='right',
            va='bottom',
            transform=ax.transAxes
        )

        outdir = './output/hf_sf'
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        outpath = pjoin(outdir, f'hf_sf_{year}.pdf')
        fig.savefig(outpath)
        plt.close(fig)

        print(f'File saved: {outpath}')

if __name__ == '__main__':
    main()