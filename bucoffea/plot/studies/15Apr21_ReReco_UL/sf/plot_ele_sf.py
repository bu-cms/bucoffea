#!/usr/bin/env python

import os
import sys
import re
import uproot
import numpy as np
import mplhep as hep

from matplotlib import pyplot as plt
from bucoffea.helpers.paths import bucoffea_path
from bucoffea.plot.util import lumi
from pprint import pprint

pjoin = os.path.join

labels = {
    'reco' : 'Electron Reconstruction SF',
    'reco_ptbelow20' : 'Electron Reconstruction SF',
    'id_veto' : 'Veto ID SF',
    'id_tight' : 'Tight ID SF',
}

def plot_ele_sf(inpath, histogram, year, sftype):
    outdir = './output/ele_sf'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    f = uproot.open(inpath)
    h = f[histogram]

    xedges, yedges = h.edges[0], h.edges[1]
    
    fig, ax = plt.subplots()
    opts = {"cmap": "viridis"}
    pc = ax.pcolormesh(xedges, yedges, h.values.T, **opts)
    fig.colorbar(pc, ax=ax, label=labels[sftype])

    ax.set_xlabel(r'Electron Supercluster $\eta$')
    ax.set_ylabel(r'Electron $p_T \ (GeV)$')

    ax.text(0., 1., '$\\bf{CMS}$ internal',
            fontsize=14,
            ha='left',
            va='bottom',
            transform=ax.transAxes
            )

    ax.text(1.,1., f'{year}, {lumi(int(year)):.1f} $fb^{{-1}}$ (13 TeV)',
        fontsize=14,
        ha='right',
        va='bottom',
        transform=ax.transAxes
    )

    # Print numerical values in the plot
    xcenters = 0.5 * (xedges + np.roll(xedges, -1))[:-1]
    ycenters = 0.5 * (yedges + np.roll(yedges, -1))[:-1]

    for ix, xcenter in enumerate(xcenters):
        for iy, ycenter in enumerate(ycenters):
            opts = {
                "horizontalalignment": "center",
                "verticalalignment": "center",
            }
            opts["color"] = (
                "black" if pc.norm(h.values[ix, iy]) > 0.5 else "lightgrey"
            )
            if 'RECO' in sftype:
                if ix in [2,9]:
                    continue
            else:
                if ix in [2,7]:
                    continue
            txtformat = opts.pop("format", r"%.2f")
            ax.text(xcenter, ycenter, txtformat % h.values[ix, iy], **opts)

    outpath = pjoin(outdir, f'electron_{sftype}_{year}.pdf')
    fig.savefig(outpath)
    plt.close(fig)
    print(f'File saved: {outpath}')

def main():
    files = {
        'reco' : {
            2017: bucoffea_path('data/sf/egamma/ul/egammaEffi_ptAbove20_txt_EGM2D_UL2017.root'),
            2018: bucoffea_path('data/sf/egamma/ul/egammaEffi_ptAbove20_txt_EGM2D_UL2018.root'),
        },
        'reco_ptbelow20' : {
            2017: bucoffea_path('data/sf/egamma/ul/egammaEffi_ptBelow20_txt_EGM2D_UL2017.root'),
            2018: bucoffea_path('data/sf/egamma/ul/egammaEffi_ptBelow20_txt_EGM2D_UL2018.root'),
        },
        'id_veto' : {
            2017: bucoffea_path('data/sf/egamma/ul/egammaEffi_EGM2D_Veto_UL17_fix.root'),
            2018: bucoffea_path('data/sf/egamma/ul/egammaEffi_EGM2D_Veto_UL18_fix.root'),
        },
        'id_tight': {
            2017: bucoffea_path('data/sf/egamma/ul/egammaEffi_txt_EGM2D_Tight_UL17.root'),   
            2018: bucoffea_path('data/sf/egamma/ul/egammaEffi_txt_Ele_Tight_EGM2D_UL18.root'),   
        }
    }
    
    histoname = 'EGamma_SF2D'

    for sftype, entry in files.items():
        years = entry.keys()
        for year in years:
            inpath = entry[year]

            plot_ele_sf(inpath, histogram=histoname, year=year, sftype=sftype)


if __name__ == '__main__':
    main()