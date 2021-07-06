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

def plot_sf(file, histograms, year):
    outdir = './output/muon_sf'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    labels = {
        'NUM_LooseID_DEN_TrackerMuons_abseta_pt' : 'Muon ID SF (loose WP)',    
        'NUM_TightID_DEN_TrackerMuons_abseta_pt' : 'Muon ID SF (tight WP)', 
        'NUM_LooseRelIso_DEN_LooseID_abseta_pt' : 'Muon ISO SF (loose WP)',
        'NUM_TightRelIso_DEN_TightIDandIPCut_abseta_pt' : 'Muon ISO SF (tight WP)'
    }

    for histoname in histograms:
        fig, ax = plt.subplots()
        h = file[histoname]

        xedges, yedges = h.edges[0], h.edges[1]

        opts = {"cmap": "viridis"}
        pc = ax.pcolormesh(xedges, yedges, h.values.T, **opts)
        fig.colorbar(pc, ax=ax, label=labels[histoname])

        ax.set_xlabel(r'Muon $|\eta|$')
        ax.set_ylabel(r'Muon $p_T \ (GeV)$')

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
                txtformat = opts.pop("format", r"%.3f")
                ax.text(xcenter, ycenter, txtformat % h.values[ix, iy], **opts)

        outpath = pjoin(outdir, f'{histoname}_{year}.pdf')
        fig.savefig(outpath)
        plt.close(fig)
        print(f'File saved: {outpath}')

def main():
    files = [
        bucoffea_path('data/sf/muon/Efficiencies_muon_generalTracks_Z_Run2017_UL_ID.root'),
        bucoffea_path('data/sf/muon/Efficiencies_muon_generalTracks_Z_Run2017_UL_ISO.root'),
        bucoffea_path('data/sf/muon/Efficiencies_muon_generalTracks_Z_Run2018_UL_ID.root'),
        bucoffea_path('data/sf/muon/Efficiencies_muon_generalTracks_Z_Run2018_UL_ISO.root'),
    ]

    histos = {
        'ID' : [
            'NUM_LooseID_DEN_TrackerMuons_abseta_pt',
            'NUM_TightID_DEN_TrackerMuons_abseta_pt',
        ],
        'ISO' : [
            'NUM_LooseRelIso_DEN_LooseID_abseta_pt',
            'NUM_TightRelIso_DEN_TightIDandIPCut_abseta_pt'
        ]
    }

    for f in files:
        inrootfile = uproot.open(f)
        if f.endswith('_ID.root'):
            histograms = histos['ID']
        else:
            histograms = histos['ISO']

        year = re.findall('201\d', f)[0]

        plot_sf(uproot.open(f), histograms, year)

if __name__ == '__main__':
    main()