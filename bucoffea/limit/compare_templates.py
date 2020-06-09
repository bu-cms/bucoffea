#!/usr/bin/env python
import numpy as np
from matplotlib import pyplot as plt 
from bucoffea.plot.util import fig_ratio
import uproot
import sys
import os
import re
import tqdm
pjoin = os.path.join


def trim(key):
    """
    Trim unneeded decorators off of ROOT key names.
    """
    key = key.decode("utf-8")
    key = re.sub(";.*","",key)
    return key

def make_dict(fname):
    """
    Reads our template files and saves histograms into a convenient ditionary.
    """
    hists = {}
    f = uproot.open(fname)

    # First level of file structure is just subfolders
    # for each of the categories
    for category in map(trim, f.keys()):
        for hname in map(trim, f[category].keys()):
            # Read histogram
            h = f[f"{category}/{hname}"]

            # Save it into our dictionary
            hists[f"{category}_{hname}"] = h
    return hists

def main():
    """
    A script to easily compare template files between different runs.

    Usage: ./compare_templates.py /path/to/first/template_file.root /path/to/second/template_file.root

    The script will loop over all templates in the files and create a comparison plot for each of them.
    All plots are dumped into a folder for inspection.
    """

    # The two input files to compare are read off the commandline
    fname1 = sys.argv[1]
    fname2 = sys.argv[2]

    # Based on their locations, derive tag names to identify the files
    tag1 = os.path.basename(os.path.dirname(fname1))
    tag2 = os.path.basename(os.path.dirname(fname2))

    # Convert to dictionary
    h1 = make_dict(fname1)
    h2 = make_dict(fname2)

    # Make sure the two files have consistent keys
    assert(h1.keys()==h2.keys())

    # Create plot folder
    outdir = "./plots/"
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # Do the actual plotting
    for key in tqdm.tqdm(h1.keys()):
        fig, ax, rax = fig_ratio()
        x = np.sum(h1[key].bins,axis=1)
        
        # Top panel with asbolute comparison
        ax.plot(
                x, 
                h1[key].values,
                'o',
                color="dodgerblue",
                label=f"{tag1}, Integral={np.sum(h1[key].values):.1f}"
                )
        ax.plot(
                x,
                h2[key].values,
                's',
                color="crimson",
                fillstyle="none",
                label=f"{tag2}, Integral={np.sum(h2[key].values):.1f}"
                )
        ax.legend()

        # Bottom panel: ratio plot
        valid = h1[key].values!=0

        rax.plot(
                x, 
                h2[key].values[valid] / h1[key].values[valid],
                's',
                color="crimson",
                fillstyle="none"
                )

        # Add indicators for bins where we could not calculate the ratio
        if np.any(~valid):
            rax.plot(
                    np.sum(h1[key].bins,axis=1)[~valid],
                    np.ones(np.sum(~valid)),
                    'x',
                    color="k",
                    fillstyle="none"
                    )

        # Aesthetics
        ax.set_title(key)
        rax.set_ylim(0,2)
        rax.set_ylabel("Recoil (GeV)")
        ax.set_ylabel("Events / bin")
        fig.savefig(pjoin(outdir, f"{key}.png"))
        plt.close(fig)

if __name__ == "__main__":
    main()