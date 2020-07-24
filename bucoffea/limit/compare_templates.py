#!/usr/bin/env python
import argparse
import os
import re
import sys

import mplhep as hep
import numpy as np
import tqdm
import uproot
from matplotlib import pyplot as plt
from tabulate import tabulate

from bucoffea.plot.util import fig_ratio

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



def parse_commandline():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "fname1",
        type=str,
        help="First template file.",
    )
    parser.add_argument(
        "fname2",
        type=str,
        help="Second template file.",
    )

    parser.add_argument(
        "--outdir",
        "-o",
        type=str,
        default="./plots",
        help="The output directory to use.",
    )
    parser.add_argument(
        "--rlim",
        type=str,
        default="0.9,1.1",
        help="Lower/upper limit of ratio panel. Comma separated.",
    )

    args = parser.parse_args()
    if "INDIR" in args.outdir:
        args.outdir = args.outdir.replace("INDIR", args.indir)


    return args

def main():
    """
    A script to easily compare template files between different runs.

    Usage: ./compare_templates.py /path/to/first/template_file.root /path/to/second/template_file.root

    The script will loop over all templates in the files and create a comparison plot for each of them.
    All plots are dumped into a folder for inspection.
    """

    args = parse_commandline()

    # Based on input locations, derive tag names to identify the files
    tag1 = os.path.basename(os.path.dirname(args.fname1))
    tag2 = os.path.basename(os.path.dirname(args.fname2))

    # Convert to dictionary
    h1 = make_dict(args.fname1)
    h2 = make_dict(args.fname2)

    # Make sure the two files have consistent keys
    # assert(h1.keys()==h2.keys())

    # Create plot folder
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    # Do the actual plotting
    table = []
    for key in tqdm.tqdm(h1.keys()):
        if key not in h2:
            print("Found missing key ", key)
        fig, ax, rax = fig_ratio()
        x = 0.5 * np.sum(h1[key].bins,axis=1)
        edges = np.unique(h1[key].bins)

        table.append([key, np.sum(h1[key].values), np.sum(h2[key].values)])
        hep.histplot(
            h1[key].values,
            edges,
            ax=ax,
            label=f"{tag1}, Integral={np.sum(h1[key].values):.1f}",
            color='navy',
            )

        hep.histplot(
            h2[key].values,
            edges,
            yerr=np.sqrt(h2[key].variances),
            ax=ax,
            label=f"{tag2}, Integral={np.sum(h2[key].values):.1f}",
            color='crimson',
            marker='o',
            markersize=5,
            histtype='errorbar'
            )

        ax.legend()

        # Bottom panel: ratio plot
        valid = h1[key].values!=0

        rax.errorbar(
                x[valid],
                h2[key].values[valid] / h1[key].values[valid],
                np.sqrt(h2[key].variances[valid]) / h1[key].values[valid],
                linestyle='none',
                marker='o',
                color="crimson",
                )

        # Add indicators for bins where we could not calculate the ratio
        if np.any(~valid):
            rax.plot(
                    0.5*np.sum(h1[key].bins,axis=1)[~valid],
                    np.ones(np.sum(~valid)),
                    'x',
                    color="k",
                    fillstyle="none"
                    )

        # Aesthetics
        ax.set_title(key)
        rax.set_ylim(*map(float, args.rlim.split(",")))
        rax.set_xlabel("Recoil (GeV)")
        rax.set_ylabel("Ratio")
        ax.set_ylabel("Events / bin")
        ax.set_yscale("log")

        try:
            ax.set_ylim(
                0.5*min(h2[key].values[h2[key].values>0]),
                1.5*max(h2[key].values),
            )
        except ValueError:
            continue
        rax.grid(linestyle='--')
        fig.savefig(pjoin(args.outdir, f"{key}.png"))
        plt.close(fig)
    print(tabulate(table))
if __name__ == "__main__":
    main()
