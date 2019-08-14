#!/usr/bin/env python

import sys
import os
from matplotlib import pyplot as plt
from coffea import hist
from coffea.util import load
import argparse
pjoin = os.path.join

def commandline():
    parser = argparse.ArgumentParser(prog='Quick and dirty plot dumper from coffea output.')
    parser.add_argument('file', type=str, help='Input file to use.')
    parser.add_argument('--outpath', type=str, help='Path to save output under.', default='./out')
    parser.add_argument('--region', type=str, default='inclusive', help='Region to plot.')
    args = parser.parse_args()
    return args


def main():
    args = commandline()
    output = load(args.file)

    h = output['weights']
    fill_opts = {
        'edgecolor': (0,0,0,0.3),
        'alpha': 0.8
    }
    fig, ax, _ = hist.plot1d(
        h. project('dataset').project("region", args.region),
        overlay='weight_type',
        overflow='all',
        fill_opts=fill_opts
        )
    fig.suptitle("Weights")
    # ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim(0.1, 1e8)
    fig.savefig(pjoin(args.outpath, "weights.pdf"))
    plt.close(fig)

if __name__ == "__main__":
    main()