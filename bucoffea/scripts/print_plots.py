#!/usr/bin/env python

import sys
from bucoffea.plot.debug import debug_plot_output
from coffea.util import load, save
import argparse
def commandline():
    parser = argparse.ArgumentParser(prog='Quick and dirty plot dumper from coffea output.')
    parser.add_argument('file', type=str, help='Input file to use.')
    parser.add_argument('--outpath', type=str, help='Path to save output under.', default='./out')
    parser.add_argument('--region', type=str, default='inclusive', help='Region to plot.')
    args = parser.parse_args()
    return args

def main():
    args = commandline()
    acc = load(args.file)
    debug_plot_output(acc, args.region, args.outpath)

if __name__ == "__main__":
    main()
