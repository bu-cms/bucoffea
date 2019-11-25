#!/usr/bin/env python

import os
import sys
from bucoffea.plot.util import acc_from_dir
from bucoffea.limit.legacy import legacy_limit_input
from klepto.archives import dir_archive
# from bucoffea.limit.new_limit import monojet
import argparse

def parse_commandline():
    parser = argparse.ArgumentParser()
    parser.add_argument('inpath', type=str, help='Input path to use.')
    parser.add_argument('--channel', type=str, choices=['monojet','monov','vbfhinv'], help='Channel to make inputs for.', default='monojet')
    args = parser.parse_args()

    if not os.path.isdir(args.inpath):
        raise RuntimeError(f"Commandline argument is not a directory: {args.inpath}")

    return args

def main():
    args = parse_commandline()

    klepto = True
    if klepto:
        acc = dir_archive(args.inpath, serialized=True, compression=0, memsize=1e3)
        acc.load('recoil')
        acc.load('sumw')
        acc.load('sumw_pileup')
        acc.load('nevents')
    else:
        acc = acc_from_dir(args.inpath)

    if args.channel == 'monojet':
        from legacy_monojet import legacy_limit_input_monojet
        legacy_limit_input_monojet(acc)
    elif args.channel == 'monov':
        from legacy_monov import legacy_limit_input_monov
        legacy_limit_input_monov(acc)
    elif args.channel == 'vbfhinv':
        from legacy_vbf import legacy_limit_input_vbf
        legacy_limit_input_vbf(acc)


if __name__ == "__main__":
    main()
