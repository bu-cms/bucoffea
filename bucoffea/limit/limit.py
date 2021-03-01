#!/usr/bin/env python

import os
from bucoffea.plot.util import acc_from_dir
from klepto.archives import dir_archive
import argparse
pjoin = os.path.join

def parse_commandline():
    parser = argparse.ArgumentParser()
    parser.add_argument('inpath', type=str, help='Input path to use.')
    parser.add_argument('--channel', type=str, help='Channel to make inputs for.', default='monojet')
    parser.add_argument('--unblind', action='store_true', help='Include signal region data')
    parser.add_argument('--nlo', action='store_true', help='Use NLO V samples')
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
        acc.load('mjj')
        acc.load('sumw')
        acc.load('sumw_pileup')
        acc.load('nevents')
    else:
        acc = acc_from_dir(args.inpath)

    args.outdir = pjoin('./output/',list(filter(lambda x:x,args.inpath.split('/')))[-1])
    for channel in args.channel.split(','):
        print(channel)
        if channel == 'monojet':
            from legacy_monojet import legacy_limit_input_monojet
            legacy_limit_input_monojet(acc, args)
        elif channel == 'monov':
            from legacy_monov import legacy_limit_input_monov
            legacy_limit_input_monov(acc, args)
        elif channel == 'vbfhinv':
            from legacy_vbf import legacy_limit_input_vbf
            legacy_limit_input_vbf(acc, outdir=args.outdir, unblind=args.unblind)


if __name__ == "__main__":
    main()
