#!/usr/bin/env python

import os
from bucoffea.plot.util import acc_from_dir
from klepto.archives import dir_archive
import argparse
from datetime import datetime
pjoin = os.path.join

def parse_commandline():
    parser = argparse.ArgumentParser()
    parser.add_argument('inpath', type=str, help='Input path to use.')
    parser.add_argument('--channel', type=str, help='Channel to make inputs for.', default='monojet')
    parser.add_argument('--unblind', action='store_true', help='Include signal region data')
    parser.add_argument('--years', nargs='*', default=[2017, 2018], help='The years to prepare the limit input for')
    parser.add_argument('--eoyxs', action='store_true', help='Use the EOY XS for normalization, otherwise use UL XS')
    parser.add_argument('--one_fifth_unblind', action='store_true', help='1/5th unblinding: Scale the MC in signal region by 1/5')
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

    outdir = pjoin('./output/',list(filter(lambda x:x,args.inpath.split('/')))[-1])

    # Store the command line arguments in the INFO.txt file
    try:
        os.makedirs(outdir)
    except FileExistsError:
        pass
    
    infofile = pjoin(outdir, 'INFO.txt')
    with open(infofile, 'w+') as f:
        f.write(f'Limit input creation: {datetime.now().strftime("%m/%d/%Y, %H:%M:%S")}\n')
        f.write('Command line arguments:\n\n')
        cli = vars(args)
        for arg, val in cli.items():
            f.write(f'{arg}: {val}\n')

    for channel in args.channel.split(','):
        print(channel)
        if channel == 'monojet':
            from legacy_monojet import legacy_limit_input_monojet
            legacy_limit_input_monojet(acc, outdir=outdir, unblind=args.unblind)
        elif channel == 'monov':
            from legacy_monov import legacy_limit_input_monov
            legacy_limit_input_monov(acc, outdir=outdir, unblind=args.unblind)
        elif channel == 'vbfhinv':
            from legacy_vbf import legacy_limit_input_vbf
            legacy_limit_input_vbf(acc, outdir=outdir, unblind=args.unblind, years=args.years, ulxs=not args.eoyxs, one_fifth_unblind=args.one_fifth_unblind)


if __name__ == "__main__":
    main()
