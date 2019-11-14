#!/usr/bin/env python

import os
import sys
from bucoffea.plot.util import acc_from_dir
from bucoffea.limit.monovlegacy import legacy_limit_input
from klepto.archives import dir_archive


def main():
    if len(sys.argv) > 1:
        inpath = os.path.abspath(sys.argv[1])
    else:
        inpath = '../plot/input/merged'
    if not os.path.isdir(inpath):
        raise RuntimeError(f"Commandline argument is not a directory: {inpath}")

    klepto = True
    if klepto:
        acc = dir_archive(inpath, serialized=True, compression=0, memsize=1e3)
        acc.load('recoil')
        acc.load('sumw')
        acc.load('sumw_pileup')
        acc.load('nevents')
    else:
        acc = acc_from_dir(inpath)

    legacy_limit_input(acc)


if __name__ == "__main__":
    main()
