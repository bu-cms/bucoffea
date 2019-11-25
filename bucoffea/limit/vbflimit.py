#!/usr/bin/env python

import os
import sys
from bucoffea.plot.util import acc_from_dir
from bucoffea.limit.vbflegacy import legacy_limit_input_vbf
from klepto.archives import dir_archive


def main():
    inpath = os.path.abspath(sys.argv[1])
    if not os.path.isdir(inpath):
        raise RuntimeError(f"Commandline argument is not a directory: {inpath}")

    klepto = True
    if klepto:
        acc = dir_archive(inpath, serialized=True, compression=0, memsize=1e3)
        acc.load('mjj')
        acc.load('sumw')
        acc.load('sumw_pileup')
        acc.load('nevents')
    else:
        acc = acc_from_dir(inpath)

    legacy_limit_input_vbf(acc)


if __name__ == "__main__":
    main()
