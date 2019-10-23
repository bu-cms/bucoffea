#!/usr/bin/env python

import os
import sys
from bucoffea.plot.util import acc_from_dir
from bucoffea.limit.legacy import legacy_limit_input

def main():
    inpath = os.path.abspath(sys.argv[1])
    if not os.path.isdir(inpath):
        raise RuntimeError(f"Commandline argument is not a directory: {inpath}")
    acc = acc_from_dir(inpath)

    legacy_limit_input  (acc)









if __name__ == "__main__":
    main()
