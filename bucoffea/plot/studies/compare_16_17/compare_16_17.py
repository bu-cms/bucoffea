#!/usr/bin/env python

import os
import re

from coffea.util import load

from bucoffea.plot.stack_plot import make_plot
from bucoffea.plot.util import acc_from_dir


def main():

    indir = "./input/met_2016"

    acc = acc_from_dir(indir)

    mc = re.compile(f'MET_2016')
    data = re.compile(f'MET_2017')

    region='cr_2m_j'
    for distribution in ['recoil', 'ak4_pt0']:
        make_plot(acc, region=region,distribution=distribution, year=2016, data=data, mc=mc, ylim=(1e-3,1e3), outdir=f'./output/{os.path.basename(indir)}')
    
    
if __name__ == "__main__":
    main()
