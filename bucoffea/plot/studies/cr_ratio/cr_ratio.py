#!/usr/bin/env python

import os
import re
import copy
from collections import defaultdict
from pprint import pprint

import matplotlib.ticker
import numpy as np
from coffea import hist
from coffea.util import load
from matplotlib import pyplot as plt

from bucoffea.plot.util import (acc_from_dir, lumi, merge_datasets,
                                merge_extensions, scale_xs_lumi)
from bucoffea.plot.stack_plot import Style, make_plot
from bucoffea.plot.cr_ratio_plot import cr_ratio_plot
from coffea.hist.plot import clopper_pearson_interval
pjoin = os.path.join

def main():
    infile=os.path.abspath('input/2019-09-09_gen_dilep_sf/')
    acc = acc_from_dir(infile)

    for year in [2017,2018]:
        mc = {
            'cr_1m_j' : re.compile(f'(TTJets.*FXFX.*|Diboson.*|ST.*|QCD_HT.*|.*DYJetsToLL_M-50_HT_MLM.*|W.*HT.*).*{year}'),
            'cr_1e_j' : re.compile(f'(TTJets.*FXFX.*|Diboson.*|ST.*|QCD_HT.*|.*DYJetsToLL_M-50_HT_MLM.*|W.*HT.*).*{year}'),
            'cr_2m_j' : re.compile(f'(EW.*|TTJets.*FXFX.*|Diboson.*|ST.*|QCD_HT.*|.*DYJetsToLL_M-50_HT_MLM.*).*{year}'),
            'cr_2e_j' : re.compile(f'(EW.*|TTJets.*FXFX.*|Diboson.*|ST.*|QCD_HT.*|.*DYJetsToLL_M-50_HT_MLM.*).*{year}'),
            'cr_g_j' : re.compile(f'(GJets.*|QCD_HT.*|W.*HT.*).*{year}'),
        }
        cr_ratio_plot(acc, year=year,tag='losf',outdir=f'./output/{os.path.basename(infile)}', mc=mc)
    
    for year in [2017,2018]:
        mc = {
                'cr_1m_j' : re.compile(f'(TTJets.*FXFX.*|Diboson.*|ST.*|QCD_HT.*|.*DY.*FXFX.*|W.*FXFX.*).*{year}'),
                'cr_1e_j' : re.compile(f'(TTJets.*FXFX.*|Diboson.*|ST.*|QCD_HT.*|.*DY.*FXFX.*|W.*FXFX.*).*{year}'),
                'cr_2m_j' : re.compile(f'(EW.*|TTJets.*FXFX.*|Diboson.*|ST.*|QCD_HT.*|.*DY.*FXFX.*).*{year}'),
                'cr_2e_j' : re.compile(f'(EW.*|TTJets.*FXFX.*|Diboson.*|ST.*|QCD_HT.*|.*DY.*FXFX.*).*{year}'),
                'cr_g_j' : re.compile(f'(GJets.*|QCD_HT.*|W.*FXFX.*).*{year}'),
        }
        cr_ratio_plot(acc, year=year,tag='nlo',outdir=f'./output/{os.path.basename(infile)}', mc=mc)

if __name__ == "__main__":
    main()

