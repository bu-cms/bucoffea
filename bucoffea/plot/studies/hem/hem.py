#!/usr/bin/env python

import os
import re
import copy
from collections import defaultdict
from pprint import pprint

import numpy as np
from coffea import hist
from coffea.util import load
from matplotlib import pyplot as plt

from bucoffea.execute.dataset_definitions import short_name
from bucoffea.helpers.dataset import is_data
from bucoffea.helpers.paths import bucoffea_path
from bucoffea.plot.util import (acc_from_dir, lumi, merge_datasets,
                                merge_extensions, scale_xs_lumi)
from bucoffea.plot.stack_plot import make_plot

def main():
    
    indir = bucoffea_path("plot/input/21Aug19_v2_newpu")

    acc = acc_from_dir(indir)

    for year in [2018]:
        data = re.compile(f'EGamma_{year}')
        mc = re.compile(f'(EW.*|TTJets.*|ZZ.*|ST.*|QCD_HT.*|WW.*|WZ.*|.*DYJetsToLL_M-50_HT_MLM.*|WJet.*HT.*){year}')
        region='cr_1e_j'
        for distribution in ['ak4_phi']:
            make_plot(acc, region=region,distribution=distribution, year=year, data=data, mc=mc, outdir=f'./output/{os.path.basename(indir)}')
    for year in [2017,2018]:
        data = re.compile(f'MET_{year}')
        mc = re.compile(f'(EW.*|TTJets.*|ZZ.*|ST.*|QCD_HT.*|WW.*|WZ.*|.*DYJetsToLL_M-50_HT_MLM.*|WJet.*HT.*){year}')
        region='cr_1m_j'
        for distribution in ['ak4_phi','recoil']:
            make_plot(acc, region=region,distribution=distribution, year=year, data=data, mc=mc, outdir=f'./output/{os.path.basename(indir)}')

if __name__ == "__main__":
    main()














# for distribution in ['ak4_pt0_eta0']:
#     make_plot(acc, region=region,distribution=distribution, year=year, data=data, mc=mc, integrate=('jeteta',0,1))
#     make_plot(acc, region=region,distribution=distribution, year=year, data=data, mc=mc, integrate=('jeteta',1,3))
#     make_plot(acc, region=region,distribution=distribution, year=year, data=data, mc=mc, integrate=('jetpt',100,200))
#     make_plot(acc, region=region,distribution=distribution, year=year, data=data, mc=mc, integrate=('jetpt',200,500))






# for year in [2017]:
#     data = re.compile(f'MET_{year}')
#     mc = re.compile(f'(TTJets.*|ZZ.*|ST.*|QCD_HT-mg.*|WW.*|WZ.*|W.*HT.*).*{year}')
#     region='cr_1m_j'
#     for distribution in ['ak4_eta0','ak4_btag','ak4_pt','ak4_eta','ak4_ptraw0','muon_pt','met','recoil','ak4_pt0']:
#         make_plot(copy.deepcopy(acc), region=region,distribution=distribution, year=year, data=data, mc=mc)
# for year in [2017]:
#     data = re.compile(f'EGamma.*{year}')
#     mc = re.compile(f'(GJets.*HT|QCD_HT-mg.*).*{year}')
#     region='cr_g_j'
#     for distribution in ['recoil','ak4_pt0','drphotonjet']:
#         make_plot(copy.deepcopy(acc), region=region,distribution=distribution, year=year, data=data, mc=mc)
