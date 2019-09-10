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

    # indir = "../input/21Aug19_v4_pdfwgt"

    # acc = acc_from_dir(indir)

    # for year in [2017]:
    #     data = re.compile(f'MET_{year}')
    #     mc = re.compile(f'(EW.*|TTJets.*|ZZ.*|ST.*|QCD_HT.*|WW.*|WZ.*|.*DYJetsToLL_M-50_HT_MLM.*){year}')
    #     region='cr_2m_j'
    #     for distribution in ['recoil']:
    #         make_plot(acc, region=region,distribution=distribution, year=year, data=data, mc=mc, ylim=(1e-3,1e3), outdir=f'./output/{os.path.basename(indir)}')
    
    
    indir = "./input/2019-09-05_all_new_sf_neweletrig"

    acc = acc_from_dir(indir)

    def dimuon_plots():
        for year in [2017,2018]:
            data = re.compile(f'MET_{year}')
            mc = re.compile(f'(EW.*|TTJets.*FXFX.*|ST.*|QCD_HT.*|Diboson.*|.*DYJetsToLL_M-50_HT_MLM.*){year}')
            region='cr_2m_j'
            for distribution in ['recoil','ak4_pt0','dimuon_mass','muon_pt', 'dimuon_pt','met']:
                make_plot(acc, region=region,distribution=distribution, year=year, data=data, mc=mc, ylim=(1e-3,1e3), outdir=f'./output/{os.path.basename(indir)}/{region}')
            for distribution in ['ak4_phi0','muon_phi','dimuon_mass']:
                make_plot(acc, region=region,distribution=distribution, year=year, data=data, mc=mc, ylim=(1e1,1e5), outdir=f'./output/{os.path.basename(indir)}/{region}')
            for distribution in ['ak4_eta0', 'muon_eta','dimuon_eta']:
                make_plot(acc, region=region,distribution=distribution, year=year, data=data, mc=mc, xlim=(-3,3),ylim=(1e3,1e5), outdir=f'./output/{os.path.basename(indir)}/{region}')
            for distribution in ['ak4_chf0','ak4_nhf0']:
                make_plot(acc, region=region,distribution=distribution, year=year, data=data, mc=mc, xlim=(0,1),ylim=(1e2,1e6), outdir=f'./output/{os.path.basename(indir)}/{region}')
    
    def dielectron_plots():
        for year in [2017,2018]:
            data = re.compile(f'EGamma_{year}')
            mc = re.compile(f'(EW.*|TTJets.*FXFX.*|ST.*|QCD_HT.*|Diboson.*|.*DYJetsToLL_M-50_HT_MLM.*){year}')
            region='cr_2e_j'
            for distribution in ['recoil','ak4_pt0','dielectron_mass','electron_pt', 'dielectron_pt','met']:
                make_plot(acc, region=region,distribution=distribution, year=year, data=data, mc=mc, ylim=(1e-3,1e3), outdir=f'./output/{os.path.basename(indir)}/{region}')
            for distribution in ['ak4_phi0','electron_phi','dielectron_mass']:
                make_plot(acc, region=region,distribution=distribution, year=year, data=data, mc=mc, ylim=(1e1,1e5), outdir=f'./output/{os.path.basename(indir)}/{region}')
            for distribution in ['ak4_eta0', 'electron_eta','dielectron_eta']:
                make_plot(acc, region=region,distribution=distribution, year=year, data=data, mc=mc, xlim=(-3,3),ylim=(1e3,1e5), outdir=f'./output/{os.path.basename(indir)}/{region}')
            for distribution in ['ak4_chf0','ak4_nhf0']:
                make_plot(acc, region=region,distribution=distribution, year=year, data=data, mc=mc, xlim=(0,1),ylim=(1e2,1e6), outdir=f'./output/{os.path.basename(indir)}/{region}')
    
    def single_muon_plots():
        for year in [2017, 2018]:
            data = re.compile(f'MET_{year}')
            mc = re.compile(f'(TTJets.*FXFX.*|ST.*|QCD_HT.*|Diboson.*|W.*HT.*).*{year}')
            region='cr_1m_j'
            for distribution in ['recoil','ak4_pt0','muon_pt','met']:
                make_plot(acc, region=region,distribution=distribution, year=year, data=data, mc=mc, ylim=(1e-3,1e5), outdir=f'./output/{os.path.basename(indir)}/{region}')
            for distribution in ['ak4_phi0','muon_phi']:
                make_plot(acc, region=region,distribution=distribution, year=year, data=data, mc=mc, ylim=(1e4,1e5), outdir=f'./output/{os.path.basename(indir)}/{region}')
            for distribution in ['ak4_eta0', 'muon_eta']:
                make_plot(acc, region=region,distribution=distribution, year=year, data=data, mc=mc, xlim=(-3,3),ylim=(1e4,1e6), outdir=f'./output/{os.path.basename(indir)}/{region}')
            for distribution in ['ak4_chf0','ak4_nhf0']:
                make_plot(acc, region=region,distribution=distribution, year=year, data=data, mc=mc, xlim=(0,1),ylim=(1e3,1e8), outdir=f'./output/{os.path.basename(indir)}/{region}')
            for distribution in ['muon_mt']:
                make_plot(acc, region=region,distribution=distribution, year=year, data=data, mc=mc, xlim=(0,180),ylim=(1e1,1e5), outdir=f'./output/{os.path.basename(indir)}/{region}')

    def single_electron_plots():
        for year in [2017, 2018]:
            data = re.compile(f'EGamma_{year}')
            mc = re.compile(f'(TTJets.*FXFX.*|ST.*|QCD_HT.*|Diboson.*|W.*HT.*).*{year}')
            region='cr_1e_j'
            for distribution in ['recoil','ak4_pt0','electron_pt','met']:
                make_plot(acc, region=region,distribution=distribution, year=year, data=data, mc=mc, ylim=(1e-3,1e5), outdir=f'./output/{os.path.basename(indir)}/{region}')
            for distribution in ['ak4_phi0','electron_phi']:
                make_plot(acc, region=region,distribution=distribution, year=year, data=data, mc=mc, ylim=(1e4,1e5), outdir=f'./output/{os.path.basename(indir)}/{region}')
            for distribution in ['ak4_eta0', 'electron_eta']:
                make_plot(acc, region=region,distribution=distribution, year=year, data=data, mc=mc, xlim=(-3,3),ylim=(1e4,1e6), outdir=f'./output/{os.path.basename(indir)}/{region}')
            for distribution in ['ak4_chf0','ak4_nhf0']:
                make_plot(acc, region=region,distribution=distribution, year=year, data=data, mc=mc, xlim=(0,1),ylim=(1e3,1e8), outdir=f'./output/{os.path.basename(indir)}/{region}')
            for distribution in ['electron_mt']:
                make_plot(acc, region=region,distribution=distribution, year=year, data=data, mc=mc, xlim=(0,180),ylim=(1e1,1e5), outdir=f'./output/{os.path.basename(indir)}/{region}')

    def photon_plots():
        for year in [2017, 2018]:
            data = re.compile(f'EGamma_{year}')
            mc = re.compile(f'(GJets.*|QCD_HT.*|W.*HT.*).*{year}')
            region='cr_g_j'
            for distribution in ['recoil','ak4_pt0','photon_pt0','met']:
                make_plot(acc, region=region,distribution=distribution, year=year, data=data, mc=mc, ylim=(1e-3,1e5), outdir=f'./output/{os.path.basename(indir)}/{region}')
            for distribution in ['ak4_phi0','photon_phi0']:
                make_plot(acc, region=region,distribution=distribution, year=year, data=data, mc=mc, ylim=(4e4,5e5), outdir=f'./output/{os.path.basename(indir)}/{region}')
            for distribution in ['ak4_eta0', 'photon_eta0']:
                make_plot(acc, region=region,distribution=distribution, year=year, data=data, mc=mc, xlim=(-3,3),ylim=(1e4,1e6), outdir=f'./output/{os.path.basename(indir)}/{region}')
            for distribution in ['ak4_chf0','ak4_nhf0']:
                make_plot(acc, region=region,distribution=distribution, year=year, data=data, mc=mc, xlim=(0,1),ylim=(1e3,1e8), outdir=f'./output/{os.path.basename(indir)}/{region}')

    # single_muon_plots()
    # dimuon_plots()
    # single_electron_plots()
    photon_plots()
    dielectron_plots()
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
