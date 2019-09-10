#!/usr/bin/env python

import os
import re
from pprint import pprint
from bucoffea.plot.util import acc_from_dir
from bucoffea.plot.stack_plot import Style, make_plot
from bucoffea.plot.cr_ratio_plot import cr_ratio_plot

from collections import defaultdict
plot_settings = defaultdict(lambda: defaultdict(lambda : None),
{
    'cr_2m_j' : {
        'recoil' : {
            'ylim' : (1e-3,1e3)
        },
        'ak4_pt0' : {
            'ylim' : (1e-3,1e3)
        },
        'muon_pt' : {
            'ylim' : (1e-3,1e3)
        },
        'met' : {
            'ylim' : (1e-3,1e3)
        },
        'ak4_phi0' : {
            'ylim' : (1e1,1e5)
        },
        'muon_phi' : {
            'ylim' : (1e1,1e5)
        },
        'dimuon_mass' : {
            'ylim' : (1e1,1e5)
        },
        'ak4_eta0' : {
            'xlim' : (-3,3),
            'ylim' : (1e3,1e5)
        },
        'muon_eta' : {
            'xlim' : (-3,3),
            'ylim' : (1e3,1e5)
        },
        'dimuon_eta' : {
            'xlim' : (-3,3),
            'ylim' : (1e3,1e5)
        },
        'ak4_chf0' : {
            'xlim' : (0,1),
            'ylim' : (1e2,1e6)
        },
        'ak4_nhf0' : {
            'xlim' : (0,1),
            'ylim' : (1e2,1e6)
        },
    },
    'cr_2e_j' : {
        'recoil' : {
            'ylim' : (1e-3,1e3)
        },
        'ak4_pt0' : {
            'ylim' : (1e-3,1e3)
        },
        'electron_pt' : {
            'ylim' : (1e-3,1e3)
        },
        'met' : {
            'ylim' : (1e-3,1e3)
        },
        'ak4_phi0' : {
            'ylim' : (1e1,1e5)
        },
        'electron_phi' : {
            'ylim' : (1e1,1e5)
        },
        'dielectron_mass' : {
            'ylim' : (1e1,1e5)
        },
        'ak4_eta0' : {
            'xlim' : (-3,3),
            'ylim' : (1e3,1e5)
        },
        'electron_eta' : {
            'xlim' : (-3,3),
            'ylim' : (1e3,1e5)
        },
        'dielectron_eta' : {
            'xlim' : (-3,3),
            'ylim' : (1e3,1e5)
        },
        'ak4_chf0' : {
            'xlim' : (0,1),
            'ylim' : (1e2,1e6)
        },
        'ak4_nhf0' : {
            'xlim' : (0,1),
            'ylim' : (1e2,1e6)
        },
    },
    'cr_1m_j' : {
        'recoil' : {
            'ylim' : (1e-3,1e5)
        },
        'ak4_pt0' : {
            'ylim' : (1e-3,1e5)
        },
        'muon_pt' : {
            'ylim' : (1e-3,1e5)
        },
        'met' : {
            'ylim' : (1e-3,1e5)
        },
        'ak4_phi0' : {
            'ylim' : (1e4,1e5)
        },
        'muon_phi' : {
            'ylim' : (1e4,1e5)
        },
        'ak4_eta0' : {
            'xlim' : (-3,3),
            'ylim' : (1e4,1e6)
        },
        'muon_eta' : {
            'xlim' : (-3,3),
            'ylim' : (1e4,1e6)
        },
        'ak4_chf0' : {
            'xlim' : (0,1),
            'ylim' : (1e3,1e8)
        },
        'ak4_nhf0' : {
            'xlim' : (0,1),
            'ylim' : (1e3,1e8)
        },
        'muon_mt' : {
            'xlim' : (0,180),
            'ylim' : (1e1,1e5)
        },
    },
    'cr_1e_j' : {
        'recoil' : {
            'ylim' : (1e-3,1e5)
        },
        'ak4_pt0' : {
            'ylim' : (1e-3,1e5)
        },
        'electron_pt' : {
            'ylim' : (1e-3,1e5)
        },
        'met' : {
            'ylim' : (1e-3,1e5)
        },
        'ak4_phi0' : {
            'ylim' : (1e4,1e5)
        },
        'electron_phi' : {
            'ylim' : (1e4,1e5)
        },
        'ak4_eta0' : {
            'xlim' : (-3,3),
            'ylim' : (1e4,1e6)
        },
        'electron_eta' : {
            'xlim' : (-3,3),
            'ylim' : (1e4,1e6)
        },
        'ak4_chf0' : {
            'xlim' : (0,1),
            'ylim' : (1e3,1e8)
        },
        'ak4_nhf0' : {
            'xlim' : (0,1),
            'ylim' : (1e3,1e8)
        },
        'electron_mt' : {
            'xlim' : (0,180),
            'ylim' : (1e1,1e5)
        },
    }
    }
)
def cr_ratio():
        indir=os.path.abspath('input/2019-09-09_gen_dilep_sf/')
        acc = acc_from_dir(indir)
        

        
        for year in [2017,2018]:
            data = {
                'cr_1m_j' : f'MET_{year}',
                'cr_2m_j' : f'MET_{year}',
                'cr_1e_j' : f'EGamma_{year}',
                'cr_2e_j' : f'EGamma_{year}',
                'cr_g_j' : f'EGamma_{year}',
            }
            mc_lo = {
                'cr_1m_j' : re.compile(f'(TTJets.*FXFX.*|Diboson.*|ST.*|QCD_HT.*|.*DYJetsToLL_M-50_HT_MLM.*|W.*HT.*).*{year}'),
                'cr_1e_j' : re.compile(f'(TTJets.*FXFX.*|Diboson.*|ST.*|QCD_HT.*|.*DYJetsToLL_M-50_HT_MLM.*|W.*HT.*).*{year}'),
                'cr_2m_j' : re.compile(f'(EW.*|TTJets.*FXFX.*|Diboson.*|ST.*|QCD_HT.*|.*DYJetsToLL_M-50_HT_MLM.*).*{year}'),
                'cr_2e_j' : re.compile(f'(EW.*|TTJets.*FXFX.*|Diboson.*|ST.*|QCD_HT.*|.*DYJetsToLL_M-50_HT_MLM.*).*{year}'),
                'cr_g_j' : re.compile(f'(GJets.*|QCD_HT.*|W.*HT.*).*{year}'),
            }
            mc_nlo = {
                    'cr_1m_j' : re.compile(f'(TTJets.*FXFX.*|Diboson.*|ST.*|QCD_HT.*|.*DY.*FXFX.*|W.*FXFX.*).*{year}'),
                    'cr_1e_j' : re.compile(f'(TTJets.*FXFX.*|Diboson.*|ST.*|QCD_HT.*|.*DY.*FXFX.*|W.*FXFX.*).*{year}'),
                    'cr_2m_j' : re.compile(f'(EW.*|TTJets.*FXFX.*|Diboson.*|ST.*|QCD_HT.*|.*DY.*FXFX.*).*{year}'),
                    'cr_2e_j' : re.compile(f'(EW.*|TTJets.*FXFX.*|Diboson.*|ST.*|QCD_HT.*|.*DY.*FXFX.*).*{year}'),
                    'cr_g_j' : re.compile(f'(GJets.*|QCD_HT.*|W.*FXFX.*).*{year}'),
            }
            outdir = f'./output/{os.path.basename(indir)}/ratios'
            cr_ratio_plot(acc, year=year,tag='losf',outdir=outdir, mc=mc_lo)
            cr_ratio_plot(acc, year=year,tag='nlo',outdir=outdir, mc=mc_nlo)


            for region in mc_lo.keys():
                outdir = f'./output/{os.path.basename(indir)}/{region}'
                plotset = plot_settings[region]
                for distribution in plotset.keys():
                    make_plot(acc, 
                            region=region,
                            distribution=distribution, 
                            year=year, 
                            data=data[region], 
                            mc=mc_lo[region], 
                            ylim=plot_settings['ylim'], 
                            tag = 'losf',
                            outdir=f'./output/{os.path.basename(indir)}/{region}')
                    make_plot(acc, 
                            region=region,
                            distribution=distribution, 
                            year=year, 
                            data=data[region], 
                            mc=mc_nlo[region], 
                            ylim=plot_settings['ylim'], 
                            tag = 'nlo',
                            outdir=f'./output/{os.path.basename(indir)}/{region}')


def main():
    cr_ratio()
    

if __name__ == "__main__":
    main()
