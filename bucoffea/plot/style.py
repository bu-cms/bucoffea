from matplotlib import pyplot as plt

def markers(tag):
    if tag =='data':
        ret = {
            'linestyle':'none',
            'marker': '.',
            'markersize': 10.,
            'color':'k',
            'elinewidth': 1,
            'emarker': '_'
        }
    return ret

def matplotlib_rc():
    plt.rc('mathtext',rm='Helvetica')
    plt.rc('mathtext',it='Helvetica')
    plt.rc('mathtext',bf='Helvetica')
    params = {'font.size':14, 'lines.linewidth' : 1}
    plt.rcParams.update(params)

from collections import defaultdict
def plot_settings():
    plot_settings = defaultdict(lambda: defaultdict(lambda : None),
    {
        'cr_2m_vbf' : {
            'recoil' : {
                'ylim' : (1e-3,1e3)
            },
            'dimuon_pt' : {
                'ylim' : (1e-3,1e4)
            },
            'ak4_pt0' : {
                'ylim' : (1e-3,1e3)
            },
            'ak4_pt1' : {
                'ylim' : (1e-3,1e3)
            },
            'ak4_pt' : {
                'ylim' : (1e-3,1e3)
            },
            'muon_pt' : {
                'ylim' : (1e-3,1e3)
            },
            'muon_pt0' : {
                'ylim' : (1e-3,1e3)
            },
            'muon_pt1' : {
                'ylim' : (1e-3,1e3)
            },
            'met' : {
                'ylim' : (1e-3,1e3)
            },
            'ak4_phi0' : {
                'ylim' : (1e0,1e6)
            },
            'ak4_phi' : {
                'ylim' : (1e0,1e7)
            },
            'muon_phi' : {
                'ylim' : (1e0,1e6)
            },
            'muon_phi0' : {
                'ylim' : (1e0,1e6)
            },
            'muon_phi1' : {
                'ylim' : (1e0,1e6)
            },
            'dimuon_mass' : {
                'ylim' : (1e-1,1e7)
            },
            'ak4_eta0' : {
                'xlim' : (-3,3),
                'ylim' : (1e0,1e6)
            },
            'ak4_eta1' : {
                'xlim' : (-3,3),
                'ylim' : (1e0,1e6)
            },
            'ak4_eta' : {
                'xlim' : (-3,3),
                'ylim' : (1e0,1e8)
            },
            'muon_eta' : {
                'xlim' : (-3,3),
                'ylim' : (1e0,1e7)
            },
            'muon_eta0' : {
                'xlim' : (-3,3),
                'ylim' : (1e-1,1e8)
            },
            'muon_eta1' : {
                'xlim' : (-3,3),
                'ylim' : (1e-1,1e8)
            },
            'dimuon_eta' : {
                'xlim' : (-3,3),
                'ylim' : (1e0,1e7)
            },
            'dimuon_dr' : {
                'xlim' : (0,2),
                'ylim' : (1e0,1e6)
            },
            'ak4_chf0' : {
                'xlim' : (0,1),
                'ylim' : (1e0,1e6)
            },
            'ak4_nhf0' : {
                'xlim' : (0,1),
                'ylim' : (1e0,1e6)
            },
            'drmuonjet' : {
                'xlim' : (0,1),
                'ylim' : (1e0,1e6)
            },
            'dpfcalo' : {
                'xlim' : (-0.75,0.75),
                'ylim' : (1e1,1e7)
            },
            'dphijr' : {
                'xlim' : (0,2),
                'ylim' : (1e0,1e6)
            },
            'dphijm' : {
                'xlim' : (0,2),
                'ylim' : (1e0,1e6)
            },
            'gen_dilepton_mult' : {
                'xlim' : (0,5),
                'ylim' : (0.1,1e7)
            },
            'ak4_mult' : {
                'xlim' : (0,10),
                'ylim' : (1e0,1e6)
            },
            'mjj' : {
                'ylim' : (1e-3,1e5)
            },
            'detajj' : {
                'ylim' : (1e-1,1e6)
            },
            'dphijj' : {
                'ylim' : (1e-1,1e6)
            }
        },
        'cr_2e_vbf' : {
            'recoil' : {
                'ylim' : (1e-3,1e3)
            },
            'recoil_nopu' : {
                'ylim' : (1e-3,1e3)
            },
            'recoil_nopog' : {
                'ylim' : (1e-3,1e3)
            },
            'recoil_nopref' : {
                'ylim' : (1e-3,1e3)
            },
            'dielectron_pt' : {
                'ylim' : (1e-3,1e4)
            },
            'ak4_pt0' : {
                'ylim' : (1e-3,1e3)
            },
            'ak4_pt1' : {
                'ylim' : (1e-3,1e3)
            },
            'ak4_pt' : {
                'ylim' : (1e-3,1e3)
            },
            'electron_pt' : {
                'ylim' : (1e-3,1e3)
            },
            'electron_pt0' : {
                'ylim' : (1e-3,1e3)
            },
            'electron_pt1' : {
                'ylim' : (1e-3,1e3)
            },
            'met' : {
                'ylim' : (1e-3,1e3)
            },
            'ak4_phi0' : {
                'ylim' : (1e0,1e6)
            },
            'ak4_phi' : {
                'ylim' : (1e0,1e6)
            },
            'electron_phi' : {
                'ylim' : (1e-1,1e6)
            },
            'electron_phi0' : {
                'ylim' : (1e1,1e5)
            },
            'electron_phi1' : {
                'ylim' : (1e1,1e5)
            },
            'electron_tightid1' : {
                'ylim' : None
            },
            'dielectron_mass' : {
                'ylim' : (1e-1,1e7)
            },
            'electron_dxy' : {
                'ylim' : (1e1,1e7),
                'xlim' : (0,0.15)
            },
            'electron_dz' : {
                'ylim' : (1e1,1e5)
            },
            'ak4_eta0' : {
                'xlim' : (-3,3),
                'ylim' : (1e0,1e6)
            },
            'ak4_eta1' : {
                'xlim' : (-3,3),
                'ylim' : (1e0,1e6)
            },
            'ak4_eta' : {
                'xlim' : (-3,3),
                'ylim' : (1e0,1e6)
            },
            'electron_eta' : {
                'xlim' : (-3,3),
                'ylim' : (1e-1,1e6)
            },
            'electron_eta0' : {
                'xlim' : (-3,3),
                'ylim' : (1e-1,1e6)
            },
            'electron_eta1' : {
                'xlim' : (-3,3),
                'ylim' : (1e-1,1e6)
            },
            'dielectron_eta' : {
                'xlim' : (-3,3),
                'ylim' : (1e-1,1e5)
            },
            'dielectron_dr' : {
                'xlim' : (0,2),
                'ylim' : (1e0,1e5)
            },
            'ak4_chf0' : {
                'xlim' : (0,1),
                'ylim' : (1e0,1e6)
            },
            'ak4_nhf0' : {
                'xlim' : (0,1),
                'ylim' : (1e0,1e6)
            },
            'drelejet' : {
                'xlim' : (0,1),
                'ylim' : (1e1,1e5)
            },
            'dpfcalo' : {
                'xlim' : (-0.75,0.75),
                'ylim' : (1e0,1e7)
            },
            'dphijr' : {
                'xlim' : (0,2),
                'ylim' : (1e0,1e6)
            },
            'dphijm' : {
                'xlim' : (0,2),
                'ylim' : (1e0,1e6)
            },
            'gen_dilepton_mult' : {
                'xlim' : (0,5),
                'ylim' : (0.1,1e7)
            },
            'ak4_mult' : {
                'xlim' : (0,10),
                'ylim' : (1e0,1e5)
            },
            'mjj' : {
                'ylim' : (1e-3,1e5)
            },
            'detajj' : {
                'ylim' : (1e-1,1e6)
            },
            'dphijj' : {
                'ylim' : (1e-1,1e6)
            }
        },
        'cr_1m_vbf' : {
            'recoil' : {
                'ylim' : (1e-3,1e5)
            },
            'ak4_pt0' : {
                'ylim' : (1e-3,1e6)
            },
            'ak4_pt1' : {
                'ylim' : (1e-3,1e6)
            },
            'ak4_pt' : {
                'ylim' : (1e-3,1e6)
            },
            'muon_pt' : {
                'ylim' : (1e-3,1e5)
            },
            'met' : {
                'ylim' : (1e-3,1e5)
            },
            'ak4_phi0' : {
                'ylim' : (1e0,1e7)
            },
            'ak4_phi' : {
                'ylim' : (1e0,1e8)
            },
            'muon_phi' : {
                'ylim' : (1e0,1e7)
            },
            'ak4_eta0' : {
                'xlim' : (-3,3),
                'ylim' : (1e0,1e9)
            },
            'ak4_eta1' : {
                'xlim' : (-3,3),
                'ylim' : (1e0,1e9)
            },
            'ak4_eta' : {
                'xlim' : (-3,3),
                'ylim' : (1e0,1e9)
            },
            'muon_eta' : {
                'xlim' : (-3,3),
                'ylim' : (1e0,1e7)
            },
            'ak4_chf0' : {
                'xlim' : (0,1),
                'ylim' : (1e0,1e8)
            },
            'ak4_nhf0' : {
                'xlim' : (0,1),
                'ylim' : (1e0,1e8)
            },
            'muon_mt' : {
                'xlim' : (0,180),
                'ylim' : (1e-1,1e4)
            },
            'drmuonjet' : {
                'xlim' : (0,1),
                'ylim' : (1e1,1e5)
            },
            'dpfcalo' : {
                'xlim' : (-0.75,0.75),
                'ylim' : (1e0,1e9)
            },
            'dphijr' : {
                'xlim' : (0,2),
                'ylim' : (1e0,1e8)
            },
            'dphijm' : {
                'xlim' : (0,2),
                'ylim' : (1e0,1e8)
            },
            'gen_dilepton_mult' : {
                'xlim' : (0,5),
                'ylim' : (0.1,1e7)
            },
            'ak4_mult' : {
                'xlim' : (0,10),
                'ylim' : (1e0,1e5)
            },
            'mjj' : {
                'ylim' : (1e-3,1e4)
            },
            'detajj' : {
                'ylim' : (1e-1,1e8)
            },
            'dphijj' : {
                'ylim' : (1e0,1e8)
            }
        },
        'cr_1e_vbf' : {
            'recoil' : {
                'ylim' : (1e-3,1e5)
            },
            'recoil' : {
                'ylim' : (1e-3,1e5)
            },
            'recoil_nopu' : {
                'ylim' : (1e-3,1e5)
            },
            'recoil_nopog' : {
                'ylim' : (1e-3,1e5)
            },
            'recoil_nopref' : {
                'ylim' : (1e-3,1e5)
            },
            'ak4_pt0' : {
                'ylim' : (1e-3,1e6)
            },
            'ak4_pt1' : {
                'ylim' : (1e-3,1e6)
            },
            'ak4_pt' : {
                'ylim' : (1e-3,1e6)
            },
            'electron_pt' : {
                'ylim' : (1e-3,1e5)
            },
            'met' : {
                'ylim' : (1e-3,1e5)
            },
            'ak4_phi0' : {
                'ylim' : (1e0,1e7)
            },
            'ak4_phi' : {
                'ylim' : (1e0,1e8)
            },
            'electron_phi' : {
                'ylim' : (1e0,1e7)
            },
            'electron_dxy' : {
                'ylim' : (1e1,1e7),
                'xlim' : (0,0.15)
            },
            'electron_dz' : {
                'ylim' : (1e1,1e5)
            },
            'ak4_eta0' : {
                'xlim' : (-3,3),
                'ylim' : (1e0,1e8)
            },
            'ak4_eta1' : {
                'xlim' : (-3,3),
                'ylim' : (1e0,1e8)
            },
            'ak4_eta' : {
                'xlim' : (-3,3),
                'ylim' : (1e0,1e8)
            },
            'electron_eta' : {
                'xlim' : (-3,3),
                'ylim' : (1e0,1e6)
            },
            'ak4_chf0' : {
                'xlim' : (0,1),
                'ylim' : (1e0,1e8)
            },
            'ak4_nhf0' : {
                'xlim' : (0,1),
                'ylim' : (1e0,1e8)
            },
            'electron_mt' : {
                'xlim' : (0,180),
                'ylim' : (1e-1,1e5)
            },
            'drelejet' : {
                'xlim' : (0,1),
                'ylim' : (1e0,1e5)
            },
            'dpfcalo' : {
                'xlim' : (-0.75,0.75),
                'ylim' : (1e0,1e8)
            },
            'dphijr' : {
                'xlim' : (0,2),
                'ylim' : (1e0,1e8)
            },
            'dphijm' : {
                'xlim' : (0,2),
                'ylim' : (1e0,1e8)
            },
            'gen_dilepton_mult' : {
                'xlim' : (0,5),
                'ylim' : (0.1,1e7)
            },
            'ak4_mult' : {
                'xlim' : (0,10),
                'ylim' : (1e0,1e5)
            },
            'mjj' : {
                'ylim' : (1e-3,1e5)
            },
            'detajj' : {
                'ylim' : (1e-1,1e8)
            },
            'dphijj' : {
                'ylim' : (1e0,1e7)
            }
        },
        'cr_g_vbf' : {
            'recoil' : {
                'ylim' : (1e-3,1e5)
            },
            'ak4_pt0' : {
                'ylim' : (1e-3,1e5)
            },
            'ak4_pt1' : {
                'ylim' : (1e-3,1e5)
            },
            'ak4_pt' : {
                'ylim' : (1e-3,1e5)
            },
            'photon_pt0' : {
                'ylim' : (1e-3,1e5)
            },
            'met' : {
                'ylim' : (1e-3,1e5)
            },
            'ak4_phi0' : {
                'ylim' : (1e0,1e5)
            },
            'ak4_phi' : {
                'ylim' : (1e0,1e5)
            },
            'photon_phi' : {
                'ylim' : (1e0,1e5)
            },
            'ak4_eta0' : {
                'xlim' : (-3,3),
                'ylim' : (1e0,1e7)
            },
            'ak4_eta1' : {
                'xlim' : (-3,3),
                'ylim' : (1e0,1e7)
            },
            'ak4_eta' : {
                'xlim' : (-3,3),
                'ylim' : (1e0,1e7)
            },
            'photon_eta' : {
                'xlim' : (-3,3),
                'ylim' : (1e0,1e6)
            },
            'ak4_chf0' : {
                'xlim' : (0,1),
                'ylim' : (1e0,1e8)
            },
            'ak4_nhf0' : {
                'xlim' : (0,1),
                'ylim' : (1e0,1e8)
            },
            'drphotonjet' : {
                'xlim' : (0,1),
                'ylim' : (1e0,1e5)
            },
            'dpfcalo' : {
                'xlim' : (-0.75,0.75),
                'ylim' : (1e0,1e7)
            },
            'dphijr' : {
                'xlim' : (0,2),
                'ylim' : (1e0,1e5)
            },
            'dphijm' : {
                'xlim' : (0,2),
                'ylim' : (1e0,1e5)
            },
            'ak4_mult' : {
                'xlim' : (0,10),
                'ylim' : (1e2,1e5)
            },
            'mjj' : {
                'ylim' : (1e-3,1e5)
            },
            'detajj' : {
                'ylim' : (1e-1,1e8)
            },
            'dphijj' : {
                'ylim' : (1e0,1e7)
            }
        },
        'cr_2m_j' : {
            'recoil' : {
                'ylim' : (1e-3,1e3)
            },
            'dimuon_pt' : {
                'ylim' : (1e-3,1e3)
            },
            'ak4_pt0' : {
                'ylim' : (1e-3,1e3)
            },
            'ak4_pt' : {
                'ylim' : (1e-3,1e3)
            },
            'muon_pt' : {
                'ylim' : (1e-3,1e3)
            },
            'muon_pt0' : {
                'ylim' : (1e-3,1e3)
            },
            'muon_pt1' : {
                'ylim' : (1e-3,1e3)
            },
            'met' : {
                'ylim' : (1e-3,1e3)
            },
            'ak4_phi0' : {
                'ylim' : (1e1,1e5)
            },
            'ak4_phi' : {
                'ylim' : (1e1,1e5)
            },
            'muon_phi' : {
                'ylim' : (1e1,1e5)
            },
            'muon_phi0' : {
                'ylim' : (1e1,1e5)
            },
            'muon_phi1' : {
                'ylim' : (1e1,1e5)
            },
            'dimuon_mass' : {
                'ylim' : (1e1,1e5)
            },
            'ak4_eta0' : {
                'xlim' : (-3,3),
                'ylim' : (1e3,1e5)
            },
            'ak4_eta' : {
                'xlim' : (-3,3),
                'ylim' : (1e3,1e5)
            },
            'muon_eta' : {
                'xlim' : (-3,3),
                'ylim' : (1e3,1e5)
            },
            'muon_eta0' : {
                'xlim' : (-3,3),
                'ylim' : (1e3,1e5)
            },
            'muon_eta1' : {
                'xlim' : (-3,3),
                'ylim' : (1e3,1e5)
            },
            'dimuon_eta' : {
                'xlim' : (-3,3),
                'ylim' : (1e3,1e5)
            },
            'dimuon_dr' : {
                'xlim' : (0,2),
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
            'drmuonjet' : {
                'xlim' : (0,1),
                'ylim' : (1e1,1e5)
            },
            'dpfcalo' : {
                'xlim' : (-0.75,0.75),
                'ylim' : (1e1,1e7)
            },
            'dphijr' : {
                'xlim' : (0,2),
                'ylim' : (1e1,1e5)
            },
            'dphijm' : {
                'xlim' : (0,2),
                'ylim' : (1e1,1e5)
            },
            'gen_dilepton_mult' : {
                'xlim' : (0,5),
                'ylim' : (0.1,1e7)
            },
            'ak4_mult' : {
                'xlim' : (0,10),
                'ylim' : (1e2,1e5)
            }
        },
        'cr_2e_j' : {
            'recoil' : {
                'ylim' : (1e-3,1e3)
            },
            'recoil_nopu' : {
                'ylim' : (1e-3,1e3)
            },
            'recoil_nopog' : {
                'ylim' : (1e-3,1e3)
            },
            'recoil_nopref' : {
                'ylim' : (1e-3,1e3)
            },
            'dielectron_pt' : {
                'ylim' : (1e-3,1e3)
            },
            'ak4_pt0' : {
                'ylim' : (1e-3,1e3)
            },
            'ak4_pt' : {
                'ylim' : (1e-3,1e3)
            },
            'electron_pt' : {
                'ylim' : (1e-3,1e3)
            },
            'electron_pt0' : {
                'ylim' : (1e-3,1e3)
            },
            'electron_pt1' : {
                'ylim' : (1e-3,1e3)
            },
            'met' : {
                'ylim' : (1e-3,1e3)
            },
            'ak4_phi0' : {
                'ylim' : (1e1,1e5)
            },
            'ak4_phi' : {
                'ylim' : (1e1,1e5)
            },
            'electron_phi' : {
                'ylim' : (1e1,1e5)
            },
            'electron_phi0' : {
                'ylim' : (1e1,1e5)
            },
            'electron_phi1' : {
                'ylim' : (1e1,1e5)
            },
            'electron_tightid1' : {
                'ylim' : None
            },
            'dielectron_mass' : {
                'ylim' : (1e1,1e5)
            },
            'electron_dxy' : {
                'ylim' : (1e1,1e7),
                'xlim' : (0,0.15)
            },
            'electron_dz' : {
                'ylim' : (1e1,1e5)
            },
            'ak4_eta0' : {
                'xlim' : (-3,3),
                'ylim' : (1e3,1e5)
            },
            'ak4_eta' : {
                'xlim' : (-3,3),
                'ylim' : (1e3,1e5)
            },
            'electron_eta' : {
                'xlim' : (-3,3),
                'ylim' : (1e3,1e5)
            },
            'electron_eta0' : {
                'xlim' : (-3,3),
                'ylim' : (1e3,1e5)
            },
            'electron_eta1' : {
                'xlim' : (-3,3),
                'ylim' : (1e3,1e5)
            },
            'dielectron_eta' : {
                'xlim' : (-3,3),
                'ylim' : (1e3,1e5)
            },
            'dielectron_dr' : {
                'xlim' : (0,2),
                'ylim' : (1e1,1e5)
            },
            'ak4_chf0' : {
                'xlim' : (0,1),
                'ylim' : (1e2,1e6)
            },
            'ak4_nhf0' : {
                'xlim' : (0,1),
                'ylim' : (1e2,1e6)
            },
            'drelejet' : {
                'xlim' : (0,1),
                'ylim' : (1e1,1e5)
            },
            'dpfcalo' : {
                'xlim' : (-0.75,0.75),
                'ylim' : (1e1,1e7)
            },
            'dphijr' : {
                'xlim' : (0,2),
                'ylim' : (1e1,1e5)
            },
            'dphijm' : {
                'xlim' : (0,2),
                'ylim' : (1e1,1e5)
            },
            'gen_dilepton_mult' : {
                'xlim' : (0,5),
                'ylim' : (0.1,1e7)
            },
            'ak4_mult' : {
                'xlim' : (0,10),
                'ylim' : (1e2,1e5)
            }
        },
        'cr_1m_j' : {
            'recoil' : {
                'ylim' : (1e-3,1e5)
            },
            'ak4_pt0' : {
                'ylim' : (1e-3,1e5)
            },
            'ak4_pt' : {
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
            'ak4_phi' : {
                'ylim' : (1e4,1e5)
            },
            'muon_phi' : {
                'ylim' : (1e4,1e5)
            },
            'ak4_eta0' : {
                'xlim' : (-3,3),
                'ylim' : (1e4,1e6)
            },
            'ak4_eta' : {
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
            'drmuonjet' : {
                'xlim' : (0,1),
                'ylim' : (1e1,1e5)
            },
            'dpfcalo' : {
                'xlim' : (-0.75,0.75),
                'ylim' : (1e1,1e7)
            },
            'dphijr' : {
                'xlim' : (0,2),
                'ylim' : (1e1,1e5)
            },
            'dphijm' : {
                'xlim' : (0,2),
                'ylim' : (1e1,1e5)
            },
            'gen_dilepton_mult' : {
                'xlim' : (0,5),
                'ylim' : (0.1,1e7)
            },
            'ak4_mult' : {
                'xlim' : (0,10),
                'ylim' : (1e2,1e5)
            }
        },
        'cr_1e_j' : {
            'recoil' : {
                'ylim' : (1e-3,1e5)
            },
            'recoil' : {
                'ylim' : (1e-3,1e5)
            },
            'recoil_nopu' : {
                'ylim' : (1e-3,1e5)
            },
            'recoil_nopog' : {
                'ylim' : (1e-3,1e5)
            },
            'recoil_nopref' : {
                'ylim' : (1e-3,1e5)
            },
            'ak4_pt0' : {
                'ylim' : (1e-3,1e5)
            },
            'ak4_pt' : {
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
            'ak4_phi' : {
                'ylim' : (1e4,1e5)
            },
            'electron_phi' : {
                'ylim' : (1e4,1e5)
            },
            'electron_dxy' : {
                'ylim' : (1e1,1e7),
                'xlim' : (0,0.15)
            },
            'electron_dz' : {
                'ylim' : (1e1,1e5)
            },
            'ak4_eta0' : {
                'xlim' : (-3,3),
                'ylim' : (1e4,1e6)
            },
            'ak4_eta' : {
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
            'drelejet' : {
                'xlim' : (0,1),
                'ylim' : (1e1,1e5)
            },
            'dpfcalo' : {
                'xlim' : (-0.75,0.75),
                'ylim' : (1e1,1e7)
            },
            'dphijr' : {
                'xlim' : (0,2),
                'ylim' : (1e1,1e5)
            },
            'dphijm' : {
                'xlim' : (0,2),
                'ylim' : (1e1,1e5)
            },
            'gen_dilepton_mult' : {
                'xlim' : (0,5),
                'ylim' : (0.1,1e7)
            },
            'ak4_mult' : {
                'xlim' : (0,10),
                'ylim' : (1e2,1e5)
            }
        },
        'cr_g_j' : {
            'recoil' : {
                'ylim' : (1e-3,1e5)
            },
            'ak4_pt0' : {
                'ylim' : (1e-3,1e5)
            },
            'ak4_pt' : {
                'ylim' : (1e-3,1e5)
            },
            'photon_pt0' : {
                'ylim' : (1e-3,1e5)
            },
            'met' : {
                'ylim' : (1e-3,1e5)
            },
            'ak4_phi0' : {
                'ylim' : (1e4,1e5)
            },
            'ak4_phi' : {
                'ylim' : (1e4,1e5)
            },
            'photon_phi' : {
                'ylim' : (1e4,1e5)
            },
            'ak4_eta0' : {
                'xlim' : (-3,3),
                'ylim' : (1e4,1e6)
            },
            'ak4_eta' : {
                'xlim' : (-3,3),
                'ylim' : (1e4,1e6)
            },
            'photon_eta' : {
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
            'drphotonjet' : {
                'xlim' : (0,1),
                'ylim' : (1e1,1e5)
            },
            'dpfcalo' : {
                'xlim' : (-0.75,0.75),
                'ylim' : (1e1,1e7)
            },
            'dphijr' : {
                'xlim' : (0,2),
                'ylim' : (1e1,1e5)
            },
            'dphijm' : {
                'xlim' : (0,2),
                'ylim' : (1e1,1e5)
            },
            'ak4_mult' : {
                'xlim' : (0,10),
                'ylim' : (1e2,1e5)
            }
        }
        }
    )
    return plot_settings
