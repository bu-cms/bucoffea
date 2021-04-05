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
        'sr_vbf' : {
            'recoil' : {
                'ylim' : (1e-3,1e3)
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
            'met' : {
                'ylim' : (1e-6,1e3)
            },
            'puppimet' : {
                'ylim' : (1e-6,1e3)
            },
            'tkmet' : {
                'ylim' : (1e-6,1e3)
            },
            'calomet' : {
                'ylim' : (1e-6,1e3)
            },
            'ak4_phi0' : {
                'ylim' : (1e0,1e7)
            },
            'ak4_phi1' : {
                'ylim' : (1e0,1e7)
            },
            'ak4_phi' : {
                'ylim' : (1e0,1e7)
            },
            'ak4_eta0' : {
                'ylim' : (1e-3,1e8)
            },
            'ak4_eta1' : {
                'ylim' : (1e-3,1e8)
            },
            'hem_ht' : {
                'ylim' : (1e-1,1e5)
            },
            'hem_pt_max' : {
                'ylim' : (1e-1,1e5)
            },
            'photon_raw_hem_pt_sum' : {
                'xlim' : (0,300)
            },
            'photon_raw_anti_hem_pt_sum' : {
                'xlim' : (0,300)
            },
            'photon_raw_anti_hem_mult' : {
                'ylim' : (1e-3,1e4)
            },
            'photon_raw' : {
                'ylim' : (1e-3,1e4)
            },
            'photon_raw_hem_mult' : {
                'ylim' : (1e-3,1e4)
            },

            'ak4_eta' : {
                'ylim' : (1e0,1e8)
            },
            'ak4_chf0' : {
                'xlim' : (0,1),
                'ylim' : (1e1,1e8)
            },
            'ak4_nhf0' : {
                'xlim' : (0,1),
                'ylim' : (1e1,1e8)
            },
            'ak4_chf1' : {
                'xlim' : (0,1),
                'ylim' : (1e1,1e8)
            },
            'ak4_nhf1' : {
                'xlim' : (0,1),
                'ylim' : (1e1,1e8)
            },
            'dpfcalo' : {
                'xlim' : (-0.75,0.75),
                'ylim' : (1e1,1e7)
            },
            'dphijr' : {
                'xlim' : (0,3.2),
                'ylim' : (1e0,1e6)
            },
            'dphijm' : {
                'xlim' : (0,3.2),
                'ylim' : (1e0,1e6)
            },
            'ak4_mult' : {
                'xlim' : (0,10),
                'ylim' : (1e0,1e6)
            },
            'mjj' : {
                'ylim' : (1e-4,1e4)
            },
            'detajj' : {
                'ylim' : (1e-1,1e6)
            },
            'dphijj' : {
                'ylim' : (1e-1,1e6)
            },
            'recoil_phi' : {
                'ylim' : (1e-1,1e6)
            },
            'met_phi' : {
                'ylim' : (1e-1,1e6)
            },
            'ak4_sigma_eta_eta0' : {
                'ylim' : (1e-1,1e5)
            },
            'ak4_sigma_eta_eta1' : {
                'ylim' : (1e-1,1e5)
            },
            'ak4_sigma_phi_phi0' : {
                'ylim' : (1e-1,1e5)
            },
            'ak4_sigma_phi_phi1' : {
                'ylim' : (1e-1,1e5)
            },
            'dPFTkMET' : {
                'ylim' : (1e-1,1e5)
            },
            'vecb' : {
                'ylim' : (1e-1,1e8)
            },
            'vecdphi' : {
                'ylim' : (1e-1,1e8)
            },
            'dphitkpf' : {
                'ylim' : (1e-1,1e8)
            },
            'ak4_nef0' : {
                'ylim' : (1e-2,1e6)
            },
            'ak4_nef1' : {
                'ylim' : (1e-2,1e6)
            },
        },
        'sr_vbf_trk_ee' : {
            'ak4_nef0' : {
                'ylim' : (1e-2,1e6)
            },
            'ak4_nef1' : {
                'ylim' : (1e-2,1e6)
            },
            'vecb' : {
                'ylim' : (1e-1,1e8)
            },
            'vecdphi' : {
                'ylim' : (1e-1,1e8)
            },
            'dphitkpf' : {
                'ylim' : (1e-1,1e8)
            },
            'dPFTkMET' : {
                'ylim' : (1e-1,1e5)
            },

        },
        'sr_vbf_ee_ee' : {
            'ak4_nef0' : {
                'ylim' : (1e-2,1e6)
            },
            'ak4_nef1' : {
                'ylim' : (1e-2,1e6)
            },
            'vecb' : {
                'ylim' : (1e-1,1e8)
            },
            'vecdphi' : {
                'ylim' : (1e-1,1e8)
            },
            'dphitkpf' : {
                'ylim' : (1e-1,1e8)
            },
        },
        'sr_j' : {
            'recoil' : {
                'ylim' : (1e-3,1e4),
                'xlim' : (250,1550)
                },
            'recoil_phi' : {
                'ylim' : (1e3,1e6)
                },
            'ak4_eta0' : {
                'xlim' : (-2.5,2.5)
            },
            'dpfcalo' : {
                'xlim' : (-0.5,0.5)
            }
            },
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
                'ylim' : (1e0,1e7)
            },
            'ak4_phi1' : {
                'ylim' : (1e0,1e7)
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
                'ylim' : (1e0,1e6)
            },
            'ak4_eta1' : {
                'ylim' : (1e0,1e6)
            },
            'ak4_eta' : {
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
                'xlim' : (0,3.2),
                'ylim' : (1e0,1e6)
            },
            'dphijm' : {
                'xlim' : (0,3.2),
                'ylim' : (1e0,1e6)
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
                'ylim' : (1e0,1e7)
            },
            'ak4_phi1' : {
                'ylim' : (1e0,1e7)
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
                'ylim' : (1e0,1e6)
            },
            'ak4_eta1' : {
                'ylim' : (1e0,1e6)
            },
            'ak4_eta' : {
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
            'ak4_chf1' : {
                'xlim' : (0,1),
                'ylim' : (1e0,1e6)
            },
            'ak4_nhf1' : {
                'xlim' : (0,1),
                'ylim' : (1e0,1e6)
            },
            'drelejet' : {
                'xlim' : (0,2),
                'ylim' : (1e1,1e5)
            },
            'dpfcalo' : {
                'xlim' : (-0.75,0.75),
                'ylim' : (1e0,1e7)
            },
            'dphijr' : {
                'xlim' : (0,3.2),
                'ylim' : (1e0,1e6)
            },
            'dphijm' : {
                'xlim' : (0,3.2),
                'ylim' : (1e0,1e6)
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
            'ak4_phi1' : {
                'ylim' : (1e0,1e7)
            },
            'ak4_phi' : {
                'ylim' : (1e0,1e8)
            },
            'muon_phi' : {
                'ylim' : (1e0,1e7)
            },
            'ak4_eta0' : {
                'ylim' : (1e0,1e9)
            },
            'ak4_eta1' : {
                'ylim' : (1e0,1e9)
            },
            'ak4_eta' : {
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
                'xlim' : (0,3.2),
                'ylim' : (1e0,1e8)
            },
            'dphijm' : {
                'xlim' : (0,3.2),
                'ylim' : (1e0,1e8)
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
            'ak4_phi1' : {
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
                'ylim' : (1e0,1e8)
            },
            'ak4_eta1' : {
                'ylim' : (1e0,1e8)
            },
            'ak4_eta' : {
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
                'xlim' : (0,2),
                'ylim' : (1e0,1e5)
            },
            'dpfcalo' : {
                'xlim' : (-0.75,0.75),
                'ylim' : (1e0,1e8)
            },
            'dphijr' : {
                'xlim' : (0,3.2),
                'ylim' : (1e0,1e8)
            },
            'dphijm' : {
                'xlim' : (0,3.2),
                'ylim' : (1e0,1e8)
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
            'ak4_phi1' : {
                'ylim' : (1e0,1e7)
            },
            'ak4_phi' : {
                'ylim' : (1e0,1e5)
            },
            'photon_phi0' : {
                'ylim' : (1e0,1e5)
            },
            'ak4_eta0' : {
                'ylim' : (1e0,1e7)
            },
            'ak4_eta1' : {
                'ylim' : (1e0,1e7)
            },
            'ak4_eta' : {
                'ylim' : (1e0,1e7)
            },
            'photon_eta0' : {
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
                'xlim' : (0,3.2),
                'ylim' : (1e0,1e5)
            },
            'dphijm' : {
                'xlim' : (0,3.2),
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
            'ak4_pt0_over_recoil' : {
                'ylim' : (1e1,1e7)
            },
            'recoil' : {
                'ylim' : (1e-3,1e3)
            },
            'recoil_nopog' : {
                'ylim' : (1e-3,1e3)
            },
            'recoil_nopu' : {
                'ylim' : (1e-3,1e3)
            },
            'recoil_nopref' : {
                'ylim' : (1e-3,1e3)
            },
            'dimuon_pt' : {
                'ylim' : (1e-3,1e3)
            },
            'ak4_pt0' : {
                'ylim' : (1e-3,1e3)
            },
            'ak4_pt' : {
                'ylim' : (1e-3,1e4)
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
                'ylim' : (1e2,1e5)
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
                'ylim' : (1e1,1e6)
            },
            'ak4_chf0' : {
                'xlim' : (0,1),
                'ylim' : (1e2,1e6)
            },
            'ak4_nhf0' : {
                'xlim' : (0,1),
                'ylim' : (1e2,1e6)
            },
            'ak4_nef0' : {
                'xlim' : (0,1),
                'ylim' : (1e0,1e4)
            },
            'drmuonjet' : {
                'xlim' : (0,2),
                'ylim' : (1e1,1e5)
            },
            'dpfcalo' : {
                'xlim' : (-0.75,0.75),
                'ylim' : (1e1,1e7)
            },
            'dphijr' : {
                'xlim' : (0,3.2),
                'ylim' : (1e1,1e5)
            },
            'dphijm' : {
                'xlim' : (0,3.2),
                'ylim' : (1e1,1e5)
            },
            'ak4_mult' : {
                'xlim' : (0,10),
                'ylim' : (1e2,1e5)
            }
        },
        'cr_2e_j' : {
            'ak4_pt0_over_recoil' : {
                'ylim' : (1e1,1e7)
            },
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
                'ylim' : (1e2,1e5)
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
                'ylim' : (1e1,1e6)
            },
            'ak4_chf0' : {
                'xlim' : (0,1),
                'ylim' : (1e2,1e6)
            },
            'ak4_nhf0' : {
                'xlim' : (0,1),
                'ylim' : (1e2,1e6)
            },
            'ak4_nef0' : {
                'xlim' : (0,1),
                'ylim' : (1e0,1e4)
            },
            'drelejet' : {
                'xlim' : (0,2),
                'ylim' : (1e1,1e5)
            },
            'dpfcalo' : {
                'xlim' : (-0.75,0.75),
                'ylim' : (1e1,1e7)
            },
            'dphijr' : {
                'xlim' : (0,3.2),
                'ylim' : (1e1,1e5)
            },
            'dphijm' : {
                'xlim' : (0,3.2),
                'ylim' : (1e1,1e5)
            },
            'ak4_mult' : {
                'xlim' : (0,10),
                'ylim' : (1e2,1e5)
            }
        },
        'cr_1m_j' : {
            'ak4_pt0_over_recoil' : {
                'ylim' : (1e1,1e7)
            },
            'recoil' : {
                'ylim' : (1e-3,1e5)
            },
            'recoil_nopog' : {
                'ylim' : (1e-3,1e5)
            },
            'recoil_nopu' : {
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
            'muon_pt' : {
                'ylim' : (1e-3,1e5)
            },
            'met' : {
                'ylim' : (1e-3,1e5)
            },
            'ak4_phi0' : {
                'ylim' : (1e4,1e7)
            },
            'ak4_phi' : {
                'ylim' : (1e4,1e7)
            },
            'muon_phi' : {
                'ylim' : (1e4,1e5)
            },
            'ak4_eta0' : {
                'xlim' : (-3,3),
                'ylim' : (1e2,1e6)
            },
            'ak4_eta' : {
                'xlim' : (-3,3),
                'ylim' : (1e2,1e6)
            },
            'muon_eta' : {
                'xlim' : (-3,3),
                'ylim' : (1e4,1e6)
            },
            'ak4_chf0' : {
                'xlim' : (0,1),
                'ylim' : (1e3,1e8)
            },
            'ak4_nef0' : {
                'xlim' : (0,1),
                'ylim' : (1e0,1e5)
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
                'xlim' : (0,3.2),
                'ylim' : (1e3,1e7)
            },
            'dphijm' : {
                'xlim' : (0,3.2),
                'ylim' : (1e3,1e7)
            },
            'ak4_mult' : {
                'xlim' : (0,10),
                'ylim' : (1e2,1e5)
            }
        },
        'cr_1e_j' : {
            'ak4_pt0_over_recoil' : {
                'ylim' : (1e1,1e7)
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
                'ylim' : (1e4,1e7)
            },
            'electron_phi' : {
                'ylim' : (1e4,1e5)
            },
            'electron_dxy' : {
                'ylim' : (1e3,1e9),
                'xlim' : (0,0.15)
            },
            'electron_dz' : {
                'ylim' : (1e3,1e0),
                'xlim' : (0,0.25)
            },
            'ak4_eta0' : {
                'xlim' : (-3,3),
                'ylim' : (1e2,1e6)
            },
            'ak4_eta' : {
                'xlim' : (-3,3),
                'ylim' : (1e2,1e6)
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
            'ak4_nef0' : {
                'xlim' : (0,1),
                'ylim' : (1e0,1e5)
            },
            'electron_mt' : {
                'xlim' : (0,180),
                'ylim' : (1e1,1e5)
            },
            'drelejet' : {
                'xlim' : (0,2),
                'ylim' : (1e1,1e5)
            },
            'dpfcalo' : {
                'xlim' : (-0.75,0.75),
                'ylim' : (1e1,1e7)
            },
            'dphijr' : {
                'xlim' : (0,3.2),
                'ylim' : (1e1,1e5)
            },
            'dphijm' : {
                'xlim' : (0,3.2),
                'ylim' : (1e1,1e5)
            },
            'ak4_mult' : {
                'xlim' : (0,10),
                'ylim' : (1e2,1e5)
            }
        },
        'cr_g_j' : {
            'ak4_pt0_over_recoil' : {
                'ylim' : (1e1,1e7)
            },
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
            'photon_eta0' : {
                'xlim' : (-1.6, 1.6),
                'ylim' : (1e4,1e6)
            },
            'met' : {
                'ylim' : (1e-3,1e5)
            },
            'ak4_phi0' : {
                'ylim' : (1e4,1e5)
            },
            'ak4_phi' : {
                'ylim' : (1e4,1e7)
            },
            'photon_phi' : {
                'ylim' : (1e4,1e5)
            },
            'ak4_eta0' : {
                'xlim' : (-3,3),
                'ylim' : (1e2,1e6)
            },
            'ak4_eta' : {
                'xlim' : (-3,3),
                'ylim' : (1e2,1e6)
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
            'ak4_nef0' : {
                'xlim' : (0,1),
                'ylim' : (1e0,1e5)
            },
            'drphotonjet' : {
                'xlim' : (0,2),
                'ylim' : (1e1,1e5)
            },
            'dpfcalo' : {
                'xlim' : (-0.75,0.75),
                'ylim' : (1e1,1e7)
            },
            'dphijr' : {
                'xlim' : (0,3.2),
                'ylim' : (1e3,1e7)
            },
            'dphijm' : {
                'xlim' : (0,3.2),
                'ylim' : (1e3,1e7)
            },
            'ak4_mult' : {
                'xlim' : (0,10),
                'ylim' : (1e3,1e6)
            }
        },
        'cr_2m_v' : {
            'ak8_mass0' : {
                'xlim' : (60,110),
                'ylim' : (1e-1,1e5)
            },
            'ak8_tau210' : {
                'xlim' : (0,1),
                'ylim' : (1e-1,1e6)
            },
            'recoil' : {
                'xlim' : (250,1000),
                'ylim' : (1e-4,2e2)
            },
            'ak8_pt0' : {
                'xlim' : (200,1000),
                'ylim' : (1e-3,2e3)
            },
            'met' : {
                'xlim' : (0,1100),
                'ylim' : (1e-5,2e2)
            },
            'ak8_phi0' : {
                'ylim' : (1e1,1e5)
            },
            'ak8_eta0' : {
                'xlim' : (-3,3),
                'ylim' : (1e-1,1e5)
            }
        },
        'cr_2e_v' : {
            'ak8_mass0' : {
                'xlim' : (60,110),
                'ylim' : (1e-1,1e5)
            },
            'ak8_tau210' : {
                'xlim' : (0,1),
                'ylim' : (1e-1,1e6)
            },
            'recoil' : {
                'xlim' : (250,1000),
                'ylim' : (1e-4,2e2)
            },
            'ak8_pt0' : {
                'xlim' : (200,1000),
                'ylim' : (1e-3,2e3)
            },
            'met' : {
                'xlim' : (0,1100),
                'ylim' : (1e-5,2e2)
            },
            'ak8_phi0' : {
                'ylim' : (1e1,1e5)
            },
            'ak8_eta0' : {
                'xlim' : (-3,3),
                'ylim' : (1e-1,1e5)
            }
        },
        'cr_1m_v' : {
            'ak8_mass0' : {
                'xlim' : (60,110),
                'ylim' : (1e-1,1e5)
            },
            'ak8_tau210' : {
                'xlim' : (0,1),
                'ylim' : (1e-1,1e6)
            },
            'recoil' : {
                'xlim' : (250,1000),
                'ylim' : (1e-3,1e3)
            },
            'ak8_pt0' : {
                'xlim' : (200,1000),
                'ylim' : (1e-2,1e4)
            },
            'met' : {
                'ylim' : (1e-4,1e3)
            },
            'ak8_phi0' : {
                'ylim' : (1e1,1e5)
            },
            'ak8_eta0' : {
                'xlim' : (-3,3),
                'ylim' : (1e-1,1e5)
            }
        },
        'cr_1e_v' : {
            'ak8_mass0' : {
                'xlim' : (60,110),
                'ylim' : (1e-1,1e5)
            },
            'ak8_tau210' : {
                'xlim' : (0,1),
                'ylim' : (1e-1,1e6)
            },
            'recoil' : {
                'xlim' : (250,1000),
                'ylim' : (1e-3,1e3)
            },
            'ak8_pt0' : {
                'xlim' : (200,1000),
                'ylim' : (1e-2,1e4)
            },
            'met' : {
                'ylim' : (1e-4,1e3)
            },
            'ak8_phi0' : {
                'ylim' : (1e1,1e5)
            },
            'ak8_eta0' : {
                'xlim' : (-3,3),
                'ylim' : (1e-1,1e5)
            }
        },
        'cr_g_v' : {
            'ak8_mass0' : {
                'xlim' : (60,110),
                'ylim' : (1e-1,1e5)
            },
            'ak8_tau210' : {
                'xlim' : (0,1),
                'ylim' : (1e-1,1e6)
            },
            'recoil' : {
                'xlim' : (250,1000),
                'ylim' : (1e-3,1e3)
            },
            'ak8_pt0' : {
                'xlim' : (200,1000),
                'ylim' : (1e-2,1e4)
            },
            'met' : {
                'xlim' : (0,1100),
                'ylim' : (1e-4,1e3)
            },
            'ak8_phi0' : {
                'ylim' : (1e1,1e5)
            },
            'ak8_eta0' : {
                'xlim' : (-3,3),
                'ylim' : (1e-1,1e5)
            }
        },
        'cr_nobveto_v' : {
            'ak8_mass0' : {
                'xlim' : (60,110),
                'ylim' : (1e-1,1e5)
            },
            'ak8_tau210' : {
                'xlim' : (0,1),
                'ylim' : (1e-1,1e6)
            },
            'recoil' : {
                'xlim' : (250,1000),
                'ylim' : (1e-3,1e3)
            },
            'ak8_pt0' : {
                'xlim' : (200,1000),
                'ylim' : (1e-2,1e4)
            },
            'met' : {
                'xlim' : (200,2000),
                'ylim' : (1e-4,1e3)
            },
            'ak8_phi0' : {
                'ylim' : (1e1,1e5)
            },
            'ak8_eta0' : {
                'xlim' : (-3,3),
                'ylim' : (1e-1,1e5)
            },
            'ak8_wvsqcd0' : {
                'xlim' : (0,1)
            },
            'ak8_wvsqcdmd0' : {
                'xlim' : (0,1)
            }
        },
        'sr_v' : {
            'ak8_mass0' : {
                'xlim' : (60,110),
                'ylim' : (1e-1,1e5)
            },
            'ak8_tau210' : {
                'xlim' : (0,1),
                'ylim' : (1e-1,1e6)
            },
            'recoil' : {
                'xlim' : (250,1000),
                'ylim' : (1e-3,1e4)
            },
            'ak8_pt0' : {
                'xlim' : (200,1000),
                'ylim' : (1e-2,1e4)
            },
            'met' : {
                'xlim' : (200,2000),
                'ylim' : (1e-4,1e4)
            },
            'ak8_phi0' : {
                'ylim' : (1e1,1e5)
            },
            'ak8_eta0' : {
                'xlim' : (-3,3),
                'ylim' : (1e-1,1e5)
            },
            'ak8_wvsqcd0' : {
                'xlim' : (0,1)
            },
            'ak8_wvsqcdmd0' : {
                'xlim' : (0,1)
            }
        }
        }
    )
    # add axes limit for control regions with different working points
    for raw_region in ['sr_v','cr_2m_v','cr_1m_v','cr_1e_v','cr_2e_v','cr_g_v','cr_nobveto_v']:
        # add default axes limits for these regions:
        plot_settings[raw_region]['ak8_wvsqcd0']={
            'xlim' : (0,1)
        }
        plot_settings[raw_region]['ak8_wvsqcdmd0']={
            'xlim' : (0,1)
        }
        for wp in ['inclusive','loose','loosemd']:
            region = raw_region.replace('_v','_'+wp+'_v')
            plot_settings[region] = {
                'ak8_mass0' : {
                    'xlim' : (60,110),
                    'ylim' : (1e-1,1e5)
                },
                'ak8_tau210' : {
                    'xlim' : (0,1),
                    'ylim' : (1e-1,1e6)
                },
                'recoil' : {
                    'xlim' : (250,1000),
                    'ylim' : (1e-3,1e4)
                },
                'ak8_pt0' : {
                    'xlim' : (200,1000),
                    'ylim' : (1e-2,1e4)
                },
                'met' : {
                    'xlim' : (200,2000),
                    'ylim' : (1e-4,1e4)
                },
                'ak8_phi0' : {
                    'ylim' : (1e1,1e5)
                },
                'ak8_eta0' : {
                    'xlim' : (-3,3),
                    'ylim' : (1e-1,1e5)
                },
                'ak8_wvsqcd0' : {
                    'xlim' : (0,1)
                },
                'ak8_wvsqcdmd0' : {
                    'xlim' : (0,1)
                }
            }
            if wp == 'inclusive':
                plot_settings[region]['ak8_mass0']={
                    'ylim' : (1e-1,1e5)
                }
        for wp in ['tight','tightmd']:
            region = raw_region.replace('_v','_'+wp+'_v')
            plot_settings[region] = {
                'ak8_mass0' : {
                    'xlim' : (60,110),
                    'ylim' : (1e-3,1e3)
                },
                'ak8_tau210' : {
                    'xlim' : (0,1),
                    'ylim' : (1e-3,1e4)
                },
                'recoil' : {
                    'xlim' : (250,1000),
                    'ylim' : (1e-5,1e2)
                },
                'ak8_pt0' : {
                    'xlim' : (200,1000),
                    'ylim' : (1e-4,1e2)
                },
                'met' : {
                    'xlim' : (200,2000),
                    'ylim' : (1e-6,1e2)
                },
                'ak8_phi0' : {
                    'ylim' : (1e-1,1e3)
                },
                'ak8_eta0' : {
                    'xlim' : (-3,3),
                    'ylim' : (1e-2,1e4)
                },
                'ak8_wvsqcd0' : {
                    'xlim' : (0,1)
                },
                'ak8_wvsqcdmd0' : {
                    'xlim' : (0,1)
                }
            }
    return plot_settings

def tangocolor( palettenumber, colornumber):
    cols = [ [] for x in range( 3 ) ]
    cols[0] = [ '#729fcf',     # Sky Blue
               '#ef2929',     # Scarlet Red
               '#fcaf3e',     # Orange
               '#8ae234',     # Chameleon
               '#fce94f',     # Butter
               '#e9b96e',     # Chocolate
               '#ad7fa8',     # Plum
               '#eeeeec',     # Aluminium 1
               '#888a85' ]    # Aluminium 1

    cols[1] = [ '#3465a4',     # Sky Blue
               '#cc0000',     # Scarlet Red
               '#f57900',     # Orange
               '#73d216',     # Chameleon
               '#edd400',     # Butter
               '#c17d11',     # Chocolate
               '#75507b',     # Plum
               '#d3d7cf',     # Aluminium 1
               '#555753' ]    # Aluminium 1

    cols[2] = [ '#204a87',     # Sky Blue
               '#a40000',     # Scarlet Red
               '#ce5c00',     # Orange
               '#4e9a06',     # Chameleon
               '#c4a000',     # Butter
               '#8f5902',     # Chocolate
               '#5c3566',     # Plum
               '#babdb6',     # Aluminium 1
               '#2e3436' ]    # Aluminium 1
    try:
        return cols[palettenumber][colornumber]
    except:
        return '#000000'
