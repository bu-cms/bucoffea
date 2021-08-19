#!/usr/bin/env python
import numpy as np
from coffea import hist

Bin = hist.Bin


common_distributions = [
    'mjj',
    'detajj',
    'dphijj',
    'recoil',
    'ak4_eta0',
    'ak4_eta1',
    'ak4_pt0',
    'ak4_pt1',
    'ak4_central_eta',
    'ak4_forward_eta',
    'dphijr',
    # 'ak4_phi0',
    # 'ak4_phi1',
    # 'ak4_nef0',
    # 'ak4_nef1',
    # 'ak4_nhf0',
    # 'ak4_nhf1',
    # 'ak4_chf0',
    # 'ak4_chf1',
    # 'ak4_mt0',
    # 'ak4_mt1',
    # 'vecb',
    # 'vecdphi',
    # 'dphitkpf',
    # 'met',
    # 'met_phi',
    # 'ak4_mult',
    # 'extra_ak4_mult',
    # 'calomet_pt',
    # 'calomet_phi',
]

# Distributions to plot for each region
distributions = {
    'sr_vbf' : common_distributions + ['ak4_nef0', 'ak4_nef1', 'ak4_nhf0', 'ak4_nhf1', 'ak4_chf0', 'ak4_chf1'],
    'cr_1m_vbf' : common_distributions + ['muon_pt', 'muon_eta', 'muon_phi', 'muon_mt'],
    'cr_1e_vbf' : common_distributions + ['electron_pt', 'electron_eta', 'electron_phi', 'electron_mt'],
    'cr_2m_vbf' : common_distributions + ['muon_pt0', 'muon_eta0', 'muon_phi0', 'muon_pt1', 'muon_eta1', 'muon_phi1', 'dimuon_mass'],
    'cr_2e_vbf' : common_distributions + ['electron_pt0', 'electron_eta0', 'electron_phi0', 'electron_pt1', 'electron_eta1', 'electron_phi1', 'dielectron_mass'],
    'cr_g_vbf'  : common_distributions + ['photon_pt0', 'photon_eta0', 'photon_phi0'],
} 

recoil_bins_2016 = [ 250,  280,  310,  340,  370,  400,  430,  470,  510, 550,  590,  640,  690,  740,  790,  840,  900,  960, 1020, 1090, 1160, 1250, 1400]

binnings = {
    'mjj': Bin('mjj', r'$M_{jj} \ (GeV)$', [200., 400., 600., 900., 1200., 1500., 2000., 2750., 3500., 5000.]),
    'ak4_pt0': Bin('jetpt',r'Leading AK4 jet $p_{T}$ (GeV)',list(range(80,600,20)) + list(range(600,1000,20)) ),
    'ak4_pt1': Bin('jetpt',r'Trailing AK4 jet $p_{T}$ (GeV)',list(range(40,600,20)) + list(range(600,1000,20)) ),
    'ak4_phi0' : Bin("jetphi", r"Leading AK4 jet $\phi$", 50,-np.pi, np.pi),
    'ak4_phi1' : Bin("jetphi", r"Trailing AK4 jet $\phi$", 50,-np.pi, np.pi),
    'ak4_central_eta' : Bin("jeteta", r"More Central Jet $\eta$", 50, -5, 5),
    'ak4_forward_eta' : Bin("jeteta", r"More Forward Jet $\eta$", 50, -5, 5),
    'ak4_nef0' : Bin('frac', 'Leading Jet Neutral EM Frac', 50, 0, 1),
    'ak4_nef1' : Bin('frac', 'Trailing Jet Neutral EM Frac', 50, 0, 1),
    'ak4_nhf0' : Bin('frac', 'Leading Jet Neutral Hadronic Frac', 50, 0, 1),
    'ak4_nhf1' : Bin('frac', 'Trailing Jet Neutral Hadronic Frac', 50, 0, 1),
    'ak4_chf0' : Bin('frac', 'Leading Jet Charged Hadronic Frac', 50, 0, 1),
    'ak4_chf1' : Bin('frac', 'Trailing Jet Charged Hadronic Frac', 50, 0, 1),
    # 'dphitkpf' : Bin('dphi', r'$\Delta\phi_{TK,PF}$', 50, 0, 3.5),
    'met' : Bin('met',r'$p_{T}^{miss}$ (GeV)',list(range(0,500,50)) + list(range(500,1000,100)) + list(range(1000,2000,250))),
    'met_phi' : Bin("phi", r"$\phi_{MET}$", 50, -np.pi, np.pi),
    'ak4_mult' : Bin("multiplicity", r"AK4 multiplicity", 10, -0.5, 9.5),
    'ak4_mt0' : Bin("mt", r"Leading AK4 $M_{T}$ (GeV)", 50, 0, 1000),
    'ak4_mt1' : Bin("mt", r"Trailing AK4 $M_{T}$ (GeV)", 50, 0, 1000),
    'dphijr' : Bin("dphi", r"$min\Delta\phi(j,recoil)$", 50, 0, 3.5),
    'extra_ak4_mult' : Bin("multiplicity", r"Additional AK4 Jet Multiplicity", 10, -0.5, 9.5),
    'recoil' : Bin('recoil','Recoil (GeV)', recoil_bins_2016),
}

ylims = {
    'ak4_eta0' : (1e-3,1e8),
    'ak4_eta1' : (1e-3,1e8),
    'ak4_nef0' : (1e0,1e8),
    'ak4_nef1' : (1e0,1e8),
    'ak4_nhf0' : (1e0,1e8),
    'ak4_nhf1' : (1e0,1e8),
    'ak4_chf0' : (1e0,1e8),
    'ak4_chf1' : (1e0,1e8),
    'vecb' : (1e-1,1e9),
    'vecdphi' : (1e0,1e9),
    'dphitkpf' : (1e0,1e9),
    'met' : (1e-3,1e5),
    'ak4_mult' : (1e-1,1e8),
}
