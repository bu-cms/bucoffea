#!/usr/bin/env python
import numpy as np
from coffea import hist

Bin = hist.Bin

distributions = [
    'mjj',
    'ak4_eta0',
    'ak4_eta1',
    'ak4_phi0',
    'ak4_phi1',
    'ak4_pt0',
    'ak4_pt1',
    'ak4_nef0',
    'ak4_nef1',
    'ak4_nhf0',
    'ak4_nhf1',
    'ak4_chf0',
    'ak4_chf1',
    # 'ak4_mt0',
    # 'ak4_mt1',
    'vecb',
    'vecdphi',
    'dphitkpf',
    'met',
    'met_phi',
    'ak4_mult',
]

binnings = {
    'mjj': Bin('mjj', r'$M_{jj} \ (GeV)$', [200., 400., 600., 900., 1200., 1500., 2000., 2750., 3500., 5000.]),
    'ak4_pt0': Bin('jetpt',r'Leading AK4 jet $p_{T}$ (GeV)',list(range(80,600,20)) + list(range(600,1000,20)) ),
    'ak4_pt1': Bin('jetpt',r'Trailing AK4 jet $p_{T}$ (GeV)',list(range(40,600,20)) + list(range(600,1000,20)) ),
    'ak4_phi0' : Bin("jetphi", r"Leading AK4 jet $\phi$", 50,-np.pi, np.pi),
    'ak4_phi1' : Bin("jetphi", r"Trailing AK4 jet $\phi$", 50,-np.pi, np.pi),
    'ak4_nef0' : Bin('frac', 'Leading Jet Neutral EM Frac', 50, 0, 1),
    'ak4_nef1' : Bin('frac', 'Trailing Jet Neutral EM Frac', 50, 0, 1),
    'ak4_nhf0' : Bin('frac', 'Leading Jet Neutral Hadronic Frac', 50, 0, 1),
    'ak4_nhf1' : Bin('frac', 'Trailing Jet Neutral Hadronic Frac', 50, 0, 1),
    'ak4_chf0' : Bin('frac', 'Leading Jet Charged Hadronic Frac', 50, 0, 1),
    'ak4_chf1' : Bin('frac', 'Trailing Jet Charged Hadronic Frac', 50, 0, 1),
    'dphitkpf' : Bin('dphi', r'$\Delta\phi_{TK,PF}$', 50, 0, 3.5),
    'met' : Bin('met',r'$p_{T}^{miss}$ (GeV)',list(range(0,500,50)) + list(range(500,1000,100)) + list(range(1000,2000,250))),
    'met_phi' : Bin("phi", r"$\phi_{MET}$", 50, -np.pi, np.pi),
    'ak4_mult' : Bin("multiplicity", r"AK4 multiplicity", 10, -0.5, 9.5),
    'ak4_mt0' : Bin("mt", r"Leading AK4 $M_{T}$ (GeV)", 50, 0, 1000),
    'ak4_mt1' : Bin("mt", r"Trailing AK4 $M_{T}$ (GeV)", 50, 0, 1000),
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
