#!/usr/bin/env python
from coffea import hist

Bin = hist.Bin

distributions = [
    'mjj',
    'ak4_eta0',
    'ak4_eta1',
    'ak4_pt0',
    'ak4_pt1',
    'ak4_nef0',
    'ak4_nef1',
    'vecb',
    'vecdphi',
    'dphitkpf',
    'met'
]

binnings = {
    'mjj': Bin('mjj', r'$M_{jj} \ (GeV)$', [200., 400., 600., 900., 1200., 1500., 2000., 2750., 3500.]),
    'ak4_pt0': Bin('jetpt',r'Leading AK4 jet $p_{T}$ (GeV)',list(range(80,600,20)) + list(range(600,1000,20)) ),
    'ak4_pt1': Bin('jetpt',r'Trailing AK4 jet $p_{T}$ (GeV)',list(range(80,600,20)) + list(range(600,1000,20)) ),
    'ak4_nef0' : Bin('frac', 'Leading Jet Neutral EM Frac', 50, 0, 1),
    'ak4_nef1' : Bin('frac', 'Trailing Jet Neutral EM Frac', 50, 0, 1),
    'dphitkpf' : Bin('dphi', r'$\Delta\phi_{TK,PF}$', 50, 0, 3.5),
    'met' : Bin('met',r'$p_{T}^{miss}$ (GeV)',list(range(0,500,50)) + list(range(500,1000,100)) + list(range(1000,2000,250)))
}

ylims = {
    'ak4_nef0' : (1e0,1e8),
    'ak4_nef1' : (1e0,1e8),
    'vecb' : (1e-1,1e9),
    'vecdphi' : (1e0,1e9),
    'dphitkpf' : (1e0,1e9),
    'met' : (1e-3,1e5)
}