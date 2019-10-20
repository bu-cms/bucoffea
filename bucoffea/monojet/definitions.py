import copy
from coffea import hist

Hist = hist.Hist
Bin = hist.Bin
Cat = hist.Cat

from coffea.analysis_objects import JaggedCandidateArray
import coffea.processor as processor
from awkward import JaggedArray
import numpy as np
from bucoffea.helpers import object_overlap
from bucoffea.helpers.paths import bucoffea_path
from bucoffea.helpers.gen import find_first_parent
from bucoffea.helpers.dataset import extract_year
from pprint import pprint

def monojet_accumulator(cfg):
    dataset_ax = Cat("dataset", "Primary dataset")
    region_ax = Cat("region", "Selection region")
    type_ax = Cat("type", "Type")

    vpt_ax = Bin("vpt",r"$p_{T}^{V}$ (GeV)", 100, 0, 2000)

    met_ax = Bin("met", r"$p_{T}^{miss}$ (GeV)", 200, 0, 2000)
    recoil_ax = Bin("recoil", r"Recoil (GeV)", 200, 0, 2000)

    jet_pt_ax = Bin("jetpt", r"$p_{T}$ (GeV)", 100, 0, 1000)
    jet_pt_ax_coarse = Bin("jetpt", r"$p_{T}$ (GeV)", 5, 0, 500)
    jet_eta_ax = Bin("jeteta", r"$\eta$", 50, -5, 5)
    jet_eta_ax_coarse = Bin("jeteta", r"$\eta$", 10, -5, 5)
    jet_phi_ax = Bin("jetphi", r"$\phi$", 50,-np.pi, np.pi)

    jet_mass_ax = Bin("mass", r"$M_{jet}$ (GeV)", 100,0,300)

    dpfcalo_ax = Bin("dpfcalo", r"$(CaloMET-PFMET) / Recoil$", 20, -1, 1)
    btag_ax = Bin("btag", r"B tag discriminator", 20, 0, 1)
    multiplicity_ax = Bin("multiplicity", r"multiplicity", 10, -0.5, 9.5)
    nconst_ax = Bin("nconst", r"Number of constituents", 25, -0.5, 99.5)
    dphi_ax = Bin("dphi", r"$\Delta\phi$", 50, 0, 3.5)
    dr_ax = Bin("dr", r"$\Delta R$", 50, 0, 2)

    dxy_ax = Bin("dxy", r"$d_{xy}$", 20, 0, 0.5)
    dz_ax = Bin("dz", r"$d_{z}$", 20, 0, 0.5)
    id_ax = Bin("id", r"ID bool", 2,-0.5,1.5)

    pt_ax = Bin("pt", r"$p_{T}$ (GeV)", 100, 0, 1000)
    ht_ax = Bin("ht", r"$H_{T}$ (GeV)", 100, 0, 4000)
    mt_ax = Bin("mt", r"$M_{T}$ (GeV)", 100, 0, 1000)
    eta_ax = Bin("eta", r"$\eta$", 50, -5, 5)
    phi_ax = Bin("phi", r"$\phi$", 50,-np.pi, np.pi)
    
    ratio_ax = Bin("ratio", "ratio", 50,0,2)

    tau21_ax = Bin("tau21", r"Tagger", 50,-5,5)
    tagger_ax = Bin("tagger", r"Tagger", 50,-5,5)

    dilepton_mass_ax = Bin("dilepton_mass", r"$M(\ell\ell)$ (GeV)", 100,50,150)

    weight_type_ax = Cat("weight_type", "Weight type")
    weight_ax = Bin("weight_value", "Weight",100,0.5,1.5)
    weight_wide_ax = Bin("weight_value", "Weight",100,-10,10)

    nvtx_ax = Bin('nvtx','Number of vertices',100,-0.5,99.5)
    rho_ax = Bin('rho','Energy density',100, 0, 100)
    frac_ax = Bin('frac','Fraction', 50, 0, 1)
    Hist = hist.Hist
    items = {}
    items["genvpt_check"] = Hist("Counts", dataset_ax, type_ax, vpt_ax)
    items["lhe_njets"] = Hist("Counts", dataset_ax, multiplicity_ax)
    items["lhe_ht"] = Hist("Counts", dataset_ax, ht_ax)
    items["lhe_htinc"] = Hist("Counts", dataset_ax, ht_ax)
    items["met"] = Hist("Counts", dataset_ax, region_ax, met_ax)
    items["met_phi"] = Hist("Counts", dataset_ax, region_ax, phi_ax)
    items["met_noweight"] = Hist("Counts", dataset_ax, region_ax, met_ax)
    items["recoil"] = Hist("Counts", dataset_ax, region_ax, recoil_ax)
    items["recoil_nopog"] = Hist("Counts", dataset_ax, region_ax, recoil_ax)
    items["recoil_nopu"] = Hist("Counts", dataset_ax, region_ax, recoil_ax)
    items["recoil_notrg"] = Hist("Counts", dataset_ax, region_ax, recoil_ax)
    items["recoil_nopref"] = Hist("Counts", dataset_ax, region_ax, recoil_ax)
    items["recoil_phi"] = Hist("Counts", dataset_ax, region_ax, phi_ax)
    items["recoil_noweight"] = Hist("Counts", dataset_ax, region_ax, recoil_ax)
    items["ak4_pt0_over_recoil"] = Hist("Counts", dataset_ax, region_ax, ratio_ax)
    items["ak4_pt0"] = Hist("Counts", dataset_ax, region_ax, jet_pt_ax)
    items["ak4_ptraw0"] = Hist("Counts", dataset_ax, region_ax, jet_pt_ax)
    items["ak4_eta0"] = Hist("Counts", dataset_ax, region_ax, jet_eta_ax)
    items["ak4_phi0"] = Hist("Counts", dataset_ax, region_ax, jet_phi_ax)
    items["ak4_chf0"] = Hist("Counts", dataset_ax, region_ax, frac_ax)
    items["ak4_nhf0"] = Hist("Counts", dataset_ax, region_ax, frac_ax)
    items["ak4_nconst0"] = Hist("Counts", dataset_ax, region_ax, nconst_ax)

    items["ak4_pt0_chf0"] = Hist("Counts", dataset_ax, region_ax, jet_pt_ax_coarse, frac_ax)
    items["ak4_pt0_nhf0"] = Hist("Counts", dataset_ax, region_ax, jet_pt_ax_coarse, frac_ax)
    items["ak4_pt0_nconst0"] = Hist("Counts", dataset_ax, region_ax, jet_pt_ax_coarse, nconst_ax)
    items["ak4_pt0_eta0"] = Hist("Counts", dataset_ax, region_ax, jet_pt_ax,jet_eta_ax_coarse)

    items["ak4_pt"] = Hist("Counts", dataset_ax, region_ax, jet_pt_ax)
    items["ak4_eta"] = Hist("Counts", dataset_ax, region_ax, jet_eta_ax)
    items["ak4_phi"] = Hist("Counts", dataset_ax, region_ax, jet_phi_ax)
    items["ak4_eta_phi"] = Hist("Counts", dataset_ax, region_ax, eta_ax, phi_ax)
    items["ak4_pt_nopref"] = Hist("Counts", dataset_ax, region_ax, jet_pt_ax)
    items["ak4_eta_nopref"] = Hist("Counts", dataset_ax, region_ax, jet_eta_ax)
    items["ak4_phi_nopref"] = Hist("Counts", dataset_ax, region_ax, jet_phi_ax)
    items["ak4_btag"] = Hist("Counts", dataset_ax, region_ax, btag_ax)

    items["ak8_pt0"] = Hist("Counts", dataset_ax, region_ax, jet_pt_ax)
    items["ak8_eta0"] = Hist("Counts", dataset_ax, region_ax, jet_eta_ax)
    items["ak8_phi0"] = Hist("Counts", dataset_ax, region_ax, jet_phi_ax)
    items["ak8_mass0"] = Hist("Counts", dataset_ax, region_ax, jet_mass_ax)
    items["ak8_tau210"] = Hist("Counts", dataset_ax, region_ax, tau21_ax)
    items["ak8_wvsqcd0"] = Hist("Counts", dataset_ax, region_ax, tagger_ax)
    items["ak8_wvsqcdmd0"] = Hist("Counts", dataset_ax, region_ax, tagger_ax)
    items["ak8_zvsqcd0"] = Hist("Counts", dataset_ax, region_ax, tagger_ax)
    items["ak8_zvsqcdmd0"] = Hist("Counts", dataset_ax, region_ax, tagger_ax)

    items["ak8_pt"] = Hist("Counts", dataset_ax, region_ax, jet_pt_ax)
    items["ak8_eta"] = Hist("Counts", dataset_ax, region_ax, jet_eta_ax)
    items["ak8_phi"] = Hist("Counts", dataset_ax, region_ax, jet_phi_ax)
    items["ak8_mass"] = Hist("Counts", dataset_ax, region_ax, jet_mass_ax)

    items["dpfcalo"] = Hist("Counts", dataset_ax, region_ax, dpfcalo_ax)
    items["dphijm"] = Hist("min(4 leading jets, MET)", dataset_ax, region_ax, dphi_ax)
    items["dphijr"] = Hist("min(4 leading jets, Recoil)", dataset_ax, region_ax, dphi_ax)

    # Multiplicity histograms
    for cand in ['ak4', 'ak8', 'bjet', 'loose_ele', 'loose_muo', 'tight_ele', 'tight_muo', 'tau', 'photon','gen_dilepton']:
        items[f"{cand}_mult"] = Hist(cand, dataset_ax, region_ax, multiplicity_ax)

    items["muon_pt"] = Hist("Counts", dataset_ax, region_ax, pt_ax)
    items["muon_eta"] = Hist("Counts", dataset_ax, region_ax, eta_ax)
    items["muon_phi"] = Hist("Counts", dataset_ax, region_ax, phi_ax)
    items["muon_eta_phi"] = Hist("Counts", dataset_ax, region_ax, eta_ax, phi_ax)
    items["muon_dxy"] = Hist("Counts", dataset_ax, region_ax, dxy_ax)
    items["muon_dz"] = Hist("Counts", dataset_ax, region_ax, dz_ax)
    items["muon_pt0"] = Hist("Counts", dataset_ax, region_ax, pt_ax)
    items["muon_eta0"] = Hist("Counts", dataset_ax, region_ax, eta_ax)
    items["muon_phi0"] = Hist("Counts", dataset_ax, region_ax, phi_ax)
    items["muon_pt1"] = Hist("Counts", dataset_ax, region_ax, pt_ax)
    items["muon_eta1"] = Hist("Counts", dataset_ax, region_ax, eta_ax)
    items["muon_phi1"] = Hist("Counts", dataset_ax, region_ax, phi_ax)
    items["muon_dxy0"] = Hist("Counts", dataset_ax, region_ax, dxy_ax)
    items["muon_dz0"] = Hist("Counts", dataset_ax, region_ax, dz_ax)
    items["muon_mt"] = Hist("Counts", dataset_ax, region_ax, mt_ax)
    items["dimuon_pt"] = Hist("Counts", dataset_ax, region_ax, pt_ax)
    items["dimuon_eta"] = Hist("Counts", dataset_ax, region_ax, eta_ax)
    items["dimuon_mass"] = Hist("Counts", dataset_ax, region_ax, dilepton_mass_ax)
    items["dimuon_dr"] = Hist("Counts", dataset_ax, region_ax, dr_ax)

    items["electron_pt"] = Hist("Counts", dataset_ax, region_ax, pt_ax)
    items["electron_eta"] = Hist("Counts", dataset_ax, region_ax, eta_ax)
    items["electron_phi"] = Hist("Counts", dataset_ax, region_ax, phi_ax)
    items["electron_eta_phi"] = Hist("Counts", dataset_ax, region_ax, eta_ax, phi_ax)
    items["electron_dxy"] = Hist("Counts", dataset_ax, region_ax, dxy_ax)
    items["electron_dz"] = Hist("Counts", dataset_ax, region_ax, dz_ax)
    items["electron_pt0"] = Hist("Counts", dataset_ax, region_ax, pt_ax)
    items["electron_eta0"] = Hist("Counts", dataset_ax, region_ax, eta_ax)
    items["electron_phi0"] = Hist("Counts", dataset_ax, region_ax, phi_ax)
    items["electron_pt1"] = Hist("Counts", dataset_ax, region_ax, pt_ax)
    items["electron_eta1"] = Hist("Counts", dataset_ax, region_ax, eta_ax)
    items["electron_phi1"] = Hist("Counts", dataset_ax, region_ax, phi_ax)
    items["electron_tightid1"] = Hist("Counts", dataset_ax, region_ax, id_ax)
    items["electron_mt"] = Hist("Counts", dataset_ax, region_ax, mt_ax)
    items["dielectron_pt"] = Hist("Counts", dataset_ax, region_ax, pt_ax)
    items["dielectron_eta"] = Hist("Counts", dataset_ax, region_ax, eta_ax)
    items["dielectron_mass"] = Hist("Counts", dataset_ax, region_ax, dilepton_mass_ax)
    items["dielectron_dr"] = Hist("Counts", dataset_ax, region_ax, dr_ax)

    items['photon_pt0'] = Hist("Counts", dataset_ax, region_ax, pt_ax)
    items['photon_pt0_noweight'] = Hist("Counts", dataset_ax, region_ax, pt_ax)
    items['photon_eta0'] = Hist("Counts", dataset_ax, region_ax, eta_ax)
    items['photon_phi0'] = Hist("Counts", dataset_ax, region_ax, phi_ax)
    items["photon_eta_phi"] = Hist("Counts", dataset_ax, region_ax, eta_ax, phi_ax)
    items['photon_pt0_recoil'] = Hist("Counts", dataset_ax, region_ax, pt_ax, recoil_ax)

    items['drphotonjet'] = Hist("Counts", dataset_ax, region_ax, dr_ax)
    items['drelejet'] = Hist("Counts", dataset_ax, region_ax, dr_ax)
    items['drmuonjet'] = Hist("Counts", dataset_ax, region_ax, dr_ax)

    # One cutflow counter per region
    regions = monojet_regions(cfg).keys()
    for region in regions:
        if region=="inclusive":
            continue
        items[f'cutflow_{region}']  = processor.defaultdict_accumulator(int)

    items['nevents'] = processor.defaultdict_accumulator(float)
    items['sumw'] = processor.defaultdict_accumulator(float)
    items['sumw2'] = processor.defaultdict_accumulator(float)
    items['sumw_pileup'] = processor.defaultdict_accumulator(float)

    items['selected_events'] = processor.defaultdict_accumulator(list)
    items['kinematics'] = processor.defaultdict_accumulator(list)

    for region in ['sr_j','cr_2m_j','cr_1m_j','cr_2e_j','cr_1e_j','cr_g_j']:
        for variable in ['recoil','weight','gen_v_pt']:
            items[f'tree_{region}_{variable}'] = processor.dict_accumulator()

    items['weights'] = Hist("Weights", dataset_ax, region_ax, weight_type_ax, weight_ax)
    items['weights_wide'] = Hist("Weights", dataset_ax, region_ax, weight_type_ax, weight_wide_ax)
    items['npv'] = Hist('Number of primary vertices', dataset_ax, region_ax, nvtx_ax)
    items['npvgood'] = Hist('Number of good primary vertices', dataset_ax, region_ax, nvtx_ax)
    items['npv_nopu'] = Hist('Number of primary vertices (No PU weights)', dataset_ax, region_ax, nvtx_ax)
    items['npvgood_nopu'] = Hist('Number of good primary vertices (No PU weights)', dataset_ax, region_ax, nvtx_ax)

    items['rho_all'] = Hist(r'$\rho$ for all PF candidates', dataset_ax, region_ax, rho_ax)
    items['rho_central'] = Hist(r'$\rho$ for central PF candidates', dataset_ax, region_ax, rho_ax)
    items['rho_all_nopu'] = Hist(r'$\rho$ for all PF candidates (No PU weights)', dataset_ax, region_ax, rho_ax)
    items['rho_central_nopu'] = Hist(r'$\rho$ for central PF candidates (No PU weights)', dataset_ax, region_ax, rho_ax)

    return  processor.dict_accumulator(items)


def setup_candidates(df, cfg):
    if df['is_data'] and extract_year(df['dataset']) != 2018:
        jes_suffix = ''
    else:
        jes_suffix = '_nom'

    muons = JaggedCandidateArray.candidatesfromcounts(
        df['nMuon'],
        pt=df['Muon_pt'],
        eta=df['Muon_eta'],
        abseta=np.abs(df['Muon_eta']),
        phi=df['Muon_phi'],
        mass=0 * df['Muon_pt'],
        charge=df['Muon_charge'],
        looseId=df['Muon_looseId'],
        iso=df["Muon_pfRelIso04_all"],
        tightId=df['Muon_tightId'],
        dxy=df['Muon_dxy'],
        dz=df['Muon_dz']
    )

    # All muons must be at least loose
    muons = muons[muons.looseId \
                    & (muons.iso < cfg.MUON.CUTS.LOOSE.ISO) \
                    & (muons.pt > cfg.MUON.CUTS.LOOSE.PT) \
                    & (muons.abseta<cfg.MUON.CUTS.LOOSE.ETA) \
                    ]


    electrons = JaggedCandidateArray.candidatesfromcounts(
        df['nElectron'],
        pt=df['Electron_pt'],
        eta=df['Electron_eta'],
        abseta=np.abs(df['Electron_eta']),
        phi=df['Electron_phi'],
        mass=0 * df['Electron_pt'],
        charge=df['Electron_charge'],
        looseId=(df[cfg.ELECTRON.BRANCH.ID]>=1),
        tightId=(df[cfg.ELECTRON.BRANCH.ID]==4),
        dxy=np.abs(df['Electron_dxy']),
        dz=np.abs(df['Electron_dz']),
        barrel=np.abs(df['Electron_eta']) <= 1.479
    )
    # All electrons must be at least loose
    pass_dxy = (electrons.barrel & (np.abs(electrons.dxy) < cfg.ELECTRON.CUTS.LOOSE.DXY.BARREL)) \
    | (~electrons.barrel & (np.abs(electrons.dxy) < cfg.ELECTRON.CUTS.LOOSE.DXY.ENDCAP))

    pass_dz = (electrons.barrel & (np.abs(electrons.dz) < cfg.ELECTRON.CUTS.LOOSE.DZ.BARREL)) \
    | (~electrons.barrel & (np.abs(electrons.dz) < cfg.ELECTRON.CUTS.LOOSE.DZ.ENDCAP))

    electrons = electrons[electrons.looseId \
                                    & (electrons.pt>cfg.ELECTRON.CUTS.LOOSE.PT) \
                                    & (electrons.abseta<cfg.ELECTRON.CUTS.LOOSE.ETA) \
                                    & pass_dxy \
                                    & pass_dz
                                    ]


    taus = JaggedCandidateArray.candidatesfromcounts(
        df['nTau'],
        pt=df['Tau_pt'],
        eta=df['Tau_eta'],
        abseta=np.abs(df['Tau_eta']),
        phi=df['Tau_phi'],
        mass=0 * df['Tau_pt'],
        decaymode=df['Tau_idDecayMode'],
        iso=df['Tau_idMVAoldDM2017v2'],
    )

    taus = taus[ object_overlap(taus, muons) \
                & object_overlap(taus, electrons) \
                & (taus.decaymode) \
                & (taus.pt > cfg.TAU.CUTS.PT)\
                & (taus.abseta < cfg.TAU.CUTS.ETA) \
                & ((taus.iso&2)==2)]

    photons = JaggedCandidateArray.candidatesfromcounts(
        df['nPhoton'],
        pt=df['Photon_pt'],
        eta=df['Photon_eta'],
        abseta=np.abs(df['Photon_eta']),
        phi=df['Photon_phi'],
        mass=0*df['Photon_pt'],
        looseId=(df[cfg.PHOTON.BRANCH.ID]>=1) & df['Photon_electronVeto'],
        mediumId=(df[cfg.PHOTON.BRANCH.ID]>=2) & df['Photon_electronVeto'],
        r9=df['Photon_r9'],
        barrel=np.abs(df['Photon_eta']) < 1.479,
    )
    photons = photons[photons.looseId \
              & (photons.pt > cfg.PHOTON.CUTS.LOOSE.pt) \
              & (photons.abseta < cfg.PHOTON.CUTS.LOOSE.eta) \
              & object_overlap(photons, muons, dr=0.5) \
              & object_overlap(photons, electrons, dr=0.5)
              ]

    ak4 = JaggedCandidateArray.candidatesfromcounts(
        df['nJet'],
        pt=df[f'Jet_pt{jes_suffix}'],
        eta=df['Jet_eta'],
        abseta=np.abs(df['Jet_eta']),
        phi=df['Jet_phi'],
        mass=np.zeros_like(df['Jet_pt']),
        looseId=(df['Jet_jetId']&4) == 4, # bitmask: 1 = loose, 2 = tight, 3 = tight + lep veto
        tightId=(df['Jet_jetId']&4) == 4, # bitmask: 1 = loose, 2 = tight, 3 = tight + lep veto
        puid=((df['Jet_puId']&2>0) | (df[f'Jet_pt{jes_suffix}']>50)), # medium pileup jet ID
        csvv2=df["Jet_btagCSVV2"],
        deepcsv=df['Jet_btagDeepB'],
        # nef=df['Jet_neEmEF'],
        nhf=df['Jet_neHEF'],
        chf=df['Jet_chHEF'],
        ptraw=df['Jet_pt']*(1-df['Jet_rawFactor']),
        nconst=df['Jet_nConstituents']
        # clean=df['Jet_cleanmask']
        # cef=df['Jet_chEmEF'],
    )
    # Before cleaning, apply HEM veto
    hem_ak4 = ak4[ (ak4.pt>30) &
        (-3.0 < ak4.eta) &
        (ak4.eta < -1.3) &
        (-1.57 < ak4.phi) &
        (ak4.phi < -0.87)
        ]
    df['hemveto'] = hem_ak4.counts == 0

    ak4 = ak4[ak4.looseId & object_overlap(ak4, muons) & object_overlap(ak4, electrons) & object_overlap(ak4, photons)]

    ak8 = JaggedCandidateArray.candidatesfromcounts(
        df['nFatJet'],
        pt=df[f'FatJet_pt{jes_suffix}'],
        eta=df['FatJet_eta'],
        abseta=np.abs(df['FatJet_eta']),
        phi=df['FatJet_phi'],
        mass=df[f'FatJet_msoftdrop{jes_suffix}'],
        tightId=(df['FatJet_jetId']&2) == 2, # Tight
        csvv2=df["FatJet_btagCSVV2"],
        deepcsv=df['FatJet_btagDeepB'],
        tau1=df['FatJet_tau1'],
        tau2=df['FatJet_tau2'],
        tau21=df['FatJet_tau2']/df['FatJet_tau1'],
        wvsqcd=df['FatJet_deepTag_WvsQCD'],
        wvsqcdmd=df['FatJet_deepTagMD_WvsQCD'],
        zvsqcd=df['FatJet_deepTag_ZvsQCD'],
        zvsqcdmd=df['FatJet_deepTagMD_ZvsQCD']
    )
    ak8 = ak8[ak8.tightId & object_overlap(ak8, muons) & object_overlap(ak8, electrons) & object_overlap(ak8, photons)]

    met_pt = df[f'MET_pt{jes_suffix}']
    met_phi = df[f'MET_phi{jes_suffix}']

    return met_pt, met_phi, ak4, ak8, muons, electrons, taus, photons


def monojet_regions(cfg):
    common_cuts = [
        'filt_met',
        'hemveto',
        'veto_ele',
        'veto_muo',
        'veto_photon',
        'veto_tau',
        'veto_b',
        'mindphijr',
        'dpfcalo',
        'recoil'
    ]
    j_cuts = [
        'leadak4_pt_eta',
        'leadak4_id',
        # 'veto_vtag'
    ]
    v_cuts = [
        'leadak8_pt_eta',
        'leadak8_id',
        'leadak8_mass',
        'leadak8_tau21',
    ]

    regions = {}
    regions['inclusive'] = ['inclusive']

    # Signal regions (v = mono-V, j = mono-jet)
    regions['sr_v'] = ['trig_met'] + common_cuts + v_cuts
    regions['sr_j'] = ['trig_met'] + common_cuts + j_cuts

    regions['cr_nofilt_j'] = copy.deepcopy(regions['sr_j'])
    regions['cr_nofilt_j'].remove('filt_met')

    # Dimuon CR
    cr_2m_cuts = ['trig_met','two_muons', 'at_least_one_tight_mu', 'dimuon_mass', 'dimuon_charge'] + common_cuts
    cr_2m_cuts.remove('veto_muo')

    regions['cr_2m_j'] = cr_2m_cuts + j_cuts
    regions['cr_2m_v'] = cr_2m_cuts + v_cuts

    # Single muon CR
    cr_1m_cuts = ['trig_met','one_muon', 'at_least_one_tight_mu', 'mt_mu'] + common_cuts
    cr_1m_cuts.remove('veto_muo')
    regions['cr_1m_j'] = cr_1m_cuts + j_cuts
    regions['cr_1m_v'] = cr_1m_cuts + v_cuts

    # Dielectron CR
    cr_2e_cuts = ['trig_ele','two_electrons', 'at_least_one_tight_el', 'dielectron_mass', 'dielectron_charge'] + common_cuts
    cr_2e_cuts.remove('veto_ele')
    regions['cr_2e_j'] = cr_2e_cuts + j_cuts
    regions['cr_2e_v'] = cr_2e_cuts + v_cuts

    # Single electron CR
    cr_1e_cuts = ['trig_ele','one_electron', 'at_least_one_tight_el', 'met_el','mt_el'] + common_cuts
    cr_1e_cuts.remove('veto_ele')
    regions['cr_1e_j'] =  cr_1e_cuts + j_cuts
    regions['cr_1e_v'] =  cr_1e_cuts + v_cuts

    # Photon CR
    cr_g_cuts = ['trig_photon', 'one_photon', 'at_least_one_tight_photon','photon_pt'] + common_cuts
    cr_g_cuts.remove('veto_photon')

    regions['cr_g_j'] = cr_g_cuts + j_cuts
    regions['cr_g_v'] = cr_g_cuts + v_cuts

    if not cfg.RUN.MONOV:
        keys_to_remove = [ x for x in regions.keys() if x.endswith('_v')]
        for key in keys_to_remove:
            regions.pop(key)

    if cfg.RUN.TRIGGER_STUDY:
        # Trigger studies
        # num = numerator, den = denominator
        # Single Mu region: Remove recoil cut, add SingleMu trigger, toggle MET trigger
        tr_1m_num_cuts = copy.deepcopy(cr_1m_cuts) + j_cuts
        tr_1m_num_cuts.remove('recoil')
        tr_1m_num_cuts.append('trig_mu')
        tr_1m_num_cuts.append('mu_pt_trig_safe')
        regions['tr_1m_num'] = tr_1m_num_cuts

        tr_1m_den_cuts = copy.deepcopy(tr_1m_num_cuts)
        tr_1m_den_cuts.remove('trig_met')
        regions['tr_1m_den'] = tr_1m_den_cuts

        # Add HLT muon requirement
        regions['tr_1m_hlt_num'] = regions['tr_1m_num'] + ['one_hlt_muon']
        regions['tr_1m_hlt_den'] = regions['tr_1m_den'] + ['one_hlt_muon']

        # Double Mu region: Remove recoil cut, toggle MET trigger
        tr_2m_num_cuts = copy.deepcopy(cr_2m_cuts) + j_cuts
        tr_2m_num_cuts.remove('recoil')
        tr_2m_num_cuts.append('trig_mu')
        tr_2m_num_cuts.append('mu_pt_trig_safe')
        regions['tr_2m_num'] = tr_2m_num_cuts

        tr_2m_den_cuts = copy.deepcopy(tr_2m_num_cuts)
        tr_2m_den_cuts.remove('trig_met')
        regions['tr_2m_den'] = tr_2m_den_cuts

        regions['tr_2m_hlt_den'] = regions['tr_2m_den'] + ['two_hlt_muons']
        regions['tr_2m_hlt_num'] = regions['tr_2m_num'] + ['two_hlt_muons']

        # Single Electron region: Remove recoil cut, toggle MET trigger
        tr_1e_num_cuts = copy.deepcopy(cr_1e_cuts) + j_cuts
        tr_1e_num_cuts.remove('recoil')
        tr_1e_num_cuts.append('trig_met')
        regions['tr_1e_num'] = tr_1e_num_cuts

        tr_1e_den_cuts = copy.deepcopy(tr_1e_num_cuts)
        tr_1e_den_cuts.remove('trig_met')
        regions['tr_1e_den'] = tr_1e_den_cuts

        # Double Electron region: Remove recoil cut, toggle MET trigger
        tr_2e_num_cuts = copy.deepcopy(cr_2e_cuts) + j_cuts
        tr_2e_num_cuts.remove('recoil')
        tr_2e_num_cuts.append('trig_met')
        regions['tr_2e_num'] = tr_2e_num_cuts

        tr_2e_den_cuts = copy.deepcopy(tr_2e_num_cuts)
        tr_2e_den_cuts.remove('trig_met')
        regions['tr_2e_den'] = tr_2e_den_cuts

        # Photon region
        tr_g_num_cuts = copy.deepcopy(cr_g_cuts) + j_cuts
        tr_g_num_cuts.remove('recoil')
        tr_g_num_cuts.remove('photon_pt')

        tr_g_den_cuts = copy.deepcopy(tr_g_num_cuts)
        tr_g_den_cuts.remove('trig_photon')

        regions[f'tr_g_notrig_num'] = copy.deepcopy(tr_g_num_cuts)
        regions[f'tr_g_notrig_den'] = copy.deepcopy(tr_g_den_cuts)

        for trgname in cfg.TRIGGERS.HT.GAMMAEFF:
            regions[f'tr_g_{trgname}_num'] = tr_g_num_cuts + [trgname]
            regions[f'tr_g_{trgname}_den'] = tr_g_den_cuts + [trgname]

            regions[f'tr_g_{trgname}_photon_pt_trig_cut_num'] = tr_g_num_cuts + [trgname, 'photon_pt_trig']
            regions[f'tr_g_{trgname}_photon_pt_trig_cut_den'] = tr_g_den_cuts + [trgname, 'photon_pt_trig']


    # HEM check
    for region in ['cr_2m_j','cr_1m_j','cr_2e_j','cr_1e_j']:
        tmp = copy.deepcopy(regions[region])
        tmp.remove('hemveto')
        regions[f'{region}_nohem'] = tmp

    for region in ['cr_2m_j','cr_1m_j','cr_2e_j','cr_1e_j','cr_g_j']:
        for cut in ['veto_ele', 'veto_muo', 'veto_photon', 'veto_tau', 'veto_b']:
            if cut in regions[region]:
                regions[f'{region}_no{cut}'] = copy.deepcopy(regions[region])
                regions[f'{region}_no{cut}'].remove(cut)
        regions[f'{region}_loose'] = copy.deepcopy(regions[region])
        regions[f'{region}_loose'].remove('mindphijr')
        regions[f'{region}_loose'].remove('dpfcalo')

    return regions


def theory_weights_monojet(weights, df, evaluator, gen_v_pt):
    if df['is_lo_w']:
        weights.add("theory", evaluator["qcd_nlo_w_2017"](gen_v_pt) * evaluator["qcd_nnlo_w"](gen_v_pt) * evaluator["ewk_nlo_w"](gen_v_pt))
    elif df['is_lo_z']:
        weights.add("theory", evaluator["qcd_nlo_z_2017"](gen_v_pt) * evaluator["qcd_nnlo_z"](gen_v_pt) * evaluator["ewk_nlo_z"](gen_v_pt))
    elif df['is_nlo_w']:
        weights.add("theory", evaluator["qcd_nnlo_w"](gen_v_pt) * evaluator["ewk_nlo_w"](gen_v_pt))
    elif df['is_nlo_z']:
        weights.add("theory", evaluator["qcd_nnlo_z"](gen_v_pt) * evaluator["ewk_nlo_z"](gen_v_pt))
    elif df['is_lo_g']:
        weights.add("theory", evaluator["ewk_nlo_g"](gen_v_pt) * evaluator["qcd_nlo_g"](gen_v_pt) * evaluator["qcd_nnlo_g"](gen_v_pt))
    else:
        weights.add("theory", np.ones(df.size))

    return weights

def theory_weights_vbf(weights, df, evaluator, gen_v_pt, mjj):
    if df['is_lo_w']:
        weights.add("theory", evaluator["qcd_nlo_w_2017"](gen_v_pt) * evaluator["qcd_nnlo_w"](gen_v_pt) * evaluator["ewk_nlo_w"](gen_v_pt))
    elif df['is_lo_w_ewk']:
        weights.add("theory", evaluator["qcd_nlo_w_ewk"](gen_v_pt, mjj))
    elif df['is_lo_z']:
        weights.add("theory", evaluator["qcd_nlo_z_2017"](gen_v_pt) * evaluator["qcd_nnlo_z"](gen_v_pt) * evaluator["ewk_nlo_z"](gen_v_pt))
    elif df['is_lo_z_ewk']:
        weights.add("theory", evaluator["qcd_nlo_z_ewk"](gen_v_pt, mjj))
    elif df['is_nlo_w']:
        weights.add("theory", evaluator["qcd_nnlo_w"](gen_v_pt) * evaluator["ewk_nlo_w"](gen_v_pt))
    elif df['is_nlo_z']:
        weights.add("theory", evaluator["qcd_nnlo_z"](gen_v_pt) * evaluator["ewk_nlo_z"](gen_v_pt))
    elif df['is_lo_g']:
        weights.add("theory", evaluator["ewk_nlo_g"](gen_v_pt) * evaluator["qcd_nlo_g"](gen_v_pt) * evaluator["qcd_nnlo_g"](gen_v_pt))
    else:
        weights.add("theory", np.ones(df.size))

    return weights

def pileup_weights(weights, df, evaluator, cfg):
    if cfg.SF.PILEUP.MODE == 'nano':
                weights.add("pileup", df['puWeight'])
    elif cfg.SF.PILEUP.MODE == 'manual':
        weights.add("pileup", evaluator['pileup'](df['Pileup_nTrueInt']))
    else:
        raise RuntimeError(f"Unknown value for cfg.PILEUP.MODE: {cfg.PILEUP.MODE}.")
    return weights



def candidate_weights(weights, df, evaluator, muons, electrons, photons):
    # Muon ID and Isolation for tight and loose WP
    # Function of pT, eta (Order!)
    weights.add("muon_id_tight", evaluator['muon_id_tight'](muons[df['is_tight_muon']].pt, muons[df['is_tight_muon']].abseta).prod())
    weights.add("muon_iso_tight", evaluator['muon_iso_tight'](muons[df['is_tight_muon']].pt, muons[df['is_tight_muon']].abseta).prod())
    weights.add("muon_id_loose", evaluator['muon_id_loose'](muons[~df['is_tight_muon']].pt, muons[~df['is_tight_muon']].abseta).prod())
    weights.add("muon_iso_loose", evaluator['muon_iso_loose'](muons[~df['is_tight_muon']].pt, muons[~df['is_tight_muon']].abseta).prod())

    # Electron ID and reco
    # Function of eta, pT (Other way round relative to muons!)
    weights.add("ele_reco", evaluator['ele_reco'](electrons.eta, electrons.pt).prod())
    weights.add("ele_id_tight", evaluator['ele_id_tight'](electrons[df['is_tight_electron']].eta, electrons[df['is_tight_electron']].pt).prod())
    weights.add("ele_id_loose", evaluator['ele_id_loose'](electrons[~df['is_tight_electron']].eta, electrons[~df['is_tight_electron']].pt).prod())

    # Photon ID and electron veto
    weights.add("photon_id_tight", evaluator['photon_id_tight'](photons[df['is_tight_photon']].eta, photons[df['is_tight_photon']].pt).prod())

    # CSEV not split only by EE/EB for now
    csev_sf_index = 0.5 * photons.barrel + 3.5 * ~photons.barrel + 1 * (photons.r9 > 0.94) + 2 * (photons.r9 <= 0.94)
    weights.add("photon_csev", evaluator['photon_csev'](csev_sf_index).prod())

    return weights
