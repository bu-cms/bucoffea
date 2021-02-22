import re
import copy
from bucoffea.helpers.gen import setup_gen_jets_ak8

import coffea.processor as processor
import numpy as np
from coffea import hist
from coffea.analysis_objects import JaggedCandidateArray

from bucoffea.helpers import object_overlap, sigmoid, exponential
from bucoffea.helpers.dataset import extract_year

Hist = hist.Hist
Bin = hist.Bin
Cat = hist.Cat

def empty_column_accumulator_int():
    return processor.column_accumulator(np.array([],dtype=np.uint64))

def empty_column_accumulator_int64():
    return processor.column_accumulator(np.array([],dtype=np.int64))
def empty_column_accumulator_float64():
    return processor.column_accumulator(np.array([],dtype=np.float64))
def empty_column_accumulator_float32():
    return processor.column_accumulator(np.array([],dtype=np.float32))
def empty_column_accumulator_float16():
    return processor.column_accumulator(np.array([],dtype=np.float16))
def empty_column_accumulator_bool():
    return processor.column_accumulator(np.array([],dtype=np.bool))


def accu_int():
    return processor.defaultdict_accumulator(int)

def defaultdict_accumulator_of_empty_column_accumulator_int64():
    return processor.defaultdict_accumulator(empty_column_accumulator_int64)
def defaultdict_accumulator_of_empty_column_accumulator_float64():
    return processor.defaultdict_accumulator(empty_column_accumulator_float64)
def defaultdict_accumulator_of_empty_column_accumulator_float32():
    return processor.defaultdict_accumulator(empty_column_accumulator_float32)
def defaultdict_accumulator_of_empty_column_accumulator_float16():
    return processor.defaultdict_accumulator(empty_column_accumulator_float16)
def defaultdict_accumulator_of_empty_column_accumulator_bool():
    return processor.defaultdict_accumulator(empty_column_accumulator_bool)

def monojet_accumulator(cfg):
    dataset_ax = Cat("dataset", "Primary dataset")
    region_ax = Cat("region", "Selection region")
    type_ax = Cat("type", "Type")
    variation_ax = Cat("variation","Variation")
    wppass_ax = Bin("wppass", "WP Pass",2,-0.5,1.5)
    vpt_ax = Bin("vpt",r"$p_{T}^{V}$ (GeV)", 50, 0, 2000)

    met_ax = Bin("met", r"$p_{T}^{miss}$ (GeV)", 40, 0, 2000)
    recoil_ax = Bin("recoil", r"Recoil (GeV)", 200, 0, 2000)
    recoil_ax_coarse = Bin("recoil", r"Recoil (GeV)", 20, 0, 2000)
    gen_v_pt_ax = Bin("pt", r"pt (GeV)", 200, 0, 2000)

    jet_pt_ax = Bin("jetpt", r"$p_{T}$ (GeV)", 50, 0, 1000)
    jet_eta_ax = Bin("jeteta", r"$\eta$", 50, -5, 5)
    jet_eta_ax_coarse = Bin("jeteta", r"$\eta$", 10, -5, 5)
    jet_phi_ax = Bin("jetphi", r"$\phi$", 50,-np.pi, np.pi)

    jet_mass_ax = Bin("mass", r"$M_{jet}$ (GeV)", 100,0,300)

    dpfcalo_ax = Bin("dpfcalo", r"$(CaloMET-PFMET) / Recoil$", 20, -1, 1)
    dpftk_ax = Bin("dpfcalo", r"$(TKMet-PFMET) / Recoil$", 20, -1, 1)
    btag_ax = Bin("btag", r"B tag discriminator", 20, 0, 1)
    multiplicity_ax = Bin("multiplicity", r"multiplicity", 10, -0.5, 9.5)
    dphi_ax = Bin("dphi", r"$\Delta\phi$", 50, 0, 3.5)
    dr_ax = Bin("dr", r"$\Delta R$", 50, 0, 2)
    vec_b_ax = Bin("vec_b", "vec_b", 50, 0,1)
    vec_dphi_ax = Bin("vec_dphi", "vec_dphi", 50, 0, np.pi)

    dxy_ax = Bin("dxy", r"$d_{xy}$", 20, 0, 0.5)
    dz_ax = Bin("dz", r"$d_{z}$", 20, 0, 0.5)

    tmp =[0, 0.25, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.92,0.94,0.95,0.96,0.97,0.98,0.99, 1.0]
    tmp = tmp + [2-x for x in reversed(tmp[:-1])]
    response_ax = Bin("response", r"response", tmp)
    id_ax = Bin("id", r"ID bool", 2,-0.5,1.5)

    pt_ax = Bin("pt", r"$p_{T}$ (GeV)", 50, 0, 1000)
    ht_ax = Bin("ht", r"$H_{T}$ (GeV)", 50, 0, 4000)
    mt_ax = Bin("mt", r"$M_{T}$ (GeV)", 50, 0, 1000)
    eta_ax = Bin("eta", r"$\eta$", 50, -5, 5)
    eta_ax_coarse = Bin("eta", r"$\eta$", 25, -5, 5)
    phi_ax = Bin("phi", r"$\phi$", 50,-np.pi, np.pi)
    phi_ax_coarse = Bin("phi", r"$\phi$", 20,-np.pi, np.pi)

    ratio_ax = Bin("ratio", "ratio", 50,0,2)

    tau21_ax = Bin("tau21", r"Tagger", 50,0,1)
    tagger_ax = Bin("tagger", r"Tagger", 50,0,1)

    dilepton_mass_ax = Bin("dilepton_mass", r"$M(\ell\ell)$ (GeV)", 50,50,150)

    weight_type_ax = Cat("weight_type", "Weight type")
    weight_ax = Bin("weight_value", "Weight",20,0.5,1.5)
    sf_ax = Bin("sf", "sf",50,0.8,1.8)

    nvtx_ax = Bin('nvtx','Number of vertices',50,-0.5,99.5)
    rho_ax = Bin('rho','Energy density',50, 0, 100)
    frac_ax = Bin('frac','Fraction', 50, 0, 1)
    deepcsv_ax = Bin('deepcsv','DeepCSV discriminator', 50, 0, 1)
    Hist = hist.Hist
    items = {}
    items["genvpt_check"] = Hist("Counts", dataset_ax, type_ax, vpt_ax)
    items["lhe_ht"] = Hist("Counts", dataset_ax, ht_ax)
    items["met"] = Hist("Counts", dataset_ax, region_ax, met_ax)
    items["met_phi"] = Hist("Counts", dataset_ax, region_ax, phi_ax)
    items["recoil"] = Hist("Counts", dataset_ax, region_ax, recoil_ax)
    items["recoil_veto_weight"] = Hist("Counts", dataset_ax, region_ax, recoil_ax,variation_ax)
    items["recoil_nopog"] = Hist("Counts", dataset_ax, region_ax, recoil_ax)
    items["recoil_nopu"] = Hist("Counts", dataset_ax, region_ax, recoil_ax)
    items["recoil_notrg"] = Hist("Counts", dataset_ax, region_ax, recoil_ax)
    items["recoil_nopref"] = Hist("Counts", dataset_ax, region_ax, recoil_ax)

    items["recoil_nodibosonnlo"] = Hist("Counts", dataset_ax, region_ax, recoil_ax)
    items["recoil_dibosonnlo_up"] = Hist("Counts", dataset_ax, region_ax, recoil_ax)
    items["recoil_dibosonnlo_dn"] = Hist("Counts", dataset_ax, region_ax, recoil_ax)

    items["gen_v_pt"] = Hist("Counts", dataset_ax, region_ax, gen_v_pt_ax)
    items["gen_v_pt_notheory"] = Hist("Counts", dataset_ax, region_ax, gen_v_pt_ax)

    if cfg and cfg.RUN.BTAG_STUDY:
        items["recoil_hardbveto"] = Hist("Counts", dataset_ax, region_ax, recoil_ax)
        items["recoil_bveto_up"] = Hist("Counts", dataset_ax, region_ax, recoil_ax)
        items["recoil_bveto_down"] = Hist("Counts", dataset_ax, region_ax, recoil_ax)
    if cfg and cfg.RUN.PHOTON_ID_STUDY:
        items["recoil_photon_id_up"] = Hist("Counts", dataset_ax, region_ax, recoil_ax)
        items["recoil_photon_id_dn"] = Hist("Counts", dataset_ax, region_ax, recoil_ax)
        items["recoil_photon_id_extrap_up"] = Hist("Counts", dataset_ax, region_ax, recoil_ax)
        items["recoil_photon_id_extrap_dn"] = Hist("Counts", dataset_ax, region_ax, recoil_ax)
    if cfg and cfg.RUN.ELE_ID_STUDY:
        items["recoil_ele_id_up"] = Hist("Counts", dataset_ax, region_ax, recoil_ax)
        items["recoil_ele_id_dn"] = Hist("Counts", dataset_ax, region_ax, recoil_ax)
        items["recoil_ele_id_nm"] = Hist("Counts", dataset_ax, region_ax, recoil_ax)
        items["recoil_ele_reco_up"] = Hist("Counts", dataset_ax, region_ax, recoil_ax)
        items["recoil_ele_reco_dn"] = Hist("Counts", dataset_ax, region_ax, recoil_ax)
    items["recoil_phi"] = Hist("Counts", dataset_ax, region_ax, phi_ax)
    items["ak4_pt0_over_recoil"] = Hist("Counts", dataset_ax, region_ax, ratio_ax)
    items["ak4_pt0"] = Hist("Counts", dataset_ax, region_ax, jet_pt_ax)
    items["ak4_ptraw0"] = Hist("Counts", dataset_ax, region_ax, jet_pt_ax)
    items["ak4_eta0"] = Hist("Counts", dataset_ax, region_ax, jet_eta_ax)
    items["ak4_phi0"] = Hist("Counts", dataset_ax, region_ax, jet_phi_ax)
    items["ak4_chf0"] = Hist("Counts", dataset_ax, region_ax, frac_ax)
    items["ak4_nhf0"] = Hist("Counts", dataset_ax, region_ax, frac_ax)
    items["ak4_nef0"] = Hist("Counts", dataset_ax, region_ax, frac_ax)
    items["ak4_muf0"] = Hist("Counts", dataset_ax, region_ax, frac_ax)
    items["ak4_cef0"] = Hist("Counts", dataset_ax, region_ax, frac_ax)
    items["ak4_chf"] = Hist("Counts", dataset_ax, region_ax, frac_ax)
    items["ak4_nhf"] = Hist("Counts", dataset_ax, region_ax, frac_ax)
    items["ak4_nef"] = Hist("Counts", dataset_ax, region_ax, frac_ax)
    items["ak4_muf"] = Hist("Counts", dataset_ax, region_ax, frac_ax)
    items["ak4_cef"] = Hist("Counts", dataset_ax, region_ax, frac_ax)

    items["ak4_pt0_eta0"] = Hist("Counts", dataset_ax, region_ax, jet_pt_ax,jet_eta_ax_coarse)
    items["ak4_pt0_recoil"] = Hist("Counts", dataset_ax, region_ax, jet_pt_ax,recoil_ax_coarse)

    items["ak4_pt"] = Hist("Counts", dataset_ax, region_ax, jet_pt_ax)
    items["ak4_eta"] = Hist("Counts", dataset_ax, region_ax, jet_eta_ax)
    items["ak4_phi"] = Hist("Counts", dataset_ax, region_ax, jet_phi_ax)
    items["ak4_deepcsv"] = Hist("Counts", dataset_ax, region_ax, deepcsv_ax)
    items["bjet_pt"] = Hist("Counts", dataset_ax, region_ax, jet_pt_ax)
    items["bjet_eta"] = Hist("Counts", dataset_ax, region_ax, jet_eta_ax)
    items["bjet_phi"] = Hist("Counts", dataset_ax, region_ax, jet_phi_ax)
    items["bjet_sf"] = Hist("Counts", dataset_ax, region_ax, sf_ax)
    items["ak4_eta_phi"] = Hist("Counts", dataset_ax, region_ax, eta_ax_coarse, phi_ax_coarse)
    items["ak4_eta0_phi0"] = Hist("Counts", dataset_ax, region_ax, eta_ax_coarse, phi_ax_coarse)
    items["ak4_btag"] = Hist("Counts", dataset_ax, region_ax, btag_ax)

    items["ak8_pt0"] = Hist("Counts", dataset_ax, region_ax, jet_pt_ax)
    items["ak8_pt0_recoil"] = Hist("Counts", dataset_ax, region_ax, jet_pt_ax, recoil_ax_coarse)
    items["ak8_eta0"] = Hist("Counts", dataset_ax, region_ax, jet_eta_ax)
    items["ak8_phi0"] = Hist("Counts", dataset_ax, region_ax, jet_phi_ax)
    items["ak8_mass0"] = Hist("Counts", dataset_ax, region_ax, jet_mass_ax)
    items["ak8_mass_response"] = Hist("Counts", dataset_ax, region_ax, response_ax)
    items["ak8_tau210"] = Hist("Counts", dataset_ax, region_ax, tau21_ax)
    items["ak8_wvsqcd0"] = Hist("Counts", dataset_ax, region_ax, tagger_ax)
    items["ak8_wvsqcdmd0"] = Hist("Counts", dataset_ax, region_ax, tagger_ax)
    items["ak8_zvsqcd0"] = Hist("Counts", dataset_ax, region_ax, tagger_ax)
    items["ak8_zvsqcdmd0"] = Hist("Counts", dataset_ax, region_ax, tagger_ax)
    items["ak8_tvsqcd0"] = Hist("Counts", dataset_ax, region_ax, tagger_ax)
    items["ak8_tvsqcdmd0"] = Hist("Counts", dataset_ax, region_ax, tagger_ax)
    items["ak8_wvstqcd0"] = Hist("Counts", dataset_ax, region_ax, tagger_ax)
    items["ak8_wvstqcdmd0"] = Hist("Counts", dataset_ax, region_ax, tagger_ax)
    items["ak8_eta_phi"] = Hist("Counts", dataset_ax, region_ax, eta_ax_coarse, phi_ax_coarse)
    items["ak8_eta0_phi0"] = Hist("Counts", dataset_ax, region_ax, eta_ax_coarse, phi_ax_coarse)

    items["trailak8_ak4_pt"] = Hist("Counts", dataset_ax, region_ax, jet_pt_ax)
    items["trailak8_ak4_dr_min"] = Hist("Counts", dataset_ax, region_ax, dr_ax)

    if cfg and cfg.RUN.MONOVMISTAG_STUDY:
        items["ak8_passloose_pt0"] = Hist("Counts", dataset_ax, region_ax, wppass_ax, jet_pt_ax)
        items["ak8_passmedium_pt0"] = Hist("Counts", dataset_ax, region_ax, wppass_ax, jet_pt_ax)
        items["ak8_passtight_pt0"] = Hist("Counts", dataset_ax, region_ax, wppass_ax, jet_pt_ax)
        items["ak8_passloosemd_pt0"] = Hist("Counts", dataset_ax, region_ax, wppass_ax, jet_pt_ax)
        items["ak8_passtightmd_pt0"] = Hist("Counts", dataset_ax, region_ax, wppass_ax, jet_pt_ax)
        items["ak8_passloose_mass0"] = Hist("Counts", dataset_ax, region_ax, wppass_ax, jet_mass_ax)
        items["ak8_passmedium_mass0"] = Hist("Counts", dataset_ax, region_ax, wppass_ax, jet_mass_ax)
        items["ak8_passtight_mass0"] = Hist("Counts", dataset_ax, region_ax, wppass_ax, jet_mass_ax)
        items["ak8_passloosemd_mass0"] = Hist("Counts", dataset_ax, region_ax, wppass_ax, jet_mass_ax)
        items["ak8_passtightmd_mass0"] = Hist("Counts", dataset_ax, region_ax, wppass_ax, jet_mass_ax)
    items["ak8_Vmatched_pt0"] = Hist("Counts", dataset_ax, region_ax, jet_pt_ax)

    items["ak8_pt"] = Hist("Counts", dataset_ax, region_ax, jet_pt_ax)
    items["ak8_eta"] = Hist("Counts", dataset_ax, region_ax, jet_eta_ax)
    items["ak8_phi"] = Hist("Counts", dataset_ax, region_ax, jet_phi_ax)
    items["ak8_mass"] = Hist("Counts", dataset_ax, region_ax, jet_mass_ax)
    items["lowmass_ak8_pt"] = Hist("Counts", dataset_ax, region_ax, jet_pt_ax)
    items["lowmass_ak8_eta"] = Hist("Counts", dataset_ax, region_ax, jet_eta_ax)
    items["lowmass_ak8_phi"] = Hist("Counts", dataset_ax, region_ax, jet_phi_ax)
    items["lowmass_ak8_mass"] = Hist("Counts", dataset_ax, region_ax, jet_mass_ax)
    items["vlowmass_ak8_pt"] = Hist("Counts", dataset_ax, region_ax, jet_pt_ax)
    items["vlowmass_ak8_eta"] = Hist("Counts", dataset_ax, region_ax, jet_eta_ax)
    items["vlowmass_ak8_phi"] = Hist("Counts", dataset_ax, region_ax, jet_phi_ax)
    items["vlowmass_ak8_mass"] = Hist("Counts", dataset_ax, region_ax, jet_mass_ax)

    items["dpfcalo"] = Hist("Counts", dataset_ax, region_ax, dpfcalo_ax)
    items["dpftk"] = Hist("Counts", dataset_ax, region_ax, dpftk_ax)
    items["dphijm"] = Hist("min(4 leading jets, MET)", dataset_ax, region_ax, dphi_ax)
    items["dphijr"] = Hist("min(4 leading jets, Recoil)", dataset_ax, region_ax, dphi_ax)

    items["vec_dphi"] = Hist("vec_dphi", dataset_ax, region_ax, vec_dphi_ax)
    items["vec_b"] = Hist("vec_b", dataset_ax, region_ax, vec_b_ax)

    items["dphi_pf_tk"] = Hist(r"\Delta\Phi(PF, TK)", dataset_ax, region_ax, dphi_ax)
    items["dphi_pf_calo"] = Hist(r"\Delta\Phi(PF, Calo)", dataset_ax, region_ax, dphi_ax)

    # Multiplicity histograms
    for cand in ['ak4', 'ak8', 'bjet', 'loose_ele', 'loose_muo', 'tight_ele', 'tight_muo', 'tau', 'photon','gen_dilepton']:
        items[f"{cand}_mult"] = Hist(cand, dataset_ax, region_ax, multiplicity_ax)

    items["muon_pt"] = Hist("Counts", dataset_ax, region_ax, pt_ax)
    items["muon_eta"] = Hist("Counts", dataset_ax, region_ax, eta_ax)
    items["muon_phi"] = Hist("Counts", dataset_ax, region_ax, phi_ax)
    items["muon_eta_phi"] = Hist("Counts", dataset_ax, region_ax, eta_ax_coarse, phi_ax_coarse)
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
    items["electron_eta_phi"] = Hist("Counts", dataset_ax, region_ax, eta_ax_coarse, phi_ax_coarse)
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
    items['photon_eta0'] = Hist("Counts", dataset_ax, region_ax, eta_ax)
    items['photon_phi0'] = Hist("Counts", dataset_ax, region_ax, phi_ax)
    items["photon_eta_phi"] = Hist("Counts", dataset_ax, region_ax, eta_ax_coarse, phi_ax_coarse)

    items['drphotonjet'] = Hist("Counts", dataset_ax, region_ax, dr_ax)
    items['drelejet'] = Hist("Counts", dataset_ax, region_ax, dr_ax)
    items['drmuonjet'] = Hist("Counts", dataset_ax, region_ax, dr_ax)

    # One cutflow counter per region
    regions = monojet_regions(cfg).keys()
    for region in regions:
        if region=="inclusive":
            continue
        items[f'cutflow_{region}']  = processor.defaultdict_accumulator(accu_int)

    items['nevents'] = processor.defaultdict_accumulator(float)
    items['sumw'] = processor.defaultdict_accumulator(float)
    items['sumw2'] = processor.defaultdict_accumulator(float)
    items['sumw_pileup'] = processor.defaultdict_accumulator(float)

    items['selected_events'] = processor.defaultdict_accumulator(empty_column_accumulator_int)
    items['kinematics'] = processor.defaultdict_accumulator(list)


    items['tree_bool'] = processor.defaultdict_accumulator(defaultdict_accumulator_of_empty_column_accumulator_bool)
    items['tree_int64'] = processor.defaultdict_accumulator(defaultdict_accumulator_of_empty_column_accumulator_int64)
    items['tree_float16'] = processor.defaultdict_accumulator(defaultdict_accumulator_of_empty_column_accumulator_float16)

    items['weights'] = Hist("Weights", dataset_ax, region_ax, weight_type_ax, weight_ax)
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
        # 2016, 2017 data
        jes_suffix = ''
        jes_suffix_met = ''
    elif df['is_data']:
        # 2018 data
        jes_suffix = '_nom'
        jes_suffix_met = '_nom'
    else:
        # MC, all years
        jes_suffix = '_nom'
        if cfg.MET.JER:
            jes_suffix_met = '_jer'
        else:
            jes_suffix_met = '_nom'

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
        etasc=df['Electron_eta']+df['Electron_deltaEtaSC'],
        absetasc=np.abs(df['Electron_eta']+df['Electron_deltaEtaSC']),
        phi=df['Electron_phi'],
        mass=0 * df['Electron_pt'],
        charge=df['Electron_charge'],
        looseId=(df[cfg.ELECTRON.BRANCH.ID]>=1),
        tightId=(df[cfg.ELECTRON.BRANCH.ID]==4),
        dxy=np.abs(df['Electron_dxy']),
        dz=np.abs(df['Electron_dz']),
        barrel=np.abs(df['Electron_eta']+df['Electron_deltaEtaSC']) <= 1.4442
    )
    # All electrons must be at least loose
    pass_dxy = (electrons.barrel & (np.abs(electrons.dxy) < cfg.ELECTRON.CUTS.LOOSE.DXY.BARREL)) \
    | (~electrons.barrel & (np.abs(electrons.dxy) < cfg.ELECTRON.CUTS.LOOSE.DXY.ENDCAP))

    pass_dz = (electrons.barrel & (np.abs(electrons.dz) < cfg.ELECTRON.CUTS.LOOSE.DZ.BARREL)) \
    | (~electrons.barrel & (np.abs(electrons.dz) < cfg.ELECTRON.CUTS.LOOSE.DZ.ENDCAP))

    electrons = electrons[electrons.looseId \
                                    & (electrons.pt>cfg.ELECTRON.CUTS.LOOSE.PT) \
                                    & (electrons.absetasc<cfg.ELECTRON.CUTS.LOOSE.ETA) \
                                    & pass_dxy \
                                    & pass_dz
                                    ]

    if cfg.OVERLAP.ELECTRON.MUON.CLEAN:
        electrons = electrons[object_overlap(electrons, muons, dr=cfg.OVERLAP.ELECTRON.MUON.DR)]

    taus = JaggedCandidateArray.candidatesfromcounts(
        df['nTau'],
        pt=df['Tau_pt'],
        eta=df['Tau_eta'],
        abseta=np.abs(df['Tau_eta']),
        phi=df['Tau_phi'],
        mass=0 * df['Tau_pt'],
        decaymode=df[cfg.TAU.BRANCH.ID],
        iso=df[cfg.TAU.BRANCH.ISO]
    )

    # For MC, add the matched gen-particle info for checking
    if not df['is_data']:
        kwargs = {'genpartflav' : df['Tau_genPartFlav']}
        taus.add_attributes(**kwargs)

    taus = taus[ (taus.decaymode) \
                & (taus.pt > cfg.TAU.CUTS.PT)\
                & (taus.abseta < cfg.TAU.CUTS.ETA) \
                & ((taus.iso&2)==2)]

    if cfg.OVERLAP.TAU.MUON.CLEAN:
        taus = taus[object_overlap(taus, muons, dr=cfg.OVERLAP.TAU.MUON.DR)]
    if cfg.OVERLAP.TAU.ELECTRON.CLEAN:
        taus = taus[object_overlap(taus, electrons, dr=cfg.OVERLAP.TAU.ELECTRON.DR)]

    # choose the right branch name for photon ID bitmap depending on the actual name in the file (different between nano v5 and v7)
    if cfg.PHOTON.BRANCH.ID in df.keys():
        PHOTON_BRANCH_ID = cfg.PHOTON.BRANCH.ID
    else:
        PHOTON_BRANCH_ID = cfg.PHOTON.BRANCH.IDV7
    photons = JaggedCandidateArray.candidatesfromcounts(
        df['nPhoton'],
        pt=df['Photon_pt'],
        eta=df['Photon_eta'],
        abseta=np.abs(df['Photon_eta']),
        phi=df['Photon_phi'],
        mass=0*df['Photon_pt'],
        looseId= (df[PHOTON_BRANCH_ID]>=1) & df['Photon_electronVeto'],
        mediumId=(df[PHOTON_BRANCH_ID]>=2) & df['Photon_electronVeto'],
        r9=df['Photon_r9'],
        barrel=df['Photon_isScEtaEB'],
    )
    photons = photons[photons.looseId \
              & (photons.pt > cfg.PHOTON.CUTS.LOOSE.pt) \
              & (photons.abseta < cfg.PHOTON.CUTS.LOOSE.eta)
              ]

    if cfg.OVERLAP.PHOTON.MUON.CLEAN:
        photons = photons[object_overlap(photons, muons, dr=cfg.OVERLAP.PHOTON.MUON.DR)]
    if cfg.OVERLAP.PHOTON.ELECTRON.CLEAN:
        photons = photons[object_overlap(photons, electrons, dr=cfg.OVERLAP.PHOTON.ELECTRON.DR)]

    ak4 = JaggedCandidateArray.candidatesfromcounts(
        df['nJet'],
        pt=df[f'Jet_pt{jes_suffix}'] if (df['is_data'] or cfg.AK4.JER) else df[f'Jet_pt{jes_suffix}']/df['Jet_corr_JER'],
        eta=df['Jet_eta'],
        abseta=np.abs(df['Jet_eta']),
        phi=df['Jet_phi'],
        mass=np.zeros_like(df['Jet_pt']),
        looseId=(df['Jet_jetId']&2) == 2, # bitmask: 1 = loose, 2 = tight, 3 = tight + lep veto
        tightId=(df['Jet_jetId']&2) == 2, # bitmask: 1 = loose, 2 = tight, 3 = tight + lep veto
        puid=((df['Jet_puId']&2>0) | ((df[f'Jet_pt{jes_suffix}'] if (df['is_data'] or cfg.AK4.JER) else df[f'Jet_pt{jes_suffix}']/df['Jet_corr_JER'])>50)), # medium pileup jet ID
        csvv2=df["Jet_btagCSVV2"],
        deepcsv=df['Jet_btagDeepB'],
        nef=df['Jet_neEmEF'],
        nhf=df['Jet_neHEF'],
        chf=df['Jet_chHEF'],
        cef=df['Jet_chEmEF'],
        muf=df['Jet_muEF'],
        ptraw=df['Jet_pt']*(1-df['Jet_rawFactor']),
        nconst=df['Jet_nConstituents'],
        hadflav= 0*df['Jet_pt'] if df['is_data'] else df['Jet_hadronFlavour']
    )

    # Before cleaning, apply HEM veto
    hem_ak4 = ak4[ (ak4.pt>30) &
        (-3.0 < ak4.eta) &
        (ak4.eta < -1.3) &
        (-1.57 < ak4.phi) &
        (ak4.phi < -0.87)
        ]
    df['hemveto'] = hem_ak4.counts == 0

    # B jets have their own overlap cleaning,
    # so deal with them before applying filtering to jets
    btag_discriminator = getattr(ak4, cfg.BTAG.algo)
    btag_cut = cfg.BTAG.CUTS[cfg.BTAG.algo][cfg.BTAG.wp]
    bjets = ak4[
        (ak4.pt > cfg.BTAG.PT) \
        & (ak4.abseta < cfg.BTAG.ETA) \
        & (btag_discriminator > btag_cut)
    ]

    if cfg.OVERLAP.BTAG.MUON.CLEAN:
        bjets = bjets[object_overlap(bjets, muons, dr=cfg.OVERLAP.BTAG.MUON.DR)]
    if cfg.OVERLAP.BTAG.ELECTRON.CLEAN:
        bjets = bjets[object_overlap(bjets, electrons, dr=cfg.OVERLAP.BTAG.ELECTRON.DR)]
    if cfg.OVERLAP.BTAG.PHOTON.CLEAN:
        bjets = bjets[object_overlap(bjets, photons, dr=cfg.OVERLAP.BTAG.PHOTON.DR)]

    ak4 = ak4[ak4.looseId]

    if cfg.OVERLAP.AK4.MUON.CLEAN:
        ak4 = ak4[object_overlap(ak4, muons, dr=cfg.OVERLAP.AK4.MUON.DR)]
    if cfg.OVERLAP.AK4.ELECTRON.CLEAN:
        ak4 = ak4[object_overlap(ak4, electrons, dr=cfg.OVERLAP.AK4.ELECTRON.DR)]
    if cfg.OVERLAP.AK4.PHOTON.CLEAN:
        ak4 = ak4[object_overlap(ak4, photons, dr=cfg.OVERLAP.AK4.PHOTON.DR)]


    if df['is_data']:
        msd = df[f'FatJet_msoftdrop{jes_suffix}']
    else:
        msd = df[f'FatJet_msoftdrop{jes_suffix}'] / (df['FatJet_msoftdrop_corr_JMR'] * df['FatJet_msoftdrop_corr_JMS'])
        if not cfg.AK8.JER:
            msd = msd / df['FatJet_corr_JER']

    ak8 = JaggedCandidateArray.candidatesfromcounts(
        df['nFatJet'],
        pt=df[f'FatJet_pt{jes_suffix}'] if (df['is_data'] or cfg.AK8.JER) else df[f'FatJet_pt{jes_suffix}']/df['FatJet_corr_JER'],
        eta=df['FatJet_eta'],
        abseta=np.abs(df['FatJet_eta']),
        phi=df['FatJet_phi'],
        mass=msd,
        tightId=(df['FatJet_jetId']&2) == 2, # Tight
        csvv2=df["FatJet_btagCSVV2"],
        deepcsv=df['FatJet_btagDeepB'],
        tau1=df['FatJet_tau1'],
        tau2=df['FatJet_tau2'],
        tau21=df['FatJet_tau2']/df['FatJet_tau1'],
        wvsqcd=df['FatJet_deepTag_WvsQCD'],
        wvsqcdmd=df['FatJet_deepTagMD_WvsQCD'],
        zvsqcd=df['FatJet_deepTag_ZvsQCD'],
        zvsqcdmd=df['FatJet_deepTagMD_ZvsQCD'],
        tvsqcd=df['FatJet_deepTag_TvsQCD'],
        tvsqcdmd=df['FatJet_deepTagMD_TvsQCD'],
        wvstqcd=df['FatJet_deepTag_WvsQCD']*(1-df['FatJet_deepTag_TvsQCD'])/(1-df['FatJet_deepTag_WvsQCD']*df['FatJet_deepTag_TvsQCD']),
        wvstqcdmd=df['FatJet_deepTagMD_WvsQCD']*(1-df['FatJet_deepTagMD_TvsQCD'])/(1-df['FatJet_deepTagMD_WvsQCD']*df['FatJet_deepTagMD_TvsQCD']),
    )
    ak8 = ak8[ak8.tightId & object_overlap(ak8, muons) & object_overlap(ak8, electrons) & object_overlap(ak8, photons)]

    if extract_year(df['dataset']) == 2017:
        met_branch = 'METFixEE2017'
    else:
        met_branch = 'MET'

    met_pt = df[f'{met_branch}_pt{jes_suffix_met}']
    met_phi = df[f'{met_branch}_phi{jes_suffix_met}']

    return met_pt, met_phi, ak4, bjets, ak8, muons, electrons, taus, photons

def monojet_regions(cfg):
    common_cuts = [
        'filt_met',
        'veto_ele',
        'veto_muo',
        'veto_photon',
        'veto_tau',
        'veto_b',
        'mindphijr',
        'dpfcalo',
        'recoil',
        'hemveto'
    ]
    j_cuts = [
        'leadak4_pt_eta',
        'leadak4_id',
        'veto_vtag'
    ]
    # Test out different working point for v tagging
    # the first one is the traditional one used in 2016
    v_cuts = [
        'leadak8_pt_eta',
        'leadak8_id',
        'leadak8_mass',
        'leadak8_tau21',
    ]

    regions = {}
    regions['inclusive'] = ['inclusive']
    # Signal regions (v = mono-V, j = mono-jet)
    regions['sr_v'] = ['trig_met','hemveto_metphi'] + common_cuts + v_cuts + ['prescale']
    regions['sr_j'] = ['trig_met','hemveto_metphi'] + common_cuts + j_cuts + ['prescale']

    # regions['cr_nofilt_j'] = copy.deepcopy(regions['sr_j'])
    # regions['cr_nofilt_j'].remove('filt_met')

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

    # DeepAK8 WPs for signal regions with different level of purities
    wp_to_run = ['loose', 'tight']
    # and more WPs when doing the mistag measurement
    if cfg and cfg.RUN.MONOVMISTAG_STUDY:
        wp_to_run += ['medium', 'inclusive']

    for region in ['sr_v','cr_2m_v','cr_1m_v','cr_2e_v','cr_1e_v','cr_g_v']:
        for wp in wp_to_run:
            # the new region name will be, for example, cr_2m_loose_v
            newRegionName=region.replace('_v','_'+wp+'_v')
            regions[newRegionName] = copy.deepcopy(regions[region])
            regions[newRegionName].remove('leadak8_tau21')
            if wp == 'inclusive':
                if cfg and cfg.RUN.MONOVMISTAG_STUDY:
                    regions[newRegionName].remove('leadak8_mass')

                    # add regions: cr_2m_hasmass_inclusive_v, cr_1m_hasmass_inclusive_v, cr_2e_hasmass_inclusive_v, cr_1e_hasmass_inclusive_v
                    # these are regions with mass cut but has no tagger cut
                    hasMassRegionName = region.replace('_v', '_hasmass_'+ wp + '_v')
                    regions[hasMassRegionName] = regions[newRegionName] + ['leadak8_mass']

            else:
                regions[newRegionName].append('leadak8_wvsqcd_'+wp)
            # save a copy of the v-tagged regions but not applying mistag SFs, for the sake of measuring mistag SF later

            if cfg and cfg.RUN.MONOVMISTAG_STUDY:
                if wp in ['loose','medium','tight']:
                    noMistagRegionName = region.replace('_v', '_nomistag_'+ wp + '_v')
                    regions[noMistagRegionName]=copy.deepcopy(regions[newRegionName])
                    noMistagNoMassRegionName = region.replace('_v', '_nomistag_nomass_'+ wp + '_v')
                    regions[noMistagNoMassRegionName]=copy.deepcopy(regions[newRegionName])
                    regions[noMistagNoMassRegionName].remove('leadak8_mass')

    # Veto weight region
    tmp = {}
    for region in regions.keys():
        if not region.startswith("sr_"):
            continue

        if cfg and cfg.RUN.VETO_STUDY:
            new_region = f"{region}_no_veto_ele"
            tmp[new_region] = copy.deepcopy(regions[region])
            tmp[new_region].remove("veto_ele")
            tmp[new_region].remove("mindphijr")
            tmp[new_region].remove("recoil")
            tmp[new_region].remove("dpfcalo")
            tmp[new_region].append("met_sr")
            tmp[new_region].append("mindphijm")
            tmp[new_region].append("dpfcalo_sr")

            new_region = f"{region}_no_veto_tau"
            tmp[new_region] = copy.deepcopy(regions[region])
            tmp[new_region].remove("veto_tau")
            tmp[new_region].remove("mindphijr")
            tmp[new_region].remove("recoil")
            tmp[new_region].remove("dpfcalo")
            tmp[new_region].append("met_sr")
            tmp[new_region].append("mindphijm")
            tmp[new_region].append("dpfcalo_sr")

            new_region = f"{region}_no_veto_muon"
            tmp[new_region] = copy.deepcopy(regions[region])
            tmp[new_region].remove("veto_muo")
            tmp[new_region].remove("mindphijr")
            tmp[new_region].remove("recoil")
            tmp[new_region].remove("dpfcalo")
            tmp[new_region].append("met_sr")
            tmp[new_region].append("mindphijm")
            tmp[new_region].append("dpfcalo_sr")

        new_region = f"{region}_no_veto_all"
        tmp[new_region] = copy.deepcopy(regions[region])
        tmp[new_region].remove("veto_muo")
        tmp[new_region].remove("veto_tau")
        tmp[new_region].remove("veto_ele")
        tmp[new_region].remove("mindphijr")
        tmp[new_region].remove("recoil")
        tmp[new_region].remove("dpfcalo")
        tmp[new_region].append("met_sr")
        tmp[new_region].append("mindphijm")
        tmp[new_region].append("dpfcalo_sr")

    regions.update(tmp)

    tmp = {}
    for region in regions.keys():
        if not re.match("sr_(loose|tight)_v.*",region):
            continue
        new_region = f"{region}_no_lowmass_ak8"
        tmp[new_region] = copy.deepcopy(regions[region])
        tmp[new_region].append("no_lowmass_ak8")

        new_region = f"{region}_lowmass_ak8"
        tmp[new_region] = copy.deepcopy(regions[region])
        tmp[new_region].append("lowmass_ak8")

        new_region = f"{region}_vlowmass_ak8"
        tmp[new_region] = copy.deepcopy(regions[region])
        tmp[new_region].append("vlowmass_ak8")
    regions.update(tmp)

    tmp = {}
    for region in regions.keys():
        if not re.match("sr_((loose|tight)_v|j)(_no_veto_all)?$",region):
            continue
        new_region = f"{region}_dphipftkveto"
        tmp[new_region] = copy.deepcopy(regions[region])
        tmp[new_region].append("dphipftkveto")

        new_region = f"{region}_dphipftkvetoinv"
        tmp[new_region] = copy.deepcopy(regions[region])
        tmp[new_region].append("dphipftkvetoinv")

    regions.update(tmp)

    # tmp = {}
    # for region in regions.keys():
    #     if not re.match("sr_j(_no_veto_all)?",region):
    #         continue
    #     new_region = f"{region}_no_overlap_removal"
    #     tmp[new_region] = copy.deepcopy(regions[region])
    #     tmp[new_region].remove("veto_vtag")
    # regions.update(tmp)

    if cfg and cfg.RUN.TRIGGER_STUDY:
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

    if cfg and not cfg.RUN.MONOV:
        keys_to_remove = [ x for x in regions.keys() if x.endswith('_v') or '_v_' in x]
        for key in keys_to_remove:
            regions.pop(key)
    if cfg and not cfg.RUN.MONOJ:
        keys_to_remove = [ x for x in regions.keys() if x.endswith('_j') or '_j_' in x]
        for key in keys_to_remove:
            regions.pop(key)

    return regions


def fitfun(x, a, b, c):
    return a * np.exp(-b * x) + c

def theory_weights_monojet(weights, df, evaluator, gen_v_pt, gen_ak8_mass):
    weights = processor.Weights(size=df.size, storeIndividual=True)

    # Guard against missing masses
    mass_valid = gen_ak8_mass > 0
    gen_ak8_mass[~mass_valid] = 85

    # Default: flat 1
    qcd_nlo_j         = np.ones(df.size)
    qcd_nlo_j_central = qcd_nlo_j
    qcd_nlo_v         = qcd_nlo_j
    qcd_nlo_v_central = qcd_nlo_j
    ewk_nlo           = qcd_nlo_j

    if df['is_lo_w']:
        qcd_nlo_j          = evaluator["qcd_nlo_w_j"](gen_v_pt)
        qcd_nlo_j_central  = evaluator["qcd_nlo_w_j_central"](gen_v_pt)
        qcd_nlo_v          = evaluator["qcd_nlo_w_v"](gen_v_pt, gen_ak8_mass)
        qcd_nlo_v_central  = evaluator["qcd_nlo_w_v_central"](gen_v_pt, gen_ak8_mass)
        ewk_nlo = evaluator["ewk_nlo_w"](gen_v_pt)
    elif df['is_lo_z']:
        if df['is_lo_znunu']:
            qcd_nlo_j          = evaluator["qcd_nlo_znn_j"](gen_v_pt)
            qcd_nlo_j_central  = qcd_nlo_j
            qcd_nlo_v          = evaluator["qcd_nlo_znn_v"](gen_v_pt, gen_ak8_mass)
            qcd_nlo_v_central  = qcd_nlo_v
        else:
            qcd_nlo_j          = evaluator["qcd_nlo_dy_j_central"](gen_v_pt)
            qcd_nlo_j_central  = qcd_nlo_j
            qcd_nlo_v          = evaluator["qcd_nlo_dy_v_central"](gen_v_pt, gen_ak8_mass)
            qcd_nlo_v_central  = qcd_nlo_v
        ewk_nlo = evaluator["ewk_nlo_z"](gen_v_pt)
    elif df['is_nlo_w']:
        ewk_nlo = evaluator["ewk_nlo_w"](gen_v_pt)
    elif df['is_nlo_z']:
        ewk_nlo = evaluator["ewk_nlo_z"](gen_v_pt)
    elif df['is_lo_g']:
        qcd_nlo_j = evaluator["qcd_nlo_g"](gen_v_pt)
        qcd_nlo_j_central = qcd_nlo_j
        qcd_nlo_v         = qcd_nlo_j
        qcd_nlo_v_central = qcd_nlo_j
        ewk_nlo = evaluator["ewk_nlo_g"](gen_v_pt)
    elif df['is_nlo_g']:
        ewk_nlo = evaluator["ewk_nlo_g"](gen_v_pt)
    else:
        pass


    # Guard against invalid input pt
    invalid_pt = (gen_v_pt <=0) | np.isinf(gen_v_pt) | np.isnan(gen_v_pt)
    qcd_nlo_j         = np.where(invalid_pt, 1, qcd_nlo_j)
    qcd_nlo_j_central = np.where(invalid_pt, 1, qcd_nlo_j_central)
    qcd_nlo_v         = np.where(invalid_pt, 1, qcd_nlo_v)
    qcd_nlo_v_central = np.where(invalid_pt, 1, qcd_nlo_v_central)


    # If anything is still invalid, that's an error
    is_invalid = lambda x: np.isinf(x) | np.isnan(x)
    assert(~np.any(is_invalid(qcd_nlo_j)))
    assert(~np.any(is_invalid(qcd_nlo_j_central)))
    assert(~np.any(is_invalid(qcd_nlo_v)))
    assert(~np.any(is_invalid(qcd_nlo_v_central)))
    assert(~np.any(is_invalid(ewk_nlo)))

    weights.add('sf_nlo_qcd_j',         qcd_nlo_j)
    weights.add('sf_nlo_qcd_j_central', qcd_nlo_j_central)
    weights.add('sf_nlo_qcd_v',         qcd_nlo_v)
    weights.add('sf_nlo_qcd_v_central', qcd_nlo_v_central)
    weights.add('sf_nlo_ewk',           ewk_nlo)
    return weights

def theory_weights_vbf(weights, df, evaluator, gen_v_pt, mjj):
    if df['is_lo_w']:
        theory_weights = evaluator["qcd_nlo_w_2017_2d"](mjj, gen_v_pt) * evaluator["ewk_nlo_w"](gen_v_pt)
    elif df['is_lo_w_ewk']:
        theory_weights = evaluator["qcd_nlo_w_ewk"](gen_v_pt, mjj)
    elif df['is_lo_z']:
        w_ewk = evaluator["ewk_nlo_z"](gen_v_pt)
        if df['is_lo_znunu']:
            w_qcd = evaluator["qcd_nlo_znn_2017_2d"](gen_v_pt, mjj)
        else:
            w_qcd = evaluator["qcd_nlo_z_2017_2d"](mjj, gen_v_pt)
        theory_weights = w_ewk * w_qcd
    elif df['is_lo_z_ewk']:
        theory_weights = evaluator["qcd_nlo_z_ewk"](gen_v_pt, mjj)
    elif df['is_nlo_w']:
        theory_weights = evaluator["ewk_nlo_w"](gen_v_pt)
    elif df['is_nlo_z']:
        theory_weights = evaluator["ewk_nlo_z"](gen_v_pt)
    elif df['is_lo_g']:
        theory_weights = evaluator["qcd_nlo_g_2017_2d"](mjj, gen_v_pt) * evaluator["ewk_nlo_g"](gen_v_pt)
    else:
        theory_weights = np.ones(df.size)

    # Guard against invalid input pt
    invalid = (gen_v_pt <=0) | np.isinf(gen_v_pt) | np.isnan(gen_v_pt)
    theory_weights[invalid] = 1

    weights.add('theory', theory_weights)

    return weights

def pileup_weights(weights, df, evaluator, cfg):

    if cfg.SF.PILEUP.MODE == 'nano':
        pu_weight = df['puWeight']
    elif cfg.SF.PILEUP.MODE == 'manual':
        pu_weight = evaluator['pileup'](df['Pileup_nTrueInt'])
    else:
        raise RuntimeError(f"Unknown value for cfg.PILEUP.MODE: {cfg.PILEUP.MODE}.")

    # Cap weights just in case
    pu_weight[np.abs(pu_weight)>5] = 1
    weights.add("pileup", pu_weight)
    return weights

def photon_trigger_sf(weights, photons, df):
    """MC-to-data photon trigger scale factor.

    The scale factor is obtained by separately fitting the
    trigger turn-on with a sigmoid function in data and MC.
    The scale factor is then the ratio of the two sigmoid
    functions as a function of the photon pt.

    :param weights: Weights object to write information into
    :type weights: WeightsContainer
    :param photons: Photon candidates
    :type photons: JaggedCandidateArray
    :param df: Data frame
    :type df: LazyDataFrame
    """
    year = extract_year(df['dataset'])
    x = photons.pt.max()
    if year == 2016:
        sf =  np.ones(df.size)
    elif year == 2017:
        sf = sigmoid(x,0.335,217.91,0.065,0.996) / sigmoid(x,0.244,212.34,0.050,1.000)
    elif year == 2018:
        sf = sigmoid(x,1.022, 218.39, 0.086, 0.999) / sigmoid(x, 0.301,212.83,0.062,1.000)

    sf[np.isnan(sf) | np.isinf(sf)] == 1

    weights.add("trigger_photon", sf)

def candidate_weights(weights, df, evaluator, muons, electrons, photons, cfg):
    year = extract_year(df['dataset'])
    # Muon ID and Isolation for tight and loose WP
    # Function of pT, eta (Order!)
    weight_muons_id_tight = evaluator['muon_id_tight'](muons[df['is_tight_muon']].pt, muons[df['is_tight_muon']].abseta).prod()
    weight_muons_iso_tight = evaluator['muon_iso_tight'](muons[df['is_tight_muon']].pt, muons[df['is_tight_muon']].abseta).prod()

    if cfg.SF.DIMUO_ID_SF.USE_AVERAGE:
        tight_dimuons = muons[df["is_tight_muon"]].distincts()
        t0 = (evaluator['muon_id_tight'](tight_dimuons.i0.pt, tight_dimuons.i0.abseta) \
             * evaluator['muon_iso_tight'](tight_dimuons.i0.pt, tight_dimuons.i0.abseta)).prod()
        t1 = (evaluator['muon_id_tight'](tight_dimuons.i1.pt, tight_dimuons.i1.abseta) \
             * evaluator['muon_iso_tight'](tight_dimuons.i1.pt, tight_dimuons.i1.abseta)).prod()
        l0 = (evaluator['muon_id_loose'](tight_dimuons.i0.pt, tight_dimuons.i0.abseta) \
             * evaluator['muon_iso_loose'](tight_dimuons.i0.pt, tight_dimuons.i0.abseta)).prod()
        l1 = (evaluator['muon_id_loose'](tight_dimuons.i1.pt, tight_dimuons.i1.abseta) \
             * evaluator['muon_iso_loose'](tight_dimuons.i1.pt, tight_dimuons.i1.abseta)).prod()
        weights_2m_tight = 0.5*( l0 * t1 + l1 * t0)
        weights.add("muon_id_iso_tight", weight_muons_id_tight*weight_muons_iso_tight*(tight_dimuons.counts!=1) + weights_2m_tight*(tight_dimuons.counts==1))
    else:
        weights.add("muon_id_iso_tight", weight_muons_id_tight*weight_muons_iso_tight )

    weights.add("muon_id_loose", evaluator['muon_id_loose'](muons[~df['is_tight_muon']].pt, muons[~df['is_tight_muon']].abseta).prod())
    weights.add("muon_iso_loose", evaluator['muon_iso_loose'](muons[~df['is_tight_muon']].pt, muons[~df['is_tight_muon']].abseta).prod())

    # Electron ID and reco
    # Function of eta, pT (Other way round relative to muons!)

    # For 2017, the reco SF is split below/above 20 GeV
    if year == 2017:
        high_et = electrons.pt>20
        ele_reco_sf = evaluator['ele_reco'](electrons.etasc[high_et], electrons.pt[high_et]).prod()
        ele_reco_sf *= evaluator['ele_reco_pt_lt_20'](electrons.etasc[~high_et], electrons.pt[~high_et]).prod()
    else:
        ele_reco_sf = evaluator['ele_reco'](electrons.etasc, electrons.pt).prod()
    weights.add("ele_reco", ele_reco_sf)
    # ID/iso SF is not split
    # in case of 2 tight electrons, we want to apply 0.5*(T1L2+T2L1) instead of T1T2
    weights_electrons_tight = evaluator['ele_id_tight'](electrons[df['is_tight_electron']].etasc, electrons[df['is_tight_electron']].pt).prod()
    if cfg.SF.DIELE_ID_SF.USE_AVERAGE:
        tight_dielectrons = electrons[df["is_tight_electron"]].distincts()
        l0 = evaluator['ele_id_loose'](tight_dielectrons.i0.etasc, tight_dielectrons.i0.pt).prod()
        t0 = evaluator['ele_id_tight'](tight_dielectrons.i0.etasc, tight_dielectrons.i0.pt).prod()
        l1 = evaluator['ele_id_loose'](tight_dielectrons.i1.etasc, tight_dielectrons.i1.pt).prod()
        t1 = evaluator['ele_id_tight'](tight_dielectrons.i1.etasc, tight_dielectrons.i1.pt).prod()
        weights_2e_tight = 0.5*( l0 * t1 + l1 * t0)
        weights.add("ele_id_tight", weights_electrons_tight*(tight_dielectrons.counts!=1) + weights_2e_tight*(tight_dielectrons.counts==1))
    else:
        weights.add("ele_id_tight", weights_electrons_tight)
    weights.add("ele_id_loose", evaluator['ele_id_loose'](electrons[~df['is_tight_electron']].etasc, electrons[~df['is_tight_electron']].pt).prod())

    # Photon ID and electron veto
    if cfg.SF.PHOTON.USETNP:
        weights.add("photon_id_tight", evaluator['photon_id_tight_tnp'](np.abs(photons[df['is_tight_photon']].eta)).prod())
    else:
        weights.add("photon_id_tight", evaluator['photon_id_tight'](photons[df['is_tight_photon']].eta, photons[df['is_tight_photon']].pt).prod())

    if year == 2016:
        csev_weight = evaluator["photon_csev"](photons.abseta, photons.pt).prod()
    elif year == 2017:
        csev_sf_index = 0.5 * photons.barrel + 3.5 * ~photons.barrel + 1 * (photons.r9 > 0.94) + 2 * (photons.r9 <= 0.94)
        csev_weight = evaluator['photon_csev'](csev_sf_index).prod()
    elif year == 2018:
        csev_weight = evaluator['photon_csev'](photons.pt, photons.abseta).prod()
    csev_weight[csev_weight==0] = 1
    weights.add("photon_csev", csev_weight)

    return weights

def data_driven_qcd_dataset(dataset):
    """Dataset name to use for data-driven QCD estimate"""
    year = extract_year(dataset)
    return f"QCD_data_{year}"

def photon_impurity_weights(photon_pt, year):
    """Photon impurity as a function of pt

    :param photon_pt: Photon pt
    :type photon_pt: 1D array
    :param year: Data-taking year
    :type year: int
    :return: Weights
    :rtype: 1Darray
    """
    if year == 2017:
        a = 6.35
        b = 4.61e-3
        c = 1.05
    elif year==2018:
        a = 11.92
        b = 8.28e-3
        c = 1.55
    elif year==2016:
        return np.ones(photon_pt.size)

    # Remember to multiply by 1e-2,
    # because fit is done to percentages
    return 1e-2*exponential(photon_pt, a, b, c)
