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
from bucoffea.monojet.definitions import accu_int
from pprint import pprint

def vbfhinv_accumulator(cfg, variations):
    dataset_ax = Cat("dataset", "Primary dataset")
    unc_ax = Cat("uncertainty", "Uncertainty weight variation")
    variation_ax = Cat("variation", "Uncertainty weight variation")
    region_ax = Cat("region", "Selection region")
    type_ax = Cat("type", "Type")
    var_ax = Cat("var", "JES/JER variation type")

    vpt_ax = Bin("vpt",r"$p_{T}^{V}$ (GeV)", 100, 0, 2000)

    met_ax = Bin("met", r"$p_{T}^{miss}$ (GeV)", 200, 0, 2000)
    recoil_ax = Bin("recoil", r"Recoil (GeV)", 200, 0, 2000)

    ratio_ax = Bin("ratio", "Var/Nom - 1", 50, -1, 1)

    mjj_ax = Bin("mjj", r"$M_{jj}$ (GeV)", 150, 0, 7500)
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
    deta_ax = Bin("deta", r"$\Delta\eta$", 50, 0, 10)
    dr_ax = Bin("dr", r"$\Delta R$", 50, 0, 2)

    pt_ax = Bin("pt", r"$p_{T}$ (GeV)", 100, 0, 1000)
    pt_ax_mu = Bin("pt", r"$p_{T}$ (GeV)", [20,25,30,40,50,60,120])
    pt_ax_el = Bin("pt", r"$p_{T}$ (GeV)", [10,20,35,50,100,200,500])
    pt_ax_tau = Bin("pt", r"$p_{T}$ (GeV)", [18,20,25,30,35,40,500,1000])

    ht_ax = Bin("ht", r"$H_{T}$ (GeV)", 100, 0, 4000)
    mt_ax = Bin("mt", r"$M_{T}$ (GeV)", 100, 0, 1000)
    eta_ax = Bin("eta", r"$\eta$", 50, -5, 5)
    eta_ax_el = Bin("eta", r"$\eta$", [-2.5, -2.0, -1.56, -1.44, -0.8, 0, 0.8, 1.44,1.56,2.0,2.5])
    abseta_ax_mu = Bin("abseta", r"$|\eta|$", [0,0.9,1.2,2.1,2.4])

    eta_ax_coarse = Bin("eta", r"$\eta$", 25, -5, 5)
    phi_ax = Bin("phi", r"$\phi$", 50,-np.pi, np.pi)
    phi_ax_coarse = Bin("phi", r"$\phi$", 20,-np.pi, np.pi)

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

    # Multiplicity histograms
    for cand in ['ak4', 'ak8', 'bjet', 'loose_ele', 'loose_muo', 'tight_ele', 'tight_muo', 'tau', 'photon','hlt_single_muon','muons_hltmatch']:
        items[f"{cand}_mult"] = Hist(cand, dataset_ax, region_ax, multiplicity_ax)

    items["muon_pt"] = Hist("Counts", dataset_ax, region_ax, pt_ax)
    items["muon_pt_abseta"] = Hist("Counts", dataset_ax, region_ax, pt_ax_mu, abseta_ax_mu)
    items["muon_eta"] = Hist("Counts", dataset_ax, region_ax, eta_ax)
    items["muon_phi"] = Hist("Counts", dataset_ax, region_ax, phi_ax)
    items["muon_pt0"] = Hist("Counts", dataset_ax, region_ax, pt_ax)
    items["muon_eta0"] = Hist("Counts", dataset_ax, region_ax, eta_ax)
    items["muon_phi0"] = Hist("Counts", dataset_ax, region_ax, phi_ax)
    items["muon_pt1"] = Hist("Counts", dataset_ax, region_ax, pt_ax)
    items["muon_eta1"] = Hist("Counts", dataset_ax, region_ax, eta_ax)
    items["muon_phi1"] = Hist("Counts", dataset_ax, region_ax, phi_ax)

    items["dimuon_pt"] = Hist("Counts", dataset_ax, region_ax, pt_ax)
    items["dimuon_eta"] = Hist("Counts", dataset_ax, region_ax, eta_ax)
    items["dimuon_mass"] = Hist("Counts", dataset_ax, region_ax, dilepton_mass_ax)

    items["electron_pt"] = Hist("Counts", dataset_ax, region_ax, pt_ax)
    items["electron_pt_eta"] = Hist("Counts", dataset_ax, region_ax, pt_ax_el, eta_ax_el)
    items["electron_eta"] = Hist("Counts", dataset_ax, region_ax, eta_ax)
    items["electron_phi"] = Hist("Counts", dataset_ax, region_ax, phi_ax)
    items["electron_pt0"] = Hist("Counts", dataset_ax, region_ax, pt_ax)
    items["electron_eta0"] = Hist("Counts", dataset_ax, region_ax, eta_ax)
    items["electron_phi0"] = Hist("Counts", dataset_ax, region_ax, phi_ax)
    items["electron_pt1"] = Hist("Counts", dataset_ax, region_ax, pt_ax)
    items["electron_eta1"] = Hist("Counts", dataset_ax, region_ax, eta_ax)
    items["electron_phi1"] = Hist("Counts", dataset_ax, region_ax, phi_ax)

    items["dielectron_pt"] = Hist("Counts", dataset_ax, region_ax, pt_ax)
    items["dielectron_eta"] = Hist("Counts", dataset_ax, region_ax, eta_ax)
    items["dielectron_mass"] = Hist("Counts", dataset_ax, region_ax, dilepton_mass_ax)

    items['photon_pt0'] = Hist("Counts", dataset_ax, region_ax, pt_ax)
    items['photon_eta0'] = Hist("Counts", dataset_ax, region_ax, eta_ax)
    items['photon_phi0'] = Hist("Counts", dataset_ax, region_ax, phi_ax)

    items[f"met"] = Hist("Counts", dataset_ax, region_ax, met_ax, var_ax)
    items[f"recoil"] = Hist("Counts", dataset_ax, region_ax, recoil_ax, var_ax)
    items[f"mjj"] = Hist("Counts", dataset_ax, region_ax, mjj_ax, var_ax)
    
    items['recoil_varovernom'] = Hist('Counts', dataset_ax, region_ax, ratio_ax, var_ax)
    items['mjj_varovernom'] = Hist('Counts', dataset_ax, region_ax, ratio_ax, var_ax)
    items['ak4_pt0_varovernom'] = Hist('Counts', dataset_ax, region_ax, ratio_ax, var_ax)
    items['ak4_pt1_varovernom'] = Hist('Counts', dataset_ax, region_ax, ratio_ax, var_ax)
    items['detajj_varovernom'] = Hist('Counts', dataset_ax, region_ax, ratio_ax, var_ax)
    items['dphijj_varovernom'] = Hist('Counts', dataset_ax, region_ax, ratio_ax, var_ax)

    items[f"ak4_pt0"] = Hist("Counts", dataset_ax, region_ax, jet_pt_ax, var_ax)
    items[f"ak4_pt1"] = Hist("Counts", dataset_ax, region_ax, jet_pt_ax, var_ax)
    items[f"ak4_pt"] = Hist("Counts", dataset_ax, region_ax, jet_pt_ax, var_ax)

    items["tau_pt"] = Hist("Counts", dataset_ax, region_ax, pt_ax_tau)

    # One cutflow counter per region
    regions = vbfhinv_regions(cfg, variations).keys()
    for region in regions:
        if region=="inclusive":
            continue
        items[f'cutflow_{region}']  = processor.defaultdict_accumulator(accu_int)

    items['nevents'] = processor.defaultdict_accumulator(float)
    items['sumw'] = processor.defaultdict_accumulator(float)
    items['sumw2'] = processor.defaultdict_accumulator(float)
    items['sumw_pileup'] = processor.defaultdict_accumulator(float)

    items['selected_events'] = processor.defaultdict_accumulator(list)
    items['kinematics'] = processor.defaultdict_accumulator(list)

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

def vbfhinv_regions(cfg, variations):
    regions = {}
    regions['inclusive'] = ['inclusive']
    
    for var in variations:
        common_cuts = [
            'veto_ele',
            'veto_muo',
            'filt_met',
            'hemveto',
            'veto_photon',
            'veto_tau',
            f'veto_b{var}',
            f'mindphijr{var}',
            f'recoil{var}',
            f'two_jets{var}',
            f'leadak4_pt_eta{var}',
            f'leadak4_id{var}',
            f'trailak4_pt_eta{var}',
            f'trailak4_id{var}',
            f'hemisphere{var}',
            f'mjj{var}',
            f'dphijj{var}',
            f'detajj{var}'
        ]

        # Signal regions (v = mono-V, j = mono-jet)
        regions[f'sr_vbf{var}'] = ['trig_met'] + common_cuts

        # For sync mode
        if cfg.RUN.SYNC:
            regions['cr_sync'] = [
                'trig_met',
                'veto_photon',
                'mindphijr',
                'recoil',
                'two_jets',
                'leadak4_pt_eta',
                'leadak4_id',
                'trailak4_pt_eta',
                'trailak4_id',
                'hemisphere',
                'mjj',
                'dphijj',
                'detajj'
            ]

        # Dimuon CR
        cr_2m_cuts = ['trig_met','two_muons', 'at_least_one_tight_mu', 'dimuon_mass', 'veto_ele', 'dimuon_charge'] + common_cuts[1:]
        cr_2m_cuts.remove('veto_muo')

        regions[f'cr_2m_vbf{var}'] = cr_2m_cuts

        # Single muon CR
        cr_1m_cuts = ['trig_met','one_muon', 'at_least_one_tight_mu',  'veto_ele'] + common_cuts[1:]
        cr_1m_cuts.remove('veto_muo')
        regions[f'cr_1m_vbf{var}'] = cr_1m_cuts 

        # Dielectron CR
        cr_2e_cuts = ['trig_ele','two_electrons', 'at_least_one_tight_el', 'dielectron_mass', 'veto_muo', 'dielectron_charge'] + common_cuts[2:]
        # cr_2e_cuts.remove('veto_ele')
        regions[f'cr_2e_vbf{var}'] = cr_2e_cuts 

        # Single electron CR
        cr_1e_cuts = ['trig_ele','one_electron', 'at_least_one_tight_el', 'veto_muo',f'met_el{var}'] + common_cuts[1:]
        # cr_1e_cuts.remove('veto_ele')
        regions[f'cr_1e_vbf{var}'] =  cr_1e_cuts

        # Photon CR
        cr_g_cuts = ['trig_photon', 'one_photon', 'at_least_one_tight_photon','photon_pt'] + common_cuts
        cr_g_cuts.remove('veto_photon')

        regions[f'cr_g_vbf{var}'] = cr_g_cuts

    ##########################################

    if cfg.RUN.SYNC:
        regions['sync_sr_vbf_round1'] = [
                                        'filt_met',
                                        'trig_met',
                                        'veto_photon',
                                        'mindphijr',
                                        'recoil',
                                        'two_jets',
                                        'leadak4_pt_eta',
                                        'leadak4_id',
                                        'trailak4_pt_eta',
                                        'trailak4_id',
                                        'hemisphere',
                                        'mjj',
                                        'dphijj',
                                        'detajj',
                                        ]

    tmp = {}
    for region in regions.keys():
        if not region.startswith("sr_"):
            continue
        new_region = f"{region}_no_veto_all"
        tmp[new_region] = copy.deepcopy(regions[region])
        tmp[new_region].remove("veto_muo")
        tmp[new_region].remove("veto_tau")
        tmp[new_region].remove("veto_ele")
        tmp[new_region].remove("mindphijr")
        tmp[new_region].remove("recoil")
        tmp[new_region].append("met_sr")
        tmp[new_region].append("mindphijm")

    regions.update(tmp)

    if cfg.RUN.TRIGGER_STUDY:
        # Trigger studies
        # num = numerator, den = denominator
        # Single Mu region: Remove mjj cut, add SingleMu trigger, toggle MET trigger
        tr_1m_num_cuts = copy.deepcopy(cr_1m_cuts) 
        tr_1m_num_cuts.remove('mjj')
        tr_1m_num_cuts.append('trig_mu')
        tr_1m_num_cuts.append('mu_pt_trig_safe')

        regions['tr_1m_num_two_central_jets'] = tr_1m_num_cuts + ['two_central_jets']
        regions['tr_1m_num_two_forward_jets'] = tr_1m_num_cuts + ['two_forward_jets']
        regions['tr_1m_num_one_jet_forward_one_jet_central'] = tr_1m_num_cuts + ['one_jet_forward_one_jet_central']

        tr_1m_den_cuts = copy.deepcopy(tr_1m_num_cuts)
        tr_1m_den_cuts.remove('trig_met')

        regions['tr_1m_den_two_central_jets'] = tr_1m_den_cuts + ['two_central_jets']
        regions['tr_1m_den_two_forward_jets'] = tr_1m_den_cuts + ['two_forward_jets']
        regions['tr_1m_den_one_jet_forward_one_jet_central'] = tr_1m_den_cuts + ['one_jet_forward_one_jet_central']

        # Double Mu region: Remove mjj cut, toggle MET trigger
        tr_2m_num_cuts = copy.deepcopy(cr_2m_cuts) 
        tr_2m_num_cuts.remove('mjj')
        tr_2m_num_cuts.append('trig_mu')
        tr_2m_num_cuts.append('mu_pt_trig_safe')

        regions['tr_2m_num_two_central_jets'] = tr_2m_num_cuts + ['two_central_jets']
        regions['tr_2m_num_two_forward_jets'] = tr_2m_num_cuts + ['two_forward_jets']
        regions['tr_2m_num_one_jet_forward_one_jet_central'] = tr_2m_num_cuts + ['one_jet_forward_one_jet_central']

        tr_2m_den_cuts = copy.deepcopy(tr_2m_num_cuts)
        tr_2m_den_cuts.remove('trig_met')

        regions['tr_2m_den_two_central_jets'] = tr_2m_den_cuts + ['two_central_jets']
        regions['tr_2m_den_two_forward_jets'] = tr_2m_den_cuts + ['two_forward_jets']
        regions['tr_2m_den_one_jet_forward_one_jet_central'] = tr_2m_den_cuts + ['one_jet_forward_one_jet_central']

        # Single Electron region: Remove mjj cut, toggle MET trigger
        tr_1e_num_cuts = copy.deepcopy(cr_1e_cuts)
        tr_1e_num_cuts.remove('mjj')
        tr_1e_num_cuts.append('trig_met')

        regions['tr_1e_num_two_central_jets'] = tr_1e_num_cuts + ['two_central_jets']
        regions['tr_1e_num_two_forward_jets'] = tr_1e_num_cuts + ['two_forward_jets']
        regions['tr_1e_num_one_jet_forward_one_jet_central'] = tr_1e_num_cuts + ['one_jet_forward_one_jet_central']

        tr_1e_den_cuts = copy.deepcopy(tr_1e_num_cuts)
        tr_1e_den_cuts.remove('trig_met')

        regions['tr_1e_den_two_central_jets'] = tr_1e_den_cuts + ['two_central_jets']
        regions['tr_1e_den_two_forward_jets'] = tr_1e_den_cuts + ['two_forward_jets']
        regions['tr_1e_den_one_jet_forward_one_jet_central'] = tr_1e_den_cuts + ['one_jet_forward_one_jet_central']

        # Double Electron region: Remove mjj cut, toggle MET trigger
        tr_2e_num_cuts = copy.deepcopy(cr_2e_cuts)
        tr_2e_num_cuts.remove('mjj')
        tr_2e_num_cuts.append('trig_met')

        regions['tr_2e_num_two_central_jets'] = tr_2e_num_cuts + ['two_central_jets']
        regions['tr_2e_num_two_forward_jets'] = tr_2e_num_cuts + ['two_forward_jets']
        regions['tr_2e_num_one_jet_forward_one_jet_central'] = tr_2e_num_cuts + ['one_jet_forward_one_jet_central']

        tr_2e_den_cuts = copy.deepcopy(tr_2e_num_cuts)
        tr_2e_den_cuts.remove('trig_met')

        regions['tr_2e_den_two_central_jets'] = tr_2e_den_cuts + ['two_central_jets']
        regions['tr_2e_den_two_forward_jets'] = tr_2e_den_cuts + ['two_forward_jets']
        regions['tr_2e_den_one_jet_forward_one_jet_central'] = tr_2e_den_cuts + ['one_jet_forward_one_jet_central']
        

    return regions


