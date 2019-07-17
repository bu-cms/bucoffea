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

def monojet_accumulator():
    dataset_ax = Cat("dataset", "Primary dataset")
    region_ax = Cat("region", "Selection region")
    type_ax = Cat("type", "Type")

    vpt_ax = Bin("vpt",r"$p_{T}^{V}$ (GeV)", 100, 0, 2000)

    met_ax = Bin("met", r"$p_{T}^{miss}$ (GeV)", 200, 0, 2000)
    recoil_ax = Bin("recoil", r"Recoil (GeV)", 200, 0, 2000)

    jet_pt_ax = Bin("jetpt", r"$p_{T}$ (GeV)", 100, 0, 1000)
    jet_eta_ax = Bin("jeteta", r"$\eta$ (GeV)", 50, -5, 5)

    jet_mass_ax = Bin("mass", r"$M_{jet}$ (GeV)", 100,0,300)

    dpfcalo_ax = Bin("dpfcalo", r"$1-Calo/PF$", 20, -1, 1)
    btag_ax = Bin("btag", r"B tag discriminator", 20, 0, 1)
    multiplicity_ax = Bin("multiplicity", r"multiplicity", 10, -0.5, 9.5)
    dphi_ax = Bin("dphi", r"$\Delta\phi$", 50, 0, 3.5)

    pt_ax = Bin("pt", r"$p_{T}$ (GeV)", 100, 0, 1000)
    mt_ax = Bin("mt", r"$M_{T}$ (GeV)", 100, 0, 1000)
    eta_ax = Bin("eta", r"$\eta$ (GeV)", 50, -5, 5)
    dilepton_mass_ax = Bin("dilepton_mass", r"$M(\ell\ell)$ (GeV)", 100,50,150)

    weight_type_ax = Cat("weight_type", "Weight type")
    weight_ax = Bin("weight_value", "Weight",100,0.5,1.5)

    Hist = hist.Hist
    items = {}
    items["genvpt_check"] = Hist("Counts", dataset_ax, type_ax, vpt_ax)
    items["met"] = Hist("Counts", dataset_ax, region_ax, met_ax)
    items["recoil"] = Hist("Counts", dataset_ax, region_ax, recoil_ax)

    items["ak4pt0"] = Hist("Counts", dataset_ax, region_ax, jet_pt_ax)
    items["ak4eta0"] = Hist("Counts", dataset_ax, region_ax, jet_eta_ax)
    items["ak4pt"] = Hist("Counts", dataset_ax, region_ax, jet_pt_ax)
    items["ak4eta"] = Hist("Counts", dataset_ax, region_ax, jet_eta_ax)
    items["ak4btag"] = Hist("Counts", dataset_ax, region_ax, btag_ax)

    items["ak8pt0"] = Hist("Counts", dataset_ax, region_ax, jet_pt_ax)
    items["ak8eta0"] = Hist("Counts", dataset_ax, region_ax, jet_eta_ax)
    items["ak8mass0"] = Hist("Counts", dataset_ax, region_ax, jet_mass_ax)
    items["ak8pt"] = Hist("Counts", dataset_ax, region_ax, jet_pt_ax)
    items["ak8eta"] = Hist("Counts", dataset_ax, region_ax, jet_eta_ax)
    items["ak8mass"] = Hist("Counts", dataset_ax, region_ax, jet_mass_ax)

    items["dpfcalo"] = Hist("Counts", dataset_ax, region_ax, dpfcalo_ax)
    items["dphijm"] = Hist("min(4 leading jets, MET)", dataset_ax, region_ax, dphi_ax)

    # Multiplicity histograms
    for cand in ['ak4', 'ak8', 'bjet', 'loose_ele', 'loose_muo', 'tight_ele', 'tight_muo', 'tau', 'photon']:
        items[f"{cand}_mult"] = Hist(cand, dataset_ax, region_ax, multiplicity_ax)

    items["muon_pt"] = Hist("Counts", dataset_ax, region_ax, pt_ax)
    items["muon_eta"] = Hist("Counts", dataset_ax, region_ax, eta_ax)
    items["muon_mt"] = Hist("Counts", dataset_ax, region_ax, mt_ax)
    items["dimuon_pt"] = Hist("Counts", dataset_ax, region_ax, pt_ax)
    items["dimuon_eta"] = Hist("Counts", dataset_ax, region_ax, eta_ax)
    items["dimuon_mass"] = Hist("Counts", dataset_ax, region_ax, dilepton_mass_ax)

    items["electron_pt"] = Hist("Counts", dataset_ax, region_ax, pt_ax)
    items["electron_eta"] = Hist("Counts", dataset_ax, region_ax, eta_ax)
    items["electron_mt"] = Hist("Counts", dataset_ax, region_ax, mt_ax)
    items["dielectron_pt"] = Hist("Counts", dataset_ax, region_ax, pt_ax)
    items["dielectron_eta"] = Hist("Counts", dataset_ax, region_ax, eta_ax)
    items["dielectron_mass"] = Hist("Counts", dataset_ax, region_ax, dilepton_mass_ax)

    # One cutflow counter per region
    regions = monojet_regions().keys()
    for region in regions:
        if region=="inclusive":
            continue
        items[f'cutflow_{region}']  = processor.defaultdict_accumulator(int)

    items['sumw'] = processor.defaultdict_accumulator(float)
    items['sumw2'] = processor.defaultdict_accumulator(float)

    items['selected_events'] = processor.defaultdict_accumulator(list)
    items['kinematics'] = processor.defaultdict_accumulator(list)

    items['weights'] = Hist("Weights", dataset_ax, region_ax, weight_type_ax, weight_ax)

    return  processor.dict_accumulator(items)

def setup_gen_candidates(df):
    # Find first ancestor with different PDG ID
    # before defining the gen candidates
    mothers = JaggedArray.fromcounts(df['nGenPart'],df['GenPart_genPartIdxMother'] )
    pdgids = JaggedArray.fromcounts(df['nGenPart'],df['GenPart_pdgId'] )
    parent_index =  find_first_parent(mothers, pdgids)

    gen = JaggedCandidateArray.candidatesfromcounts(
        df['nGenPart'],
        pt=df['GenPart_pt'],
        eta=df['GenPart_eta'],
        phi=df['GenPart_phi'],
        mass=df['GenPart_mass'],
        charge=df['GenPart_pdgId'],
        pdg=df['GenPart_pdgId'],
        status=df['GenPart_status'],
        flag = df['GenPart_statusFlags'],
        mother = df['GenPart_genPartIdxMother'],
        parentIndex=parent_index.flatten())
    return gen

def setup_candidates(df, cfg):
    muons = JaggedCandidateArray.candidatesfromcounts(
        df['nMuon'],
        pt=df['Muon_pt'],
        eta=df['Muon_eta'],
        phi=df['Muon_phi'],
        mass=0 * df['Muon_pt'],
        charge=df['Muon_charge'],
        looseId=df['Muon_looseId'],
        iso=df["Muon_pfRelIso04_all"],
        tightId=df['Muon_tightId']
    )

    # All muons must be at least loose
    muons = muons[muons.looseId \
                    & (muons.iso < cfg.MUON.CUTS.LOOSE.ISO) \
                    & (muons.pt > cfg.MUON.CUTS.LOOSE.PT) \
                    & (np.abs(muons.eta)<cfg.MUON.CUTS.LOOSE.ETA) \
                    ]


    electrons = JaggedCandidateArray.candidatesfromcounts(
        df['nElectron'],
        pt=df['Electron_pt'],
        eta=df['Electron_eta'],
        phi=df['Electron_phi'],
        mass=0 * df['Electron_pt'],
        charge=df['Electron_charge'],
        looseId=(df[cfg.ELECTRON.BRANCH.ID]>=1),
        tightId=(df[cfg.ELECTRON.BRANCH.ID]==4),
        dxy=np.abs(df['Electron_dxy']),
        dz=np.abs(df['Electron_dz']),
        barrel=np.abs(df['Electron_eta'])< 1.479
    )
    # All electrons must be at least loose
    pass_dxy = (electrons.barrel & (electrons.dxy < cfg.ELECTRON.CUTS.LOOSE.DXY.BARREL)) \
    | (~electrons.barrel & (electrons.dxy < cfg.ELECTRON.CUTS.LOOSE.DXY.ENDCAP))

    pass_dz = (electrons.barrel & (electrons.dz < cfg.ELECTRON.CUTS.LOOSE.DZ.BARREL)) \
    | (~electrons.barrel & (electrons.dz < cfg.ELECTRON.CUTS.LOOSE.DZ.ENDCAP))

    electrons = electrons[electrons.looseId \
                                    & (electrons.pt>cfg.ELECTRON.CUTS.LOOSE.PT) \
                                    & (np.abs(electrons.eta)<cfg.ELECTRON.CUTS.LOOSE.ETA) \
                                    & pass_dxy \
                                    & pass_dz
                                    ]


    taus = JaggedCandidateArray.candidatesfromcounts(
        df['nTau'],
        pt=df['Tau_pt'],
        eta=df['Tau_eta'],
        phi=df['Tau_phi'],
        mass=0 * df['Tau_pt'],
        decaymode=df['Tau_idDecayModeNewDMs'] & df['Tau_idDecayMode'],
        clean=df['Tau_cleanmask'],
        iso=df['Tau_idMVAnewDM2017v2'],
    )

    taus = taus[ object_overlap(taus, muons) \
                & object_overlap(taus, electrons) \
                & (taus.decaymode) \
                & (taus.pt > cfg.TAU.CUTS.PT)\
                & (np.abs(taus.eta) < cfg.TAU.CUTS.ETA) \
                & ((taus.iso&4)==4)]
    photons = JaggedCandidateArray.candidatesfromcounts(
        df['nPhoton'],
        pt=df['Photon_pt'],
        eta=df['Photon_eta'],
        phi=df['Photon_phi'],
        mass=0*df['Photon_pt'],
        looseId=(df[cfg.PHOTON.BRANCH.ID]>=1) & df['Photon_electronVeto'],
        mediumId=(df[cfg.PHOTON.BRANCH.ID]>=2) & df['Photon_electronVeto'],
        barrel=np.abs(df['Photon_eta']) < 1.479,
    )
    photons = photons[photons.looseId \
              & (photons.pt > cfg.PHOTON.CUTS.LOOSE.pt) \
              & (np.abs(photons.eta) < cfg.PHOTON.CUTS.LOOSE.eta)]

    ak4 = JaggedCandidateArray.candidatesfromcounts(
        df['nJet'],
        pt=df['Jet_pt'],
        eta=df['Jet_eta'],
        phi=df['Jet_phi'],
        mass=df['Jet_mass'],
        looseId=(df['Jet_jetId']&2) == 2, # bitmask: 1 = loose, 2 = tight
        tightId=(df['Jet_jetId']&2) == 2, # bitmask: 1 = loose, 2 = tight
        csvv2=df["Jet_btagCSVV2"],
        deepcsv=df['Jet_btagDeepB'],
        # nef=df['Jet_neEmEF'],
        nhf=df['Jet_neHEF'],
        chf=df['Jet_chHEF'],
        # clean=df['Jet_cleanmask']
        # cef=df['Jet_chEmEF'],
    )

    ak4 = ak4[ak4.looseId & object_overlap(ak4, muons) & object_overlap(ak4, electrons)]

    ak8 = JaggedCandidateArray.candidatesfromcounts(
    df['nFatJet'],
    pt=df['FatJet_pt'],
    eta=df['FatJet_eta'],
    phi=df['FatJet_phi'],
    mass=df['FatJet_msoftdrop'],

    tightId=(df['FatJet_jetId']&2) == 2, # Tight
    csvv2=df["FatJet_btagCSVV2"],
    deepcsv=df['FatJet_btagDeepB'],
    tau1=df['FatJet_tau1'],
    tau2=df['FatJet_tau2']
    )

    hlt = JaggedCandidateArray.candidatesfromcounts(
        df['nTrigObj'],
        pt=df['TrigObj_pt'],
        eta=df['TrigObj_eta'],
        phi=df['TrigObj_phi'],
        mass=0*df['TrigObj_pt'],
        id=df['TrigObj_id'],
        filter=df['TrigObj_filterBits']
    )
    return ak4, ak8, muons, electrons, taus, photons, hlt


def monojet_regions():
    common_cuts = [
        'filt_met',
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
    # v_cuts = [
        # 'leadak8_pt_eta',
        # 'leadak8_id',
        # 'leadak8_mass',
    # ]

    regions = {}
    regions['inclusive'] = ['inclusive']

    # Signal regions (v = mono-V, j = mono-jet)
    # regions['sr_v'] = ['trig_met'] + common_cuts + v_cuts
    regions['sr_j'] = ['trig_met'] + common_cuts + j_cuts

    # Dimuon CR
    cr_2m_cuts = ['trig_met','two_muons', 'at_least_one_tight_mu', 'dimuon_mass', 'dimuon_charge'] + common_cuts
    cr_2m_cuts.remove('veto_muo')

    regions['cr_2m_j'] = cr_2m_cuts + j_cuts
    # regions['cr_2m_v'] = cr_2m_cuts + v_cuts

    # Single muon CR
    cr_1m_cuts = ['trig_met','one_muon', 'at_least_one_tight_mu', 'mt_mu'] + common_cuts
    cr_1m_cuts.remove('veto_muo')
    regions['cr_1m_j'] = cr_1m_cuts + j_cuts
    # regions['cr_1m_v'] = cr_1m_cuts + v_cuts

    # Dielectron CR
    cr_2e_cuts = ['trig_ele','two_electrons', 'at_least_one_tight_el', 'dielectron_mass', 'dielectron_charge'] + common_cuts
    cr_2e_cuts.remove('veto_ele')
    regions['cr_2e_j'] = cr_2e_cuts + j_cuts
    # regions['cr_2e_v'] = cr_2e_cuts + v_cuts

    # Single electron CR
    cr_1e_cuts = ['trig_ele','one_electron', 'at_least_one_tight_el', 'mt_el'] + common_cuts
    cr_1e_cuts.remove('veto_ele')
    regions['cr_1e_j'] =  cr_1e_cuts + j_cuts
    # regions['cr_1e_v'] =  cr_1e_cuts + v_cuts

    # Photon CR
    cr_g_cuts = ['trig_photon', 'one_photon', 'at_least_one_tight_photon'] + common_cuts
    cr_g_cuts.remove('veto_photon')

    regions['cr_g_j'] = cr_g_cuts + j_cuts
    # regions['cr_g_v'] = cr_g_cuts + v_cuts

    # Trigger studies
    # num = numerator, den = denominator
    # Single Mu region: Remove recoil cut, add SingleMu trigger, toggle MET trigger
    tr_1m_num_cuts = copy.deepcopy(cr_1m_cuts) + j_cuts
    tr_1m_num_cuts.remove('recoil')
    tr_1m_num_cuts.append('trig_mu')
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
    tr_g_num_cuts = copy.deepcopy(cr_g_cuts)
    tr_g_num_cuts.remove('recoil')
    tr_g_num_cuts.append('trig_ht_for_g_eff')
    regions['tr_g_num_cuts'] = tr_g_num_cuts

    tr_g_den_cuts = copy.deepcopy(tr_g_num_cuts)
    tr_g_den_cuts.remove('trig_photon')
    regions['tr_g_den_cuts'] = tr_g_den_cuts

    return regions



from coffea.lookup_tools import extractor

def monojet_evaluator(cfg):
    """Initiates the SF evaluator and populates it with the right values

    :param cfg: Configuration
    :type cfg: DynaConf object
    :return: Ready-to-use SF evaluator
    :rtype: coffea.lookup_tools.evaluator
    """
    ext = extractor()

    for sfname, definition in cfg.SF.items():
        fpath = bucoffea_path(definition['file'])
        ext.add_weight_sets([f"{sfname} {definition['histogram']} {fpath}"])

    ext.finalize()

    evaluator = ext.make_evaluator()
    return evaluator
