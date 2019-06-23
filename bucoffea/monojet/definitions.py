from coffea import hist
Hist = hist.Hist
Bin = hist.Bin
Cat = hist.Cat

from coffea.analysis_objects import JaggedCandidateArray
import coffea.processor as processor
from awkward import JaggedArray
import numpy as np

def monojet_accumulator():
    dataset_ax = Cat("dataset", "Primary dataset")
    region_ax = Cat("region", "Selection region")
    met_ax = Bin("met", r"$p_{T}^{miss}$ (GeV)", 100, 0, 1000)
    recoil_ax = Bin("recoil", r"Recoil (GeV)", 100, 0, 1000)
    jet_pt_ax = Bin("jetpt", r"$p_{T}$ (GeV)", 100, 0, 1000)
    jet_eta_ax = Bin("jeteta", r"$\eta$ (GeV)", 50, -5, 5)
    dpfcalo_ax = Bin("dpfcalo", r"$1-Calo/PF$", 20, -1, 1)
    btag_ax = Bin("btag", r"B tag discriminator", 20, 0, 1)
    multiplicity_ax = Bin("multiplicity", r"multiplicity", 10, -0.5, 9.5)
    dphi_ax = Bin("dphi", r"$\Delta\phi$", 50, 0, 3.5)

    pt_ax = Bin("pt", r"$p_{T}$ (GeV)", 100, 0, 1000)
    eta_ax = Bin("eta", r"$\eta$ (GeV)", 50, -5, 5)
    dilepton_mass_ax = Bin("dilepton_mass", r"$M(\ell\ell)$ (GeV)", 100,50,150)

    Hist = hist.Hist
    items = {}
    items["met"] = Hist("Counts", dataset_ax, region_ax, met_ax)
    items["recoil"] = Hist("Counts", dataset_ax, region_ax, recoil_ax)
    items["jet0pt"] = Hist("Counts", dataset_ax, region_ax, jet_pt_ax)
    items["jet0eta"] = Hist("Counts", dataset_ax, region_ax, jet_eta_ax)
    items["jetpt"] = Hist("Counts", dataset_ax, region_ax, jet_pt_ax)
    items["jeteta"] = Hist("Counts", dataset_ax, region_ax, jet_eta_ax)
    items["dpfcalo"] = Hist("Counts", dataset_ax, region_ax, dpfcalo_ax)
    items["jetbtag"] = Hist("Counts", dataset_ax, region_ax, btag_ax)
    items["jet_mult"] = Hist("Jets", dataset_ax, region_ax, multiplicity_ax)
    items["bjet_mult"] = Hist("B Jets", dataset_ax, region_ax, multiplicity_ax)
    items["loose_ele_mult"] = Hist("Loose electrons", dataset_ax, region_ax, multiplicity_ax)
    items["tight_ele_mult"] = Hist("Tight electrons", dataset_ax, region_ax, multiplicity_ax)
    items["loose_muo_mult"] = Hist("Loose muons", dataset_ax, region_ax, multiplicity_ax)
    items["tight_muo_mult"] = Hist("Tight muons", dataset_ax, region_ax, multiplicity_ax)
    items["tau_mult"] = Hist("Taus", dataset_ax, region_ax, multiplicity_ax)
    items["photon_mult"] = Hist("Photons", dataset_ax, region_ax, multiplicity_ax)
    items["dphijm"] = Hist("min(4 leading jets, MET)", dataset_ax, region_ax, dphi_ax)

    items["muon_pt"] = Hist("Counts", dataset_ax, region_ax, pt_ax)
    items["muon_eta"] = Hist("Counts", dataset_ax, region_ax, eta_ax)
    items["dimuon_pt"] = Hist("Counts", dataset_ax, region_ax, pt_ax)
    items["dimuon_eta"] = Hist("Counts", dataset_ax, region_ax, eta_ax)
    items["dimuon_mass"] = Hist("Counts", dataset_ax, region_ax, dilepton_mass_ax)

    items["electron_pt"] = Hist("Counts", dataset_ax, region_ax, pt_ax)
    items["electron_eta"] = Hist("Counts", dataset_ax, region_ax, eta_ax)
    items["dielectron_pt"] = Hist("Counts", dataset_ax, region_ax, pt_ax)
    items["dielectron_eta"] = Hist("Counts", dataset_ax, region_ax, eta_ax)
    items["dielectron_mass"] = Hist("Counts", dataset_ax, region_ax, dilepton_mass_ax)

    items["cutflow_sr_j"] = processor.defaultdict_accumulator(int)
    items["cutflow_cr_2m_j"] = processor.defaultdict_accumulator(int)
    items["cutflow_cr_1m_j"] = processor.defaultdict_accumulator(int)
    items["cutflow_cr_2e_j"] = processor.defaultdict_accumulator(int)
    items["cutflow_cr_1e_j"] = processor.defaultdict_accumulator(int)
    items["cutflow_sr_v"] = processor.defaultdict_accumulator(int)
    items["cutflow_cr_2m_v"] = processor.defaultdict_accumulator(int)
    items["cutflow_cr_1m_v"] = processor.defaultdict_accumulator(int)
    items["cutflow_cr_2e_v"] = processor.defaultdict_accumulator(int)
    items["cutflow_cr_1e_v"] = processor.defaultdict_accumulator(int)
    return  processor.dict_accumulator(items)


def setup_candidates(df, cfg):
    muons = JaggedCandidateArray.candidatesfromcounts(
        df['nMuon'],
        pt=df['Muon_pt'],
        eta=df['Muon_eta'],
        phi=df['Muon_phi'],
        mass=df['Muon_mass'],
        charge=df['Muon_charge'],
        mediumId=df['Muon_mediumId'],
        iso=df["Muon_pfRelIso04_all"],
        tightId=df['Muon_tightId']
    )

    electrons = JaggedCandidateArray.candidatesfromcounts(
        df['nElectron'],
        pt=df['Electron_pt'],
        eta=df['Electron_eta'],
        phi=df['Electron_phi'],
        mass=df['Electron_mass'],
        charge=df['Electron_charge'],
        looseId=(df['Electron_cutBased']>=1),
        tightId=(df['Electron_cutBased']==4)
    )
    taus = JaggedCandidateArray.candidatesfromcounts(
        df['nTau'],
        pt=df['Tau_pt'],
        eta=df['Tau_eta'],
        phi=df['Tau_phi'],
        mass=df['Tau_mass'],
        decaymode=df['Tau_idDecayModeNewDMs'],
        clean=df['Tau_cleanmask'],
        iso=df['Tau_idMVAnewDM2017v2'],
        antimu=df['Tau_idAntiMu'],
        antiele=df['Tau_idAntiEle'],
    )

    taus = taus[ (taus.clean==1) \
                         & (taus.decaymode) \
                         & (taus.pt > cfg.TAU.CUTS.PT)\
                         & (np.abs(taus.eta) < cfg.TAU.CUTS.ETA) \
                         & ((taus.iso&2)==2) \
                         & (taus.antimu>0) \
                         & (taus.antiele>0)]

    photons = JaggedCandidateArray.candidatesfromcounts(
        df['nPhoton'],
        pt=df['Photon_pt'],
        eta=df['Photon_eta'],
        phi=df['Photon_phi'],
        mass=df['Photon_mass'],
        id=(df['Photon_cutBased']==1) & (df['Photon_electronVeto']==1),
        clean=df['Photon_cleanmask'],
    )
    photons = photons[(photons.clean==1) \
              & photons.id \
              & (photons.pt > cfg.PHOTON.CUTS.pt) \
              & (np.abs(photons.eta) < cfg.PHOTON.CUTS.eta)]
    jets = JaggedCandidateArray.candidatesfromcounts(
        df['nJet'],
        pt=df['Jet_pt'],
        eta=df['Jet_eta'],
        phi=df['Jet_phi'],
        mass=df['Jet_mass'],

        # Jet ID bit mask:
        # Bit 0 = Loose
        # Bit 1 = Tight
        tightId=(df['Jet_jetId']&2) == 2,
        csvv2=df["Jet_btagCSVV2"],
        deepcsv=df['Jet_btagDeepB'],
        # nef=df['Jet_neEmEF'],
        nhf=df['Jet_neHEF'],
        chf=df['Jet_chHEF'],
        clean=df['Jet_cleanmask']
        # cef=df['Jet_chEmEF'],
    )
    jets = jets[jets.clean==1]
    return jets, muons, electrons, taus, photons

