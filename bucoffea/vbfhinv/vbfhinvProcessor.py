import copy
import coffea.processor as processor
import re
import numpy as np
from dynaconf import settings as cfg

from bucoffea.helpers import (
                              bucoffea_path,
                              dphi,
                              evaluator_from_config,
                              mask_and,
                              mask_or,
                              min_dphi_jet_met,
                              mt,
                              recoil,
                              weight_shape,
                              candidates_in_hem,
                              electrons_in_hem,
                              calculate_vecB,
                              calculate_vecDPhi
                              )
from bucoffea.helpers.dataset import (
                                      extract_year,
                                      is_data,
                                      is_lo_g,
                                      is_lo_w,
                                      is_lo_z,
                                      is_lo_w_ewk,
                                      is_lo_z_ewk,
                                      is_nlo_w,
                                      is_nlo_z,
                                      is_lo_znunu
                                      )
from bucoffea.helpers.gen import (
                                  setup_gen_candidates,
                                  setup_dressed_gen_candidates,
                                  setup_lhe_cleaned_genjets,
                                  fill_gen_v_info
                                 )
from bucoffea.helpers.weights import (
                                  get_veto_weights,
                                  btag_weights,
                                  get_varied_ele_sf,
                                  get_varied_muon_sf
                                 )
from bucoffea.monojet.definitions import (
                                          candidate_weights,
                                          pileup_weights,
                                          setup_candidates,
                                          theory_weights_vbf,
                                          photon_trigger_sf,
                                          photon_impurity_weights,
                                          data_driven_qcd_dataset
                                          )

from bucoffea.vbfhinv.definitions import (
                                           vbfhinv_accumulator,
                                           vbfhinv_regions,
                                           ak4_em_frac_weights,
                                           met_trigger_sf,
                                           apply_hfmask_weights,
                                           apply_hf_weights_for_qcd_estimation,
                                           apply_endcap_weights,
                                           hfmask_sf,
                                           met_xy_correction,
                                           pileup_sf_variations
                                         )

def trigger_selection(selection, df, cfg):
    pass_all = np.zeros(df.size) == 0
    pass_none = ~pass_all
    dataset = df['dataset']

    if df['is_data']:
        selection.add('filt_met', mask_and(df, cfg.FILTERS.DATA))
    else:
        selection.add('filt_met', mask_and(df, cfg.FILTERS.MC))
    selection.add('trig_met', mask_or(df, cfg.TRIGGERS.MET))

    # Electron trigger overlap
    if df['is_data']:
        if "SinglePhoton" in dataset:
            # Backup photon trigger, but not main electron trigger
            trig_ele = mask_or(df, cfg.TRIGGERS.ELECTRON.SINGLE_BACKUP) & (~mask_or(df, cfg.TRIGGERS.ELECTRON.SINGLE))
        elif "SingleElectron" in dataset:
            # Main electron trigger, no check for backup
            trig_ele = mask_or(df, cfg.TRIGGERS.ELECTRON.SINGLE)
        elif "EGamma" in dataset:
            # 2018 has everything in one stream, so simple OR
            trig_ele = mask_or(df, cfg.TRIGGERS.ELECTRON.SINGLE_BACKUP) | mask_or(df, cfg.TRIGGERS.ELECTRON.SINGLE)
        else:
            trig_ele = pass_none
    else:
        trig_ele = mask_or(df, cfg.TRIGGERS.ELECTRON.SINGLE_BACKUP) | mask_or(df, cfg.TRIGGERS.ELECTRON.SINGLE)

    selection.add('trig_ele', trig_ele)

    # Photon trigger:
    if (not df['is_data']) or ('SinglePhoton' in dataset) or ('EGamma' in dataset):
        trig_photon = mask_or(df, cfg.TRIGGERS.PHOTON.SINGLE)
    else:
        trig_photon = pass_none
    selection.add('trig_photon', trig_photon)

    for trgname in cfg.TRIGGERS.HT.GAMMAEFF:
        if (not df['is_data']) or ('JetHT' in dataset):
            selection.add(trgname, mask_or(df,[trgname]))
        else:
            selection.add(trgname, np.ones(df.size)==1)

    # Muon trigger
    selection.add('trig_mu', mask_or(df, cfg.TRIGGERS.MUON.SINGLE))

    return selection

class vbfhinvProcessor(processor.ProcessorABC):
    def __init__(self, blind=False):
        self._year=None
        self._blind=blind
        self._configure()
        self._accumulator = vbfhinv_accumulator(cfg)

    @property
    def accumulator(self):
        return self._accumulator

    def _configure(self, df=None):
        cfg.DYNACONF_WORKS="merge_configs"
        cfg.MERGE_ENABLED_FOR_DYNACONF=True
        cfg.SETTINGS_FILE_FOR_DYNACONF = bucoffea_path("config/vbfhinv.yaml")

        # Reload config based on year
        if df:
            dataset = df['dataset']
            self._year = extract_year(dataset)
            df["year"] = self._year
            cfg.ENV_FOR_DYNACONF = f"era{self._year}"
        else:
            cfg.ENV_FOR_DYNACONF = f"default"
        cfg.reload()

    def process(self, df):
        if not df.size:
            return self.accumulator.identity()
        self._configure(df)
        dataset = df['dataset']
        df['is_lo_w'] = is_lo_w(dataset)
        df['is_lo_z'] = is_lo_z(dataset)
        df['is_lo_znunu'] = is_lo_znunu(dataset)
        df['is_lo_w_ewk'] = is_lo_w_ewk(dataset)
        df['is_lo_z_ewk'] = is_lo_z_ewk(dataset)
        df['is_lo_g'] = is_lo_g(dataset)
        df['is_nlo_z'] = is_nlo_z(dataset)
        df['is_nlo_w'] = is_nlo_w(dataset)
        df['has_lhe_v_pt'] = df['is_lo_w'] | df['is_lo_z'] | df['is_nlo_z'] | df['is_nlo_w'] | df['is_lo_g'] | df['is_lo_w_ewk'] | df['is_lo_z_ewk']
        df['is_data'] = is_data(dataset)

        gen_v_pt = None
        if df['is_lo_w'] or df['is_lo_z'] or df['is_nlo_z'] or df['is_nlo_w'] or df['is_lo_z_ewk'] or df['is_lo_w_ewk']:
            gen = setup_gen_candidates(df)
            dressed = setup_dressed_gen_candidates(df)
            fill_gen_v_info(df, gen, dressed)
            gen_v_pt = df['gen_v_pt_combined']
        elif df['is_lo_g']:
            gen = setup_gen_candidates(df)
            all_gen_photons = gen[(gen.pdg==22)]
            prompt_mask = (all_gen_photons.status==1)&(all_gen_photons.flag&1==1)
            stat1_mask = (all_gen_photons.status==1)
            gen_photons = all_gen_photons[prompt_mask | (~prompt_mask.any()) & stat1_mask ]
            gen_photon = gen_photons[gen_photons.pt.argmax()]

            gen_v_pt = gen_photon.pt.max()

        # Generator-level leading dijet mass
        if df['has_lhe_v_pt']:
            genjets = setup_lhe_cleaned_genjets(df)
            digenjet = genjets[:,:2].distincts()
            df['mjj_gen'] = digenjet.mass.max()
            df['mjj_gen'] = np.where(df['mjj_gen'] > 0, df['mjj_gen'], 0)


        # Candidates
        # Already pre-filtered!
        # All leptons are at least loose
        # Check out setup_candidates for filtering details
        met_pt, met_phi, ak4, bjets, _, muons, electrons, taus, photons = setup_candidates(df, cfg)

        # Remove jets in accordance with the noise recipe
        if not cfg.RUN.ULEGACYV8 and df['year'] == 2017:
            ak4   = ak4[(ak4.ptraw>50) | (ak4.abseta<2.65) | (ak4.abseta>3.139)]
            bjets = bjets[(bjets.ptraw>50) | (bjets.abseta<2.65) | (bjets.abseta>3.139)]

        # Filtering ak4 jets according to pileup ID
        ak4 = ak4[ak4.puid]

        # Recalculate MET pt and phi based on npv-corrections
        if cfg.MET.XYCORR:
            met_pt_uncorr, met_phi_uncorr = met_pt, met_phi
            met_pt, met_phi = met_xy_correction(df, met_pt, met_phi)

        # Muons
        df['is_tight_muon'] = muons.tightId \
                      & (muons.iso < cfg.MUON.CUTS.TIGHT.ISO) \
                      & (muons.pt>cfg.MUON.CUTS.TIGHT.PT) \
                      & (muons.abseta<cfg.MUON.CUTS.TIGHT.ETA)

        dimuons = muons.distincts()
        dimuon_charge = dimuons.i0['charge'] + dimuons.i1['charge']

        df['MT_mu'] = ((muons.counts==1) * mt(muons.pt, muons.phi, met_pt, met_phi)).max()

        # Electrons
        df['is_tight_electron'] = electrons.tightId \
                            & (electrons.pt > cfg.ELECTRON.CUTS.TIGHT.PT) \
                            & (electrons.absetasc < cfg.ELECTRON.CUTS.TIGHT.ETA)

        dielectrons = electrons.distincts()
        dielectron_charge = dielectrons.i0['charge'] + dielectrons.i1['charge']

        df['MT_el'] = ((electrons.counts==1) * mt(electrons.pt, electrons.phi, met_pt, met_phi)).max()

        # ak4
        leadak4_index=ak4.pt.argmax()

        elejet_pairs = ak4[:,:1].cross(electrons)
        df['dREleJet'] = np.hypot(elejet_pairs.i0.eta-elejet_pairs.i1.eta , dphi(elejet_pairs.i0.phi,elejet_pairs.i1.phi)).min()
        muonjet_pairs = ak4[:,:1].cross(muons)
        df['dRMuonJet'] = np.hypot(muonjet_pairs.i0.eta-muonjet_pairs.i1.eta , dphi(muonjet_pairs.i0.phi,muonjet_pairs.i1.phi)).min()

        # Recoil
        df['recoil_pt_uncorr'], df['recoil_phi_uncorr'] = recoil(met_pt_uncorr, met_phi_uncorr, electrons, muons, photons)
        df['recoil_pt'], df['recoil_phi'] = recoil(met_pt,met_phi, electrons, muons, photons)

        df["dPFCaloSR"] = (met_pt - df["CaloMET_pt"]) / met_pt
        df["dPFCaloCR"] = (met_pt - df["CaloMET_pt"]) / df["recoil_pt"]

        df["dPFTkSR"] = (met_pt - df["TkMET_pt"]) / met_pt

        df["minDPhiJetRecoil"] = min_dphi_jet_met(ak4, df['recoil_phi'], njet=4, ptmin=30, etamax=5.0)
        df["minDPhiJetMet"] = min_dphi_jet_met(ak4, met_phi, njet=4, ptmin=30, etamax=5.0)
        selection = processor.PackedSelection()

        # Triggers
        pass_all = np.ones(df.size)==1
        selection.add('inclusive', pass_all)
        selection = trigger_selection(selection, df, cfg)

        selection.add('mu_pt_trig_safe', muons.pt.max() > 30)

        # Common selection
        selection.add('veto_ele', electrons.counts==0)
        selection.add('veto_muo', muons.counts==0)
        selection.add('veto_photon', photons.counts==0)
        selection.add('veto_tau', taus.counts==0)
        selection.add('at_least_one_tau', taus.counts>0)
        selection.add('mindphijr',df['minDPhiJetRecoil'] > cfg.SELECTION.SIGNAL.MINDPHIJR)
        selection.add('mindphijm',df['minDPhiJetMet'] > cfg.SELECTION.SIGNAL.MINDPHIJR)

        # Inverted min DPhi(j,met) cut for QCD CR
        selection.add('mindphijr_inv', df['minDPhiJetRecoil'] <= cfg.SELECTION.SIGNAL.MINDPHIJR)

        # B jets are treated using veto weights
        # So accept them in MC, but reject in data
        if df['is_data']:
            selection.add('veto_b', bjets.counts==0)
        else:
            selection.add('veto_b', pass_all)

        selection.add('dpfcalo_sr',np.abs(df['dPFCaloSR']) < cfg.SELECTION.SIGNAL.DPFCALO)
        selection.add('dpfcalo_cr',np.abs(df['dPFCaloCR']) < cfg.SELECTION.SIGNAL.DPFCALO)

        selection.add('recoil', df['recoil_pt']>cfg.SELECTION.SIGNAL.RECOIL)
        selection.add('met_sr', met_pt>cfg.SELECTION.SIGNAL.RECOIL)
        selection.add('small_met', met_pt<50)

        # Relaxed recoil cut for Zmm region
        selection.add('recoil_zmm', df['recoil_pt']>100)

        # AK4 dijet
        diak4 = ak4[:,:2].distincts()
        leadak4_pt_eta = (diak4.i0.pt > cfg.SELECTION.SIGNAL.LEADAK4.PT) & (np.abs(diak4.i0.eta) < cfg.SELECTION.SIGNAL.LEADAK4.ETA)
        trailak4_pt_eta = (diak4.i1.pt > cfg.SELECTION.SIGNAL.TRAILAK4.PT) & (np.abs(diak4.i1.eta) < cfg.SELECTION.SIGNAL.TRAILAK4.ETA)
        hemisphere = (diak4.i0.eta * diak4.i1.eta < 0).any()
        has_track0 = np.abs(diak4.i0.eta) <= 2.5
        has_track1 = np.abs(diak4.i1.eta) <= 2.5

        leadak4_id = diak4.i0.tightId & (has_track0*((diak4.i0.chf > cfg.SELECTION.SIGNAL.LEADAK4.CHF) & (diak4.i0.nhf < cfg.SELECTION.SIGNAL.LEADAK4.NHF)) + ~has_track0)
        trailak4_id = has_track1*((diak4.i1.chf > cfg.SELECTION.SIGNAL.TRAILAK4.CHF) & (diak4.i1.nhf < cfg.SELECTION.SIGNAL.TRAILAK4.NHF)) + ~has_track1

        def get_more_central_jeteta(diak4):
            mask = diak4.i0.abseta > diak4.i1.abseta
            return (mask * diak4.i1.eta) + (~mask * diak4.i0.eta)

        def get_more_forward_jeteta(diak4):
            mask = diak4.i0.abseta > diak4.i1.abseta
            return (mask * diak4.i0.eta) + (~mask * diak4.i1.eta)

        df['mjj'] = diak4.mass.max()
        df['dphijj'] = dphi(diak4.i0.phi.min(), diak4.i1.phi.max())
        df['detajj'] = np.abs(diak4.i0.eta - diak4.i1.eta).max()

        df['ak4_mt0'] = mt(diak4.i0.pt, diak4.i0.phi, met_pt, met_phi).max()
        df['ak4_mt1'] = mt(diak4.i1.pt, diak4.i1.phi, met_pt, met_phi).max()

        df['dphi_ak40_met'] = dphi(diak4.i0.phi.min(), met_phi)
        df['dphi_ak41_met'] = dphi(diak4.i1.phi.min(), met_phi)

        leading_jet_in_horn = ((diak4.i0.abseta<3.2) & (diak4.i0.abseta>2.8)).any()
        trailing_jet_in_horn = ((diak4.i1.abseta<3.2) & (diak4.i1.abseta>2.8)).any()

        # Additional jets in the central region (|eta| < 2.5)
        extra_ak4 = ak4[:,2:]
        extra_ak4_central = extra_ak4[extra_ak4.abseta < 2.5]

        # selection.add('hornveto', (df['dPFTkSR'] < 0.8) | ~(leading_jet_in_horn | trailing_jet_in_horn))
        
        df['htmiss'] = ak4[ak4.pt>30].p4.sum().pt
        df['ht'] = ak4[ak4.pt>30].pt.sum()

        if df['year'] == 2018:
            if df['is_data']:
                metphihem_mask = ~((met_phi > -1.8) & (met_phi < -0.6) & (df['run'] > 319077))
            else:
                metphihem_mask = pass_all
            selection.add("metphihemextveto", metphihem_mask)
            selection.add('no_el_in_hem', electrons[electrons_in_hem(electrons)].counts==0)
        else:
            selection.add("metphihemextveto", pass_all)
            selection.add('no_el_in_hem', pass_all)

        # Sigma eta & phi cut (only for v8 samples because we have the info there)
        if cfg.RUN.ULEGACYV8:
            jets_for_cut = ak4[(ak4.pt > cfg.RUN.HF_PT_THRESH) & (ak4.abseta > 2.99) & (ak4.abseta < 5.0)]

            # We will only consider jets that are back to back with MET i.e. dPhi(jet,MET) > 2.5
            dphi_hfjet_met = dphi(jets_for_cut.phi, df['recoil_phi'])
            dphimask = dphi_hfjet_met > 2.5
            jets_for_cut = jets_for_cut[dphimask]

            seta_minus_phi_alljets = jets_for_cut.setaeta - jets_for_cut.sphiphi

            # Cut away the low sigma eta & phi corner (< 0.02)
            setaphi_corner_cut = ~((jets_for_cut.setaeta < 0.02) & (jets_for_cut.sphiphi < 0.02))
            # Sigma eta - phi < 0.02 requirement
            setaphi_diff_cut_alljets = (seta_minus_phi_alljets < 0.02)

            # For jets with |eta| > 4, we have a different requirement
            setaphi_cut_higheta = (jets_for_cut.setaeta < 0.1) & (jets_for_cut.sphiphi > 0.02)

            is_high_eta_jet = jets_for_cut.abseta > 4.0
            setaphi_cut_alljets = (is_high_eta_jet * setaphi_cut_higheta + ~is_high_eta_jet * (setaphi_corner_cut & setaphi_diff_cut_alljets)).all()
            
            stripsize_cut_alljets = (jets_for_cut.hfcentralstripsize < 3).all()

            fail_hf_cuts = (~setaphi_cut_alljets) | (~stripsize_cut_alljets)
            
            selection.add('sigma_eta_minus_phi', setaphi_cut_alljets)
            selection.add('central_stripsize_cut', stripsize_cut_alljets)
            selection.add('fail_hf_cuts', fail_hf_cuts)

        else:
            selection.add('sigma_eta_minus_phi', pass_all)
            selection.add('central_stripsize_cut', pass_all)
            selection.add('fail_hf_cuts', pass_all)

        selection.add('two_jets', diak4.counts>0)
        selection.add('leadak4_pt_eta', leadak4_pt_eta.any())
        selection.add('trailak4_pt_eta', trailak4_pt_eta.any())
        selection.add('hemisphere', hemisphere)
        selection.add('leadak4_id',leadak4_id.any())
        selection.add('trailak4_id',trailak4_id.any())
        selection.add('mjj', df['mjj'] > cfg.SELECTION.SIGNAL.DIJET.SHAPE_BASED.MASS)
        selection.add('dphijj', df['dphijj'] < cfg.SELECTION.SIGNAL.DIJET.SHAPE_BASED.DPHI)
        selection.add('detajj', df['detajj'] > cfg.SELECTION.SIGNAL.DIJET.SHAPE_BASED.DETA)

        # Cleaning cuts for signal region

        # NEF cut: Only for endcap jets, require NEF < 0.7
        ak40_in_endcap = (diak4.i0.abseta > 2.5) & (diak4.i0.abseta < 3.0)
        ak41_in_endcap = (diak4.i1.abseta > 2.5) & (diak4.i1.abseta < 3.0)

        max_neEmEF_ak40 = (~ak40_in_endcap) | (diak4.i0.nef < 0.7) 
        max_neEmEF_ak41 = (~ak41_in_endcap) | (diak4.i1.nef < 0.7) 

        max_neEmEF = (max_neEmEF_ak40 & max_neEmEF_ak41).any()
        # selection.add('max_neEmEF', max_neEmEF)
        
        vec_b = calculate_vecB(ak4, met_pt, met_phi)
        vec_dphi = calculate_vecDPhi(ak4, met_pt, met_phi, df['TkMET_phi'])

        dphitkpf = dphi(met_phi, df['TkMET_phi'])

        no_jet_in_trk = (diak4.i0.abseta>2.5).any() & (diak4.i1.abseta>2.5).any()
        no_jet_in_hf = (diak4.i0.abseta<3.0).any() & (diak4.i1.abseta<3.0).any()

        at_least_one_jet_in_hf = (diak4.i0.abseta>3.0).any() | (diak4.i1.abseta>3.0).any()
        at_least_one_jet_in_trk = (diak4.i0.abseta<2.5).any() | (diak4.i1.abseta<2.5).any()

        # Categorized cleaning cuts
        eemitigation = (
                        (no_jet_in_hf | at_least_one_jet_in_trk) & (vec_dphi < 1.0)
                    ) | (
                        (no_jet_in_trk & at_least_one_jet_in_hf) & (vec_b < 0.2)
                    )

        selection.add('eemitigation', eemitigation)

        # HF-HF veto in SR
        both_jets_in_hf = (diak4.i0.abseta > 3.0) & (diak4.i1.abseta > 3.0)
        selection.add('veto_hfhf', ~both_jets_in_hf.any())

        # Leading jet |eta| < 2.9
        # leadak4_not_in_hf = (diak4.i0.abseta < 2.9).any()
        # selection.add('leadak4_not_in_hf', leadak4_not_in_hf)

        # Reject events where the leading jet has momentum > 6.5 TeV
        leadak4_clean = diak4.i0.pt * np.cosh(diak4.i0.eta) < 6500
        selection.add('leadak4_clean', leadak4_clean.any())

        # Divide into three categories for trigger study
        if cfg.RUN.TRIGGER_STUDY:
            two_central_jets = (np.abs(diak4.i0.eta) <= 2.5) & (np.abs(diak4.i1.eta) <= 2.5)
            two_hf_jets = (np.abs(diak4.i0.eta) > 3.0) & (np.abs(diak4.i1.eta) > 3.0)
            one_jet_forward_one_jet_central = (~two_central_jets) & (~two_hf_jets)

            selection.add('two_central_jets', two_central_jets.any())
            selection.add('one_jet_forward_one_jet_central', one_jet_forward_one_jet_central.any())
        
        # Mask for 1/5th unlbinding
        one_fifth_mask = ~pass_all

        # Only pick each 5 entry in data
        one_fifth_mask[::5] = True
        if df['is_data']:
            selection.add('one_fifth_mask', one_fifth_mask)
        else:
            selection.add('one_fifth_mask', pass_all)

        # Dimuon CR
        leadmuon_index=muons.pt.argmax()
        selection.add('at_least_one_tight_mu', df['is_tight_muon'].any())
        selection.add('dimuon_mass', ((dimuons.mass > cfg.SELECTION.CONTROL.DOUBLEMU.MASS.MIN) \
                                    & (dimuons.mass < cfg.SELECTION.CONTROL.DOUBLEMU.MASS.MAX)).any())
        selection.add('dimuon_charge', (dimuon_charge==0).any())
        selection.add('two_muons', muons.counts==2)

        # Single muon CR
        selection.add('one_muon', muons.counts==1)
        selection.add('mt_mu', df['MT_mu'] < cfg.SELECTION.CONTROL.SINGLEMU.MT)

        # Diele CR
        leadelectron_index=electrons.pt.argmax()

        selection.add('one_electron', electrons.counts==1)
        selection.add('two_electrons', electrons.counts==2)
        selection.add('at_least_one_tight_el', df['is_tight_electron'].any())


        selection.add('dielectron_mass', ((dielectrons.mass > cfg.SELECTION.CONTROL.DOUBLEEL.MASS.MIN)  \
                                        & (dielectrons.mass < cfg.SELECTION.CONTROL.DOUBLEEL.MASS.MAX)).any())
        selection.add('dielectron_charge', (dielectron_charge==0).any())

        # Single Ele CR
        selection.add('met_el', met_pt > cfg.SELECTION.CONTROL.SINGLEEL.MET)
        selection.add('mt_el', df['MT_el'] < cfg.SELECTION.CONTROL.SINGLEEL.MT)

        # Photon CR
        leadphoton_index=photons.pt.argmax()

        df['is_tight_photon'] = photons.mediumId & photons.barrel

        selection.add('one_photon', photons.counts==1)
        selection.add('at_least_one_tight_photon', df['is_tight_photon'].any())
        selection.add('photon_pt', photons.pt.max() > cfg.PHOTON.CUTS.TIGHT.PT)
        selection.add('photon_pt_trig', photons.pt.max() > cfg.PHOTON.CUTS.TIGHT.PTTRIG)

        # Fill histograms
        output = self.accumulator.identity()

        # Gen
        if df['has_lhe_v_pt']:
            output['genvpt_check'].fill(vpt=gen_v_pt,type="Nano", dataset=dataset)

        if 'LHE_Njets' in df:
            output['lhe_njets'].fill(dataset=dataset, multiplicity=df['LHE_Njets'])
        if 'LHE_HT' in df:
            output['lhe_ht'].fill(dataset=dataset, ht=df['LHE_HT'])
        if 'LHE_HTIncoming' in df:
            output['lhe_htinc'].fill(dataset=dataset, ht=df['LHE_HTIncoming'])

        # Weights
        evaluator = evaluator_from_config(cfg)

        weights = processor.Weights(size=df.size, storeIndividual=True)
        if not df['is_data']:
            weights.add('gen', df['Generator_weight'])

            try:
                weights.add('prefire', df['PrefireWeight'])
            except KeyError:
                weights.add('prefire', np.ones(df.size))

            weights = candidate_weights(weights, df, evaluator, muons, electrons, photons, cfg)

            # EWK corrections to VBF signal
            if cfg.RUN.APPLY_EWK_CORR_TO_SIGNAL:
                if re.match('VBF_HToInv.*', df['dataset']):
                    # Get Higgs pt from GEN collection
                    gen = setup_gen_candidates(df)
                    higgs_pt = gen[(gen.pdg==25)&(gen.status==62)].pt.max()

                    def ewk_correction(a, b):
                        return (1 + a * higgs_pt + b) / 0.95
                     
                    coeff = [-0.000372, -0.0304]
                    ewk_corr_signal = ewk_correction(*coeff)
                    weights.add('ewk_corr_signal', ewk_corr_signal)

            # B jet veto weights
            bsf_variations = btag_weights(bjets,cfg)
            weights.add("bveto", (1-bsf_variations["central"]).prod())

            weights = pileup_weights(weights, df, evaluator, cfg)
            if cfg.RUN.APPLY_CLEANING_CUTS:
                weights = ak4_em_frac_weights(weights, diak4, evaluator)
            if not (gen_v_pt is None):
                weights = theory_weights_vbf(weights, df, evaluator, gen_v_pt, df['mjj_gen'])

        # Save per-event values for synchronization
        if cfg.RUN.KINEMATICS.SAVE:
            for event in cfg.RUN.KINEMATICS.EVENTS:
                mask = df['event'] == event
                if not mask.any():
                    continue
                output['kinematics']['event'] += [event]
                output['kinematics']['met'] += [met_pt[mask]]
                output['kinematics']['met_phi'] += [met_phi[mask]]
                output['kinematics']['recoil'] += [df['recoil_pt'][mask]]
                output['kinematics']['recoil_phi'] += [df['recoil_phi'][mask]]

                output['kinematics']['ak4pt0'] += [ak4[leadak4_index][mask].pt]
                output['kinematics']['ak4eta0'] += [ak4[leadak4_index][mask].eta]
                output['kinematics']['leadbtag'] += [ak4.pt.max()<0][mask]

                output['kinematics']['nLooseMu'] += [muons.counts[mask]]
                output['kinematics']['nTightMu'] += [muons[df['is_tight_muon']].counts[mask]]
                output['kinematics']['mupt0'] += [muons[leadmuon_index][mask].pt]
                output['kinematics']['mueta0'] += [muons[leadmuon_index][mask].eta]

                output['kinematics']['nLooseEl'] += [electrons.counts[mask]]
                output['kinematics']['nTightEl'] += [electrons[df['is_tight_electron']].counts[mask]]
                output['kinematics']['elpt0'] += [electrons[leadelectron_index][mask].pt]
                output['kinematics']['eleta0'] += [electrons[leadelectron_index][mask].eta]

                output['kinematics']['nLooseGam'] += [photons.counts[mask]]
                output['kinematics']['nTightGam'] += [photons[df['is_tight_photon']].counts[mask]]
                output['kinematics']['gpt0'] += [photons[leadphoton_index][mask].pt]
                output['kinematics']['geta0'] += [photons[leadphoton_index][mask].eta]


        # Sum of all weights to use for normalization
        # TODO: Deal with systematic variations
        output['nevents'][dataset] += df.size
        if not df['is_data']:
            output['sumw'][dataset] +=  df['genEventSumw']
            output['sumw2'][dataset] +=  df['genEventSumw2']
            output['sumw_pileup'][dataset] +=  weights._weights['pileup'].sum()

        regions = vbfhinv_regions(cfg)

        # Get veto weights (only for MC)
        if not df['is_data']:
            veto_weights = get_veto_weights(df, cfg, evaluator, electrons, muons, taus, do_variations=cfg.RUN.VETO_WEIGHTS_STUDY)
        
        for region, cuts in regions.items():
            if not re.match(cfg.RUN.REGIONREGEX, region):
                continue
            # Run on selected regions only
            exclude = [None]
            if region == 'sr_vbf_no_pu':
                exclude = ['pileup']
            region_weights = copy.deepcopy(weights)

            if not df['is_data']:
                ### Trigger weights
                if re.match(r'cr_(\d+)e.*', region):
                    p_pass_data = 1 - (1-evaluator["trigger_electron_eff_data"](electrons.etasc, electrons.pt)).prod()
                    p_pass_mc   = 1 - (1-evaluator["trigger_electron_eff_mc"](electrons.etasc, electrons.pt)).prod()
                    trigger_weight = p_pass_data/p_pass_mc
                    trigger_weight[np.isnan(trigger_weight)] = 1
                    region_weights.add('trigger', trigger_weight)
                elif re.match(r'cr_(\d+)m.*', region) or re.match('sr_.*', region):
                    met_trigger_sf(region_weights, diak4, df, apply_categorized=cfg.RUN.APPLY_CATEGORIZED_SF)
                elif re.match(r'cr_g.*', region):
                    photon_trigger_sf(region_weights, photons, df)

                if cfg.RUN.APPLY_HF_CUTS:
                    region_weights = hfmask_sf(ak4, region_weights, evaluator, df, cfg)
                if cfg.RUN.APPLY_WEIGHTS.HFMASK:
                    region_weights = apply_hfmask_weights(ak4, region_weights, evaluator, met_phi, cfg)
                if cfg.RUN.APPLY_WEIGHTS.ENDCAP:
                    region_weights = apply_endcap_weights(diak4, region_weights, evaluator)

                # Veto weights
                if re.match('.*no_veto.*', region):
                    exclude = [
                            "muon_id_iso_tight",
                            "muon_id_tight",
                            "muon_iso_tight",
                            "muon_id_loose",
                            "muon_iso_loose",
                            "ele_reco",
                            "ele_id_tight",
                            "ele_id_loose",
                            "tau_id"
                        ]
                    region_weights.add("veto",veto_weights.partial_weight(include=["nominal"]))

                # SR without prefiring weights applied
                if region == 'sr_vbf_no_pref':
                    exclude = ['prefire']

                # HEM-veto weights for signal region MC
                if re.match('^sr_vbf.*', region) and df['year'] == 2018 and 'no_hem_veto' not in region:
                    # Events that lie in the HEM-veto region
                    events_to_weight_mask = (met_phi > -1.8) & (met_phi < -0.6)
                    # Weight is the "good lumi fraction" for 2018
                    weight = 21.1/59.7
                    hem_weight = np.where(events_to_weight_mask, weight, 1.0)

                    region_weights.add("hem_weight", hem_weight)

            # Weights for QCD estimation
            if cfg.RUN.QCD_ESTIMATION:
                apply_hf_weights_for_qcd_estimation(ak4, region_weights, evaluator, df, cfg, region) 

            # This is the default weight for this region
            rweight = region_weights.partial_weight(exclude=exclude)

            # Blinding
            if(self._blind and df['is_data'] and region.startswith('sr')):
                continue

            # Cutflow plot for signal and control regions
            if any(x in region for x in ["sr", "cr", "tr"]):
                output['cutflow_' + region][dataset]['all']+=df.size
                for icut, cutname in enumerate(cuts):
                    output['cutflow_' + region][dataset][cutname] += selection.all(*cuts[:icut+1]).sum()

            mask = selection.all(*cuts)

            if cfg.RUN.SAVE.TREE:
                if region in ['sr_vbf', 'sr_vbf_no_veto_all']:
                    output['tree_int64'][region]["event"]             +=  processor.column_accumulator(df["event"][mask])
                    output['tree_int64'][region]["run"]               +=  processor.column_accumulator(df["run"][mask])
                    output['tree_int64'][region]["lumi"]              +=  processor.column_accumulator(df["luminosityBlock"][mask])
                    
                    output['tree_float16'][region]["leadak4_pt"]        +=  processor.column_accumulator(np.float16(diak4.i0.pt[mask]))
                    output['tree_float16'][region]["leadak4_eta"]       +=  processor.column_accumulator(np.float16(diak4.i0.eta[mask]))
                    output['tree_float16'][region]["leadak4_phi"]       +=  processor.column_accumulator(np.float16(diak4.i0.phi[mask]))
                    output['tree_float16'][region]["leadak4_nef"]       +=  processor.column_accumulator(np.float16(diak4.i0.nef[mask]))
                    output['tree_float16'][region]["leadak4_nhf"]       +=  processor.column_accumulator(np.float16(diak4.i0.nhf[mask]))
                    output['tree_float16'][region]["leadak4_chf"]       +=  processor.column_accumulator(np.float16(diak4.i0.chf[mask]))
                    output['tree_float16'][region]["leadak4_cef"]       +=  processor.column_accumulator(np.float16(diak4.i0.cef[mask]))

                    if cfg.RUN.ULEGACYV8:
                        output['tree_float16'][region]["leadak4_setaeta"]   +=  processor.column_accumulator(np.float16(diak4.i0.setaeta[mask]))
                        output['tree_float16'][region]["leadak4_sphiphi"]   +=  processor.column_accumulator(np.float16(diak4.i0.sphiphi[mask]))
                        output['tree_float16'][region]["leadak4_cssize"]    +=  processor.column_accumulator(np.float16(diak4.i0.hfcentralstripsize[mask]))
                        output['tree_float16'][region]["leadak4_btagDeepFlavQG"]    +=  processor.column_accumulator(np.float16(diak4.i0.btagdf[mask]))
                
                    output['tree_float16'][region]["trailak4_pt"]        +=  processor.column_accumulator(np.float16(diak4.i1.pt[mask]))
                    output['tree_float16'][region]["trailak4_eta"]       +=  processor.column_accumulator(np.float16(diak4.i1.eta[mask]))
                    output['tree_float16'][region]["trailak4_phi"]       +=  processor.column_accumulator(np.float16(diak4.i1.phi[mask]))
                    output['tree_float16'][region]["trailak4_nef"]       +=  processor.column_accumulator(np.float16(diak4.i1.nef[mask]))
                    output['tree_float16'][region]["trailak4_nhf"]       +=  processor.column_accumulator(np.float16(diak4.i1.nhf[mask]))
                    output['tree_float16'][region]["trailak4_chf"]       +=  processor.column_accumulator(np.float16(diak4.i1.chf[mask]))
                    output['tree_float16'][region]["trailak4_cef"]       +=  processor.column_accumulator(np.float16(diak4.i1.cef[mask]))

                    if cfg.RUN.ULEGACYV8:
                        output['tree_float16'][region]["trailak4_setaeta"]   +=  processor.column_accumulator(np.float16(diak4.i1.setaeta[mask]))
                        output['tree_float16'][region]["trailak4_sphiphi"]   +=  processor.column_accumulator(np.float16(diak4.i1.sphiphi[mask]))
                        output['tree_float16'][region]["trailak4_cssize"]    +=  processor.column_accumulator(np.float16(diak4.i1.hfcentralstripsize[mask]))
                        output['tree_float16'][region]["trailak4_btagDeepFlavQG"]    +=  processor.column_accumulator(np.float16(diak4.i1.btagdf[mask]))

                    output['tree_float16'][region]["mjj"]               +=  processor.column_accumulator(np.float16(df["mjj"][mask]))
                    output['tree_float16'][region]["detajj"]            +=  processor.column_accumulator(np.float16(df["detajj"][mask]))
                    output['tree_float16'][region]["dphijj"]            +=  processor.column_accumulator(np.float16(df["dphijj"][mask]))
                    output['tree_float16'][region]["uncorr_recoil_pt"]  +=  processor.column_accumulator(np.float16(df["recoil_pt_uncorr"][mask]))
                    output['tree_float16'][region]["uncorr_recoil_phi"] +=  processor.column_accumulator(np.float16(df["recoil_phi_uncorr"][mask]))
                    output['tree_float16'][region]["recoil_pt"]         +=  processor.column_accumulator(np.float16(df["recoil_pt"][mask]))
                    output['tree_float16'][region]["recoil_phi"]        +=  processor.column_accumulator(np.float16(df["recoil_phi"][mask]))
                    output['tree_float16'][region]["uncorr_met_pt"]     +=  processor.column_accumulator(np.float16(met_pt_uncorr[mask]))
                    output['tree_float16'][region]["uncorr_met_phi"]    +=  processor.column_accumulator(np.float16(met_phi_uncorr[mask]))
                    output['tree_float16'][region]["met_pt"]            +=  processor.column_accumulator(np.float16(met_pt[mask]))
                    output['tree_float16'][region]["met_phi"]           +=  processor.column_accumulator(np.float16(met_phi[mask]))
                    output['tree_float16'][region]["CaloMet_pt"]        +=  processor.column_accumulator(np.float16(df['CaloMET_pt'][mask]))
                    output['tree_float16'][region]["CaloMet_phi"]       +=  processor.column_accumulator(np.float16(df['CaloMET_phi'][mask]))
                    output['tree_float16'][region]["minDPhiJetMet"]     +=  processor.column_accumulator(np.float16(df["minDPhiJetMet"][mask]))
                    output['tree_float16'][region]["minDPhiJetRecoil"]  +=  processor.column_accumulator(np.float16(df["minDPhiJetRecoil"][mask]))
                    output['tree_float16'][region]["dphi_ak40_met"]     +=  processor.column_accumulator(np.float16(df["dphi_ak40_met"][mask]))
                    output['tree_float16'][region]["dphi_ak41_met"]     +=  processor.column_accumulator(np.float16(df["dphi_ak41_met"][mask]))

                    output['tree_float16'][region]["htmiss"]            +=  processor.column_accumulator(np.float16(df['htmiss'][mask]))
                    output['tree_float16'][region]["ht"]                +=  processor.column_accumulator(np.float16(df['ht'][mask]))
                    output['tree_float16'][region]["vecb"]              +=  processor.column_accumulator(np.float16(vec_b[mask]))
                    output['tree_float16'][region]["vecdphi"]           +=  processor.column_accumulator(np.float16(vec_dphi[mask]))
                    output['tree_float16'][region]["dphitkpf"]          +=  processor.column_accumulator(np.float16(dphitkpf[mask]))
                    
                    output['tree_float16'][region]["nLooseMuon"]        +=  processor.column_accumulator(np.float16(muons.counts[mask]))
                    output['tree_float16'][region]["nLooseElectron"]    +=  processor.column_accumulator(np.float16(electrons.counts[mask]))
                    output['tree_float16'][region]["nLooseTau"]         +=  processor.column_accumulator(np.float16(taus.counts[mask]))
                    
                    event_has_ele = electrons[mask].counts != 0
                    ele_pt = np.where(event_has_ele, electrons.pt.max()[mask], -999)
                    ele_eta = np.where(event_has_ele, electrons[electrons.pt.argmax()].eta.max()[mask], -999)
                    ele_phi = np.where(event_has_ele, electrons[electrons.pt.argmax()].phi.max()[mask], -999)

                    event_has_mu = muons[mask].counts != 0
                    mu_pt = np.where(event_has_mu, muons.pt.max()[mask], -999)
                    mu_eta = np.where(event_has_mu, muons[muons.pt.argmax()].eta.max()[mask], -999)
                    mu_phi = np.where(event_has_mu, muons[muons.pt.argmax()].phi.max()[mask], -999)

                    output['tree_float16'][region]["lead_ele_pt"]   +=  processor.column_accumulator(ele_pt)
                    output['tree_float16'][region]["lead_ele_eta"]   +=  processor.column_accumulator(ele_eta)
                    output['tree_float16'][region]["lead_ele_phi"]   +=  processor.column_accumulator(ele_phi)

                    output['tree_float16'][region]["lead_muon_pt"]   +=  processor.column_accumulator(mu_pt)
                    output['tree_float16'][region]["lead_muon_eta"]   +=  processor.column_accumulator(mu_eta)
                    output['tree_float16'][region]["lead_muon_phi"]   +=  processor.column_accumulator(mu_phi)

                    if not df['is_data']:                    
                        ele_genpartflav = np.where(event_has_ele, electrons[electrons.pt.argmax()].genpartflav.max()[mask], -999)
                        mu_genpartflav = np.where(event_has_mu, muons[muons.pt.argmax()].genpartflav.max()[mask], -999)
                        output['tree_float16'][region]["lead_ele_genpartflav"]   +=  processor.column_accumulator(ele_genpartflav)
                        output['tree_float16'][region]["lead_muon_genpartflav"]   +=  processor.column_accumulator(mu_genpartflav)

                    for name, w in region_weights._weights.items():
                        output['tree_float16'][region][f"weight_{name}"] += processor.column_accumulator(np.float16(w[mask]))
                    
                    output['tree_float16'][region][f"weight_total"] += processor.column_accumulator(np.float16(rweight[mask]))

            # Save the event numbers of events passing this selection
            if cfg.RUN.SAVE.PASSING:
                output['selected_events'][region] += list(df['event'][mask])


            # Multiplicities
            def fill_mult(name, candidates):
                output[name].fill(
                                  dataset=dataset,
                                  region=region,
                                  multiplicity=candidates[mask].counts,
                                  weight=rweight[mask]
                                  )

            fill_mult('ak4_mult', ak4[ak4.pt>30])
            fill_mult('bjet_mult',bjets)
            fill_mult('loose_ele_mult',electrons)
            fill_mult('tight_ele_mult',electrons[df['is_tight_electron']])
            fill_mult('loose_muo_mult',muons)
            fill_mult('tight_muo_mult',muons[df['is_tight_muon']])
            fill_mult('tau_mult',taus)
            fill_mult('photon_mult',photons)
            # Number of additional jets in the central region, |eta| < 2.5
            fill_mult('extra_ak4_mult', extra_ak4_central)

            def ezfill(name, **kwargs):
                """Helper function to make filling easier."""
                output[name].fill(
                                  dataset=dataset,
                                  region=region,
                                  **kwargs
                                  )

            # Monitor weights
            for wname, wvalue in region_weights._weights.items():
                ezfill("weights", weight_type=wname, weight_value=wvalue[mask])
                ezfill("weights_wide", weight_type=wname, weight_value=wvalue[mask])

            # All ak4
            # This is a workaround to create a weight array of the right dimension
            w_alljets = weight_shape(ak4[mask].eta, rweight[mask])
            w_alljets_nopref = weight_shape(ak4[mask].eta, region_weights.partial_weight(exclude=exclude+['prefire'])[mask])

            ezfill('ak4_eta',    jeteta=ak4[mask].eta.flatten(), weight=w_alljets)
            ezfill('ak4_phi',    jetphi=ak4[mask].phi.flatten(), weight=w_alljets)
            ezfill('ak4_pt',     jetpt=ak4[mask].pt.flatten(),   weight=w_alljets)

            ezfill('ak4_eta_phi',    jeteta=ak4[mask].eta.flatten(),   jetphi=ak4[mask].phi.flatten(),   weight=w_alljets)

            ezfill('ak4_eta_nopref',    jeteta=ak4[mask].eta.flatten(), weight=w_alljets_nopref)
            ezfill('ak4_phi_nopref',    jetphi=ak4[mask].phi.flatten(), weight=w_alljets_nopref)
            ezfill('ak4_pt_nopref',     jetpt=ak4[mask].pt.flatten(),   weight=w_alljets_nopref)

            # Leading ak4
            w_diak4 = weight_shape(diak4.pt[mask], rweight[mask])
            ezfill('ak4_eta0',      jeteta=diak4.i0.eta[mask].flatten(),    weight=w_diak4)
            ezfill('ak4_phi0',      jetphi=diak4.i0.phi[mask].flatten(),    weight=w_diak4)
            ezfill('ak4_pt0',       jetpt=diak4.i0.pt[mask].flatten(),      weight=w_diak4)
            ezfill('ak4_ptraw0',    jetpt=diak4.i0.ptraw[mask].flatten(),   weight=w_diak4)
            ezfill('ak4_chf0',      frac=diak4.i0.chf[mask].flatten(),      weight=w_diak4)
            ezfill('ak4_nhf0',      frac=diak4.i0.nhf[mask].flatten(),      weight=w_diak4)
            ezfill('ak4_nef0',      frac=diak4.i0.nef[mask].flatten(),      weight=w_diak4)
            ezfill('ak4_nconst0',   nconst=diak4.i0.nconst[mask].flatten(), weight=w_diak4)

            # ezfill('ak4_pt0_eta0',  jetpt=diak4.i0.pt[mask].flatten(),     jeteta=diak4.i0.eta[mask].flatten(),     weight=w_diak4)

            # Trailing ak4
            ezfill('ak4_eta1',      jeteta=diak4.i1.eta[mask].flatten(),    weight=w_diak4)
            ezfill('ak4_phi1',      jetphi=diak4.i1.phi[mask].flatten(),    weight=w_diak4)
            ezfill('ak4_pt1',       jetpt=diak4.i1.pt[mask].flatten(),      weight=w_diak4)
            ezfill('ak4_ptraw1',    jetpt=diak4.i1.ptraw[mask].flatten(),   weight=w_diak4)
            ezfill('ak4_chf1',      frac=diak4.i1.chf[mask].flatten(),      weight=w_diak4)
            ezfill('ak4_nhf1',      frac=diak4.i1.nhf[mask].flatten(),      weight=w_diak4)
            ezfill('ak4_nef1',      frac=diak4.i1.nef[mask].flatten(),      weight=w_diak4)
            ezfill('ak4_nconst1',   nconst=diak4.i1.nconst[mask].flatten(), weight=w_diak4)

            # Eta of more central and more forward VBF jets
            ezfill('ak4_central_eta',    jeteta=get_more_central_jeteta(diak4)[mask].flatten(),    weight=w_diak4)
            ezfill('ak4_forward_eta',    jeteta=get_more_forward_jeteta(diak4)[mask].flatten(),    weight=w_diak4)

            if cfg.RUN.ULEGACYV8:
                def is_hf_jet(_ak4, ptmin=100, etamin=2.9, etamax=5.0):
                    return (_ak4.pt > ptmin) & (_ak4.abseta > etamin) & (_ak4.abseta < etamax)

                hfmask = is_hf_jet(ak4[mask])

                # Consider only high pt jets with pt > 100 GeV
                w_hfjets = np.where(
                    hfmask.flatten(),
                    w_alljets,
                    0.
                )

                ezfill('ak4_sigma_eta_eta',   sigmaetaeta=ak4[mask].setaeta.flatten(),        jeta=ak4[mask].abseta.flatten(),   weight=w_hfjets)
                ezfill('ak4_sigma_phi_phi',   sigmaphiphi=ak4[mask].sphiphi.flatten(),        jeta=ak4[mask].abseta.flatten(),   weight=w_hfjets)

                # 2D sigma eta vs. phi
                ezfill('ak4_sigma_eta_phi',   
                        sigmaetaeta=ak4[mask].setaeta.flatten(),    
                        sigmaphiphi=ak4[mask].sphiphi.flatten(),  
                        jeta=ak4[mask].abseta.flatten(),   
                        weight=w_hfjets
                    )

                ezfill('ak4_hfcentral_adjacent_etastripsize',   
                    centraletastripsize=ak4[mask].hfcentralstripsize.flatten(),
                    adjacentetastripsize=ak4[mask].hfadjacentstripsize.flatten(),
                    jeta=ak4[mask].abseta.flatten(),
                    weight=w_hfjets
                    )

                # Leading and trailing jets
                hfmask_i0 = is_hf_jet(diak4.i0[mask])
                hfmask_i1 = is_hf_jet(diak4.i1[mask])

                w_hfjets_i0 = np.where(
                    hfmask_i0.flatten(),
                    w_diak4,
                    0.
                )

                w_hfjets_i1 = np.where(
                    hfmask_i1.flatten(),
                    w_diak4,
                    0.
                )

                ezfill('ak4_sigma_eta_phi0',   
                    sigmaetaeta=diak4.i0.setaeta[mask].flatten(),    
                    sigmaphiphi=diak4.i0.sphiphi[mask].flatten(),    
                    jeta=diak4.i0.abseta[mask].flatten(),   
                    weight=w_hfjets_i0
                )
                
                ezfill('ak4_sigma_eta_phi1',   
                    sigmaetaeta=diak4.i1.setaeta[mask].flatten(),    
                    sigmaphiphi=diak4.i1.sphiphi[mask].flatten(),    
                    jeta=diak4.i1.abseta[mask].flatten(),   
                    weight=w_hfjets_i1
                )
            
                ezfill('ak4_hfcentral_adjacent_etastripsize0',   
                    centraletastripsize=diak4.i0.hfcentralstripsize[mask].flatten(),
                    adjacentetastripsize=diak4.i0.hfadjacentstripsize[mask].flatten(),
                    jeta=diak4.i0.abseta[mask].flatten(),
                    weight=w_hfjets_i0
                )

                ezfill('ak4_hfcentral_adjacent_etastripsize1',   
                    centraletastripsize=diak4.i1.hfcentralstripsize[mask].flatten(),
                    adjacentetastripsize=diak4.i1.hfadjacentstripsize[mask].flatten(),
                    jeta=diak4.i1.abseta[mask].flatten(),
                    weight=w_hfjets_i1
                )

            # B tag discriminator
            btag = getattr(ak4, cfg.BTAG.ALGO)
            w_btag = weight_shape(btag[mask], rweight[mask])
            ezfill('ak4_btag', btag=btag[mask].flatten(), weight=w_btag )

            # MET
            ezfill('dpfcalo_cr',            dpfcalo=df["dPFCaloCR"][mask],       weight=rweight[mask] )
            ezfill('dpfcalo_sr',            dpfcalo=df["dPFCaloSR"][mask],       weight=rweight[mask] )
            ezfill('met',                met=met_pt[mask],            weight=rweight[mask] )
            ezfill('met_phi',            phi=met_phi[mask],           weight=rweight[mask] )
            ezfill('recoil',             recoil=df["recoil_pt"][mask],      weight=rweight[mask] )
            ezfill('recoil_phi',         phi=df["recoil_phi"][mask],        weight=rweight[mask] )
            ezfill('dphijm',             dphi=df["minDPhiJetMet"][mask],    weight=rweight[mask] )
            ezfill('dphijr',             dphi=df["minDPhiJetRecoil"][mask], weight=rweight[mask] )

            ezfill('dphijj',             dphi=df["dphijj"][mask],   weight=rweight[mask] )
            ezfill('detajj',             deta=df["detajj"][mask],   weight=rweight[mask] )
            ezfill('mjj',                mjj=df["mjj"][mask],      weight=rweight[mask] )

            rweight_nopref = region_weights.partial_weight(exclude=exclude+['prefire'])
            ezfill('mjj_nopref',                mjj=df["mjj"][mask],      weight=rweight_nopref[mask] )

            if not df['is_data']:
                ezfill('gen_met_mjj',        met=df["GenMET_pt"][mask],     mjj=df["mjj"][mask],      weight=rweight[mask])

            ezfill('vecdphi',     vecdphi=vec_dphi[mask],       weight=rweight[mask] )
            ezfill('vecb',        vecb=vec_b[mask],            weight=rweight[mask] )
            ezfill('dphitkpf',    dphi=dphitkpf[mask],         weight=rweight[mask] )

            if region != 'inclusive':
                ezfill('dphitkpf_ak4_eta0',  dphi=dphitkpf[mask],     jeteta=diak4.i0.abseta[mask].flatten(),     weight=rweight[mask])

            # Consider events where only (at least) one of the jets is in horn
            dpftkmet_weight = np.where(
                leading_jet_in_horn | trailing_jet_in_horn,
                rweight,
                0
            )

            ezfill('dPFTkMET',   dpftk=df['dPFTkSR'][mask],    weight=dpftkmet_weight[mask])

            # b-tag weight up and down variations
            if cfg.RUN.BTAG_STUDY:
                if not df['is_data']:
                    rw = region_weights.partial_weight(exclude=exclude+['bveto'])
                    ezfill('mjj_bveto_up',    mjj=df['mjj'][mask],  weight=(rw*(1-bsf_variations['up']).prod())[mask])
                    ezfill('mjj_bveto_down',  mjj=df['mjj'][mask],  weight=(rw*(1-bsf_variations['down']).prod())[mask])

            if gen_v_pt is not None:
                ezfill('gen_vpt', vpt=gen_v_pt[mask], weight=df['Generator_weight'][mask])
                ezfill('gen_mjj', mjj=df['mjj_gen'][mask], weight=df['Generator_weight'][mask])


            if cfg.RUN.ELE_SF_STUDY and re.match('cr_(\d)e_vbf', region) and not df['is_data']:
                eleloose_id_sf, eletight_id_sf, ele_reco_sf = get_varied_ele_sf(electrons, df, evaluator)
                rw = region_weights.partial_weight(exclude=exclude+['ele_id_tight','ele_id_loose'])
                for ele_id_variation in eletight_id_sf.keys():
                    ezfill('mjj_ele_id', 
                        mjj=df['mjj'][mask], 
                        variation=ele_id_variation,
                        weight=(rw * eletight_id_sf[ele_id_variation].prod() * eleloose_id_sf[ele_id_variation].prod())[mask],
                    )

                for ele_reco_variation, w in ele_reco_sf.items():
                    ezfill('mjj_ele_reco',
                        mjj=df['mjj'][mask],
                        variation=ele_reco_variation,
                        weight=(rw * w)[mask]
                    )

            if cfg.RUN.MUON_SF_STUDY and re.match('cr_(\d)m_vbf', region) and not df['is_data']:
                muon_looseid_sf, muon_tightid_sf, muon_looseiso_sf, muon_tightiso_sf = get_varied_muon_sf(muons, df, evaluator)
                rw = region_weights.partial_weight(exclude=exclude+['muon_id_tight','muon_id_loose'])

                for mu_id_variation in muon_looseid_sf.keys():
                    ezfill('mjj_muon_id', 
                        mjj=df['mjj'][mask], 
                        variation=mu_id_variation,
                        weight=(rw * muon_tightid_sf[mu_id_variation].prod() * muon_looseid_sf[mu_id_variation].prod())[mask],
                    )

                for mu_iso_variation in muon_looseiso_sf.keys():
                    ezfill('mjj_muon_iso', 
                        mjj=df['mjj'][mask], 
                        variation=mu_iso_variation,
                        weight=(rw * muon_tightiso_sf[mu_iso_variation].prod() * muon_looseiso_sf[mu_iso_variation].prod())[mask],
                    )

            if cfg.RUN.ELE_TRIG_STUDY and not df['is_data'] and re.match('cr_(\d)e_vbf', region):
                # Note that electrons in the gap do not count in this study
                mask_electron_nogap = (np.abs(electrons.etasc)<1.4442) | (np.abs(electrons.etasc)>1.566)
                electrons_nogap = electrons[mask_electron_nogap]

                # Up and down variations: Vary the efficiency in data
                data_eff_up = evaluator['trigger_electron_eff_data'](electrons_nogap.etasc, electrons_nogap.pt) + evaluator['trigger_electron_eff_data_error'](electrons_nogap.etasc, electrons_nogap.pt)
                data_eff_down = evaluator['trigger_electron_eff_data'](electrons_nogap.etasc, electrons_nogap.pt) - evaluator['trigger_electron_eff_data_error'](electrons_nogap.etasc, electrons_nogap.pt)

                p_pass_data = 1 - (1-evaluator["trigger_electron_eff_data"](electrons_nogap.etasc, electrons_nogap.pt)).prod()
                p_pass_data_up = 1 - (1-data_eff_up).prod()
                p_pass_data_down = 1 - (1-data_eff_down).prod()

                p_pass_mc = 1 - (1-evaluator["trigger_electron_eff_mc"](electrons_nogap.etasc, electrons_nogap.pt)).prod()

                trigger_weight_nom = p_pass_data / p_pass_mc
                trigger_weight_up = p_pass_data_up / p_pass_mc
                trigger_weight_down = p_pass_data_down / p_pass_mc

                trigger_weight_nom[np.isnan(trigger_weight_nom) | np.isinf(trigger_weight_nom)] = 1.
                trigger_weight_up[np.isnan(trigger_weight_up) | np.isinf(trigger_weight_up)] = 1.
                trigger_weight_down[np.isnan(trigger_weight_down) | np.isinf(trigger_weight_down)] = 1.

                ele_trig_sf = {
                    "nom" : trigger_weight_nom,
                    "up" : trigger_weight_up,
                    "down" : trigger_weight_down,
                }

                for variation, trigw in ele_trig_sf.items():
                    rw = region_weights.partial_weight(exclude=exclude+['trigger_ele'])
                    ezfill(
                        'mjj_ele_trig_weight',
                        mjj=df['mjj'][mask],
                        variation=variation,
                        weight=(rw*trigw)[mask] 
                    )

            if cfg.RUN.PILEUP_SF_STUDY and not df['is_data']:
                rw_nopu = region_weights.partial_weight(exclude=exclude+['pileup'])

                puweights = pileup_sf_variations(df, evaluator, cfg)
                for puvar, w in puweights.items():
                    ezfill('mjj_pu_weights',
                        mjj=df['mjj'][mask],
                        variation=puvar,
                        weight=(rw_nopu * w)[mask]
                    )

            if cfg.RUN.VETO_WEIGHTS_STUDY and 'no_veto_all' in region and not df['is_data']:
                variations = ['nominal', 'tau_id_up', 'tau_id_dn', 
                    'ele_id_up', 'ele_id_dn', 'ele_reco_up', 'ele_reco_dn',
                    'muon_id_up', 'muon_id_dn', 'muon_iso_up', 'muon_iso_dn'
                    ]
                rw_no_veto = region_weights.partial_weight(exclude=exclude+['veto'])
                for v in variations:
                    ezfill('mjj_veto_weight',
                        mjj=df['mjj'][mask],
                        variation=v,
                        weight=(rw_no_veto * veto_weights.partial_weight(include=[v]))[mask]
                    )


            # Photon CR data-driven QCD estimate
            if df['is_data'] and re.match("cr_g.*", region) and re.match("(SinglePhoton|EGamma).*", dataset):
                w_imp = photon_impurity_weights(photons[leadphoton_index].pt.max()[mask], df["year"])
                output['mjj'].fill(
                                    dataset=data_driven_qcd_dataset(dataset),
                                    region=region,
                                    mjj=df["mjj"][mask],
                                    weight=rweight[mask] * w_imp
                                )
                output['recoil'].fill(
                                    dataset=data_driven_qcd_dataset(dataset),
                                    region=region,
                                    recoil=df["recoil_pt"][mask],
                                    weight=rweight[mask] * w_imp
                                )

            # Uncertainty variations
            if df['is_lo_z'] or df['is_nlo_z'] or df['is_lo_z_ewk']:
                theory_uncs = [x for x in cfg.SF.keys() if x.startswith('unc')]
                for unc in theory_uncs:
                    reweight = evaluator[unc](gen_v_pt)
                    w = (region_weights.weight() * reweight)[mask]
                    ezfill(
                        'mjj_unc',
                        mjj=df['mjj'][mask],
                        uncertainty=unc,
                        weight=w)

            # Two dimensional
            ezfill('recoil_mjj',         recoil=df["recoil_pt"][mask], mjj=df["mjj"][mask], weight=rweight[mask] )

            # Muons
            if '_1m_' in region or '_2m_' in region or 'no_veto' in region:
                w_allmu = weight_shape(muons.pt[mask], rweight[mask])
                ezfill('muon_pt',   pt=muons.pt[mask].flatten(),    weight=w_allmu )
                ezfill('muon_pt_abseta',pt=muons.pt[mask].flatten(),abseta=muons.eta[mask].flatten(),    weight=w_allmu )
                ezfill('muon_mt',   mt=df['MT_mu'][mask],           weight=rweight[mask])
                ezfill('muon_eta',  eta=muons.eta[mask].flatten(),  weight=w_allmu)
                ezfill('muon_phi',  phi=muons.phi[mask].flatten(),  weight=w_allmu)

            # Dimuon
            if '_2m_' in region:
                w_dimu = weight_shape(dimuons.pt[mask], rweight[mask])
                ezfill('muon_pt0',      pt=dimuons.i0.pt[mask].flatten(),           weight=w_dimu)
                ezfill('muon_pt1',      pt=dimuons.i1.pt[mask].flatten(),           weight=w_dimu)
                ezfill('muon_eta0',     eta=dimuons.i0.eta[mask].flatten(),         weight=w_dimu)
                ezfill('muon_eta1',     eta=dimuons.i1.eta[mask].flatten(),         weight=w_dimu)
                ezfill('muon_phi0',     phi=dimuons.i0.phi[mask].flatten(),         weight=w_dimu)
                ezfill('muon_phi1',     phi=dimuons.i1.phi[mask].flatten(),         weight=w_dimu)
                ezfill('dimuon_pt',     pt=dimuons.pt[mask].flatten(),              weight=w_dimu)
                ezfill('dimuon_eta',    eta=dimuons.eta[mask].flatten(),            weight=w_dimu)
                ezfill('dimuon_mass',   dilepton_mass=dimuons.mass[mask].flatten(), weight=w_dimu )

            # Electrons
            if '_1e_' in region or '_2e_' in region or 'no_veto' in region:
                w_allel = weight_shape(electrons.pt[mask], rweight[mask])
                ezfill('electron_pt',   pt=electrons.pt[mask].flatten(),    weight=w_allel)
                ezfill('electron_pt_eta',   pt=electrons.pt[mask].flatten(), eta=electrons.eta[mask].flatten(),    weight=w_allel)
                ezfill('electron_mt',   mt=df['MT_el'][mask],               weight=rweight[mask])
                ezfill('electron_eta',  eta=electrons.eta[mask].flatten(),  weight=w_allel)
                ezfill('electron_phi',  phi=electrons.phi[mask].flatten(),  weight=w_allel)

            # Dielectron
            if '_2e_' in region:
                w_diel = weight_shape(dielectrons.pt[mask], rweight[mask])
                ezfill('electron_pt0',      pt=dielectrons.i0.pt[mask].flatten(),               weight=w_diel)
                ezfill('electron_pt1',      pt=dielectrons.i1.pt[mask].flatten(),               weight=w_diel)
                ezfill('electron_eta0',     eta=dielectrons.i0.eta[mask].flatten(),             weight=w_diel)
                ezfill('electron_eta1',     eta=dielectrons.i1.eta[mask].flatten(),             weight=w_diel)
                ezfill('electron_phi0',     phi=dielectrons.i0.phi[mask].flatten(),             weight=w_diel)
                ezfill('electron_phi1',     phi=dielectrons.i1.phi[mask].flatten(),             weight=w_diel)
                ezfill('dielectron_pt',     pt=dielectrons.pt[mask].flatten(),                  weight=w_diel)
                ezfill('dielectron_eta',    eta=dielectrons.eta[mask].flatten(),                weight=w_diel)
                ezfill('dielectron_mass',   dilepton_mass=dielectrons.mass[mask].flatten(),     weight=w_diel)

            # Photon
            if '_g_' in region:
                w_leading_photon = weight_shape(photons[leadphoton_index].pt[mask],rweight[mask]);
                ezfill('photon_pt0',              pt=photons[leadphoton_index].pt[mask].flatten(),    weight=w_leading_photon)
                ezfill('photon_eta0',             eta=photons[leadphoton_index].eta[mask].flatten(),  weight=w_leading_photon)
                ezfill('photon_phi0',             phi=photons[leadphoton_index].phi[mask].flatten(),  weight=w_leading_photon)
                ezfill('photon_pt0_recoil',       pt=photons[leadphoton_index].pt[mask].flatten(), recoil=df['recoil_pt'][mask&(leadphoton_index.counts>0)],  weight=w_leading_photon)
                ezfill('photon_eta_phi',          eta=photons[leadphoton_index].eta[mask].flatten(), phi=photons[leadphoton_index].phi[mask].flatten(),  weight=w_leading_photon)

                # w_drphoton_jet = weight_shape(df['dRPhotonJet'][mask], rweight[mask])

            # Tau
            if 'no_veto' in region:
                w_all_taus = weight_shape(taus.pt[mask], rweight[mask])
                ezfill("tau_pt", pt=taus.pt[mask].flatten(), weight=w_all_taus)

            # PV
            ezfill('npv', nvtx=df['PV_npvs'][mask], weight=rweight[mask])
            ezfill('npvgood', nvtx=df['PV_npvsGood'][mask], weight=rweight[mask])

            ezfill('npv_nopu', nvtx=df['PV_npvs'][mask], weight=region_weights.partial_weight(exclude=exclude+['pileup'])[mask])
            ezfill('npvgood_nopu', nvtx=df['PV_npvsGood'][mask], weight=region_weights.partial_weight(exclude=exclude+['pileup'])[mask])

            ezfill('rho_all', rho=df['fixedGridRhoFastjetAll'][mask], weight=region_weights.partial_weight(exclude=exclude)[mask])
            ezfill('rho_central', rho=df['fixedGridRhoFastjetCentral'][mask], weight=region_weights.partial_weight(exclude=exclude)[mask])
            ezfill('rho_all_nopu', rho=df['fixedGridRhoFastjetAll'][mask], weight=region_weights.partial_weight(exclude=exclude+['pileup'])[mask])
            ezfill('rho_central_nopu', rho=df['fixedGridRhoFastjetCentral'][mask], weight=region_weights.partial_weight(exclude=exclude+['pileup'])[mask])
        return output

    def postprocess(self, accumulator):
        return accumulator
