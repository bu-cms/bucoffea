import os
import numpy as np
import re
from coffea import hist
import coffea.processor as processor

from dynaconf import settings as cfg

from bucoffea.monojet.definitions import monojet_accumulator, setup_candidates, monojet_regions
from bucoffea.helpers import min_dphi_jet_met, recoil, mt, weight_shape, bucoffea_path, dphi,mask_and, mask_or, evaluator_from_config
from bucoffea.helpers.dataset import is_lo_z, is_lo_w, is_lo_g, is_nlo_z, is_nlo_w, is_data, extract_year
from bucoffea.helpers.gen import find_gen_dilepton, setup_gen_candidates

def trigger_selection(selection, df, cfg):
    pass_all = np.zeros(df.size) == 0
    pass_none = ~pass_all
    dataset = df['dataset']
    if cfg.RUN.SYNC: # Synchronization mode
        selection.add('filt_met', pass_all)
        selection.add('trig_met', pass_all)
        selection.add('trig_ele', pass_all)
        selection.add('trig_mu',  pass_all)
        selection.add('trig_photon',  pass_all)

    else:
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

class monojetProcessor(processor.ProcessorABC):
    def __init__(self, blind=True):
        self._year=None
        self._blind=blind
        self._configure()
        self._accumulator = monojet_accumulator(cfg)

    @property
    def accumulator(self):
        return self._accumulator

    def _configure(self, df=None):
        cfg.DYNACONF_WORKS="merge_configs"
        cfg.MERGE_ENABLED_FOR_DYNACONF=True
        cfg.SETTINGS_FILE_FOR_DYNACONF = bucoffea_path("config/monojet.yaml")

        # Reload config based on year
        if df:
            dataset = df['dataset']
            self._year = extract_year(dataset)
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
        df['is_lo_g'] = is_lo_g(dataset)
        df['is_nlo_z'] = is_nlo_z(dataset)
        df['is_nlo_w'] = is_nlo_w(dataset)
        df['has_lhe_v_pt'] = df['is_lo_w'] | df['is_lo_z'] | df['is_nlo_z'] | df['is_nlo_w'] | df['is_lo_g']
        df['is_data'] = is_data(dataset)

        gen_v_pt = None
        if df['is_lo_w'] or df['is_lo_z'] or df['is_nlo_z'] or df['is_nlo_w']:
            gen = setup_gen_candidates(df)
            if is_lo_z(dataset) or is_nlo_z(dataset):
                pdgsum = 0
            elif is_lo_w(dataset) or is_nlo_w(dataset):
                pdgsum = 1
            gen_dilep = find_gen_dilepton(gen, pdgsum)
            gen_v_pt = gen_dilep[gen_dilep.mass.argmax()].pt.max()
        elif df['is_lo_g']:
            gen_v_pt = df['LHE_Vpt']

        # Candidates
        # Already pre-filtered!
        # All leptons are at least loose
        # Check out setup_candidates for filtering details
        met_pt, met_phi, ak4, ak8, muons, electrons, taus, photons, hlt = setup_candidates(df, cfg)

        # Muons
        is_tight_muon = muons.tightId \
                      & (muons.iso < cfg.MUON.CUTS.TIGHT.ISO) \
                      & (muons.pt>cfg.MUON.CUTS.TIGHT.PT) \
                      & (np.abs(muons.eta)<cfg.MUON.CUTS.TIGHT.ETA)

        dimuons = muons.distincts()
        dimuon_charge = dimuons.i0['charge'] + dimuons.i1['charge']

        df['MT_mu'] = ((muons.counts==1) * mt(muons.pt, muons.phi, met_pt, met_phi)).max()

        # Electrons
        is_tight_electron = electrons.tightId \
                            & (electrons.pt > cfg.ELECTRON.CUTS.TIGHT.PT) \
                            & (np.abs(electrons.eta) < cfg.ELECTRON.CUTS.TIGHT.ETA)

        dielectrons = electrons.distincts()
        dielectron_charge = dielectrons.i0['charge'] + dielectrons.i1['charge']

        df['MT_el'] = ((electrons.counts==1) * mt(electrons.pt, electrons.phi, met_pt, met_phi)).max()

        # ak4
        jet_acceptance = np.abs(ak4.eta)<2.4
        leadak4_index=ak4.pt.argmax()

        elejet_pairs = ak4[:,:1].cross(electrons)
        df['dREleJet'] = np.hypot(elejet_pairs.i0.eta-elejet_pairs.i1.eta , dphi(elejet_pairs.i0.phi,elejet_pairs.i1.phi)).min()
        muonjet_pairs = ak4[:,:1].cross(muons)
        df['dRMuonJet'] = np.hypot(muonjet_pairs.i0.eta-muonjet_pairs.i1.eta , dphi(muonjet_pairs.i0.phi,muonjet_pairs.i1.phi)).min()

        # B tagged ak4
        btag_cut = cfg.BTAG.CUTS[cfg.BTAG.algo][cfg.BTAG.wp]
        jet_btag_val = getattr(ak4, cfg.BTAG.algo)
        jet_btagged = jet_btag_val > btag_cut
        bjets = ak4[ jet_acceptance \
                     & jet_btagged \
                     & (ak4.pt>20) ]

        # Photons
        # Angular distance leading photon - leading jet
        phojet_pairs = ak4[:,:1].cross(photons[:,:1])
        df['dRPhotonJet'] = np.hypot(phojet_pairs.i0.eta-phojet_pairs.i1.eta , dphi(phojet_pairs.i0.phi,phojet_pairs.i1.phi)).min()

        # Recoil
        df['recoil_pt'], df['recoil_phi'] = recoil(met_pt,met_phi, electrons, muons, photons)
        df["dPFCalo"] = (met_pt - df["CaloMET_pt"]) / df["recoil_pt"]
        df["minDPhiJetRecoil"] = min_dphi_jet_met(ak4, df['recoil_phi'], njet=4, ptmin=30)
        df["minDPhiJetMet"] = min_dphi_jet_met(ak4, met_phi, njet=4, ptmin=30)
        selection = processor.PackedSelection()



        # Triggers
        pass_all = np.ones(df.size)==1
        pass_none = ~pass_all
        selection.add('inclusive', pass_all)
        selection = trigger_selection(selection, df, cfg)

        # Trigger objects
        hlt_single_muons = hlt[(hlt.id==13) & (hlt.filter & 8 == 8)]
        muons_hltmatch = muons[muons.match(hlt_single_muons,deltaRCut=0.2,deltaPtCut=0.25)]
        selection.add('one_hlt_muon', muons_hltmatch.counts>=1)
        selection.add('two_hlt_muons', muons_hltmatch.counts==2)
        selection.add('mu_pt_trig_safe', muons.pt.max() > 30)

        # Common selection
        selection.add('veto_ele', electrons.counts==0)
        selection.add('veto_muo', muons.counts==0)
        selection.add('veto_photon', photons.counts==0)
        selection.add('veto_tau', taus.counts==0)
        selection.add('veto_b', bjets.counts==0)
        selection.add('mindphijr',df['minDPhiJetRecoil'] > cfg.SELECTION.SIGNAL.MINDPHIJR)
        selection.add('mindphijm',df['minDPhiJetMet'] > cfg.SELECTION.SIGNAL.MINDPHIJR)
        selection.add('dpfcalo',np.abs(df['dPFCalo']) < cfg.SELECTION.SIGNAL.DPFCALO)
        selection.add('recoil', df['recoil_pt']>cfg.SELECTION.SIGNAL.RECOIL)

        if(cfg.MITIGATION.HEM and extract_year(df['dataset']) == 2018 and not cfg.RUN.SYNC):
            selection.add('hemveto', df['hemveto'])
        else:
            selection.add('hemveto', np.ones(df.size)==1)

        # AK4 Jet
        leadak4_pt_eta = (ak4.pt.max() > cfg.SELECTION.SIGNAL.leadak4.PT) \
                         & (np.abs(ak4.eta[leadak4_index]) < cfg.SELECTION.SIGNAL.leadak4.ETA).any()
        selection.add('leadak4_pt_eta', leadak4_pt_eta)

        selection.add('leadak4_id',(ak4.tightId[leadak4_index] \
                                                    & (ak4.chf[leadak4_index] >cfg.SELECTION.SIGNAL.leadak4.CHF) \
                                                    & (ak4.nhf[leadak4_index]<cfg.SELECTION.SIGNAL.leadak4.NHF)).any())

        # AK8 Jet
        leadak8_index=ak8.pt.argmax()
        leadak8_pt_eta = (ak8.pt.max() > cfg.SELECTION.SIGNAL.leadak8.PT) \
                         & (np.abs(ak8.eta[leadak8_index]) < cfg.SELECTION.SIGNAL.leadak8.ETA).any()
        selection.add('leadak8_pt_eta', leadak8_pt_eta)

        selection.add('leadak8_id',(ak8.tightId[leadak8_index]).any())

        # Mono-V selection
        selection.add('leadak8_tau21', ((ak8.tau2[leadak8_index] / ak8.tau1[leadak8_index]) < cfg.SELECTION.SIGNAL.LEADAK8.TAU21).any())
        selection.add('leadak8_mass', ((ak8.mass[leadak8_index] > cfg.SELECTION.SIGNAL.LEADAK8.MASS.MIN) \
                                    & (ak8.mass[leadak8_index] < cfg.SELECTION.SIGNAL.LEADAK8.MASS.MAX)).any())

        selection.add('veto_vtag', ~selection.all("leadak8_pt_eta", "leadak8_id", "leadak8_tau21", "leadak8_mass"))

        # Dimuon CR
        leadmuon_index=muons.pt.argmax()
        selection.add('at_least_one_tight_mu', is_tight_muon.any())
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
        selection.add('at_least_one_tight_el', is_tight_electron.any())

        selection.add('dielectron_mass', ((dielectrons.mass > cfg.SELECTION.CONTROL.DOUBLEEL.MASS.MIN)  \
                                        & (dielectrons.mass < cfg.SELECTION.CONTROL.DOUBLEEL.MASS.MAX)).any())
        selection.add('dielectron_charge', (dielectron_charge==0).any())
        selection.add('two_electrons', electrons.counts==2)

        # Single Ele CR
        selection.add('met_el', met_pt > cfg.SELECTION.CONTROL.SINGLEEL.MET)
        selection.add('mt_el', df['MT_el'] < cfg.SELECTION.CONTROL.SINGLEEL.MT)

        # Photon CR
        leadphoton_index=photons.pt.argmax()

        is_tight_photon = photons.mediumId \
                         & (np.abs(photons.eta) < cfg.PHOTON.CUTS.TIGHT.ETA)

        selection.add('one_photon', photons.counts==1)
        selection.add('at_least_one_tight_photon', is_tight_photon.any())
        selection.add('photon_pt', photons.pt.max() > cfg.PHOTON.CUTS.TIGHT.PT)
        selection.add('photon_pt_trig', photons.pt.max() > cfg.PHOTON.CUTS.TIGHT.PTTRIG)

        # Fill histograms
        output = self.accumulator.identity()

        # Gen
        if gen_v_pt is not None:
            output['genvpt_check'].fill(vpt=gen_v_pt,type="Nano", dataset=dataset, weight=df['Generator_weight'])

        if 'LHE_Njets' in df:
            output['lhe_njets'].fill(dataset=dataset, multiplicity=df['LHE_Njets'])
        if 'LHE_HT' in df:
            output['lhe_ht'].fill(dataset=dataset, ht=df['LHE_HT'])
        if 'LHE_HTIncoming' in df:
            output['lhe_htinc'].fill(dataset=dataset, ht=df['LHE_HTIncoming'])

        # Weights
        evaluator = evaluator_from_config(cfg)


        weight = np.ones(df.size)
        weight_nopu = np.ones(df.size)
        weight_nopref = np.ones(df.size)
        all_weights = {}
        if not df['is_data']:
            all_weights['gen'] = df['Generator_weight']

            try:
                all_weights['prefire'] = df['PrefireWeight']
            except KeyError:
                all_weights['prefire'] = np.ones(df.size)

            # Muon ID and Isolation for tight and loose WP
            # Function of pT, eta (Order!)
            all_weights["muon_id_tight"] = evaluator['muon_id_tight'](muons[is_tight_muon].pt, muons[is_tight_muon].eta).prod()
            all_weights["muon_iso_tight"] = evaluator['muon_iso_tight'](muons[is_tight_muon].pt, muons[is_tight_muon].eta).prod()
            all_weights["muon_id_loose"] = evaluator['muon_id_loose'](muons[~is_tight_muon].pt, muons[~is_tight_muon].eta).prod()
            all_weights["muon_iso_loose"] = evaluator['muon_iso_loose'](muons[~is_tight_muon].pt, muons[~is_tight_muon].eta).prod()

            # Electron ID and reco
            # Function of eta, pT (Other way round relative to muons!)
            all_weights["ele_reco"] = evaluator['ele_reco'](electrons.eta, electrons.pt).prod()
            all_weights["ele_id_tight"] = evaluator['ele_id_tight'](electrons[is_tight_electron].eta, electrons[is_tight_electron].pt).prod()
            all_weights["ele_id_loose"] = evaluator['ele_id_loose'](electrons[~is_tight_electron].eta, electrons[~is_tight_electron].pt).prod()

            # Photon ID and electron veto
            all_weights["photon_id_tight"] = evaluator['photon_id_tight'](photons[is_tight_photon].eta, photons[is_tight_photon].pt).prod()

            # CSEV not split only by EE/EB for now
            csev_sf_index = 0.5 * photons.barrel + 3.5 * ~photons.barrel + 1 * (photons.r9 > 0.94) + 2 * (photons.r9 <= 0.94)
            all_weights["photon_csev"] = evaluator['photon_csev'](csev_sf_index).prod()

            if cfg.SF.PILEUP.MODE == 'nano':
                all_weights["pileup"] = df['puWeight']
            elif cfg.SF.PILEUP.MODE == 'manual':
                all_weights["pileup"] = evaluator['pileup'](df['Pileup_nTrueInt'])
            else:
                raise RuntimeError(f"Unknown value for cfg.PILEUP.MODE: {cfg.PILEUP.MODE}.")

            if df['is_lo_w']:
                all_weights["theory"] = evaluator["qcd_nlo_w_2017"](gen_v_pt) * evaluator["qcd_nnlo_w"](gen_v_pt) * evaluator["ewk_nlo_w"](gen_v_pt)
            elif df['is_lo_z']:
                all_weights["theory"] = evaluator["qcd_nlo_z_2017"](gen_v_pt) * evaluator["qcd_nnlo_z"](gen_v_pt) * evaluator["ewk_nlo_z"](gen_v_pt)
            elif df['is_nlo_w']:
                all_weights["theory"] = evaluator["qcd_nnlo_w"](gen_v_pt) * evaluator["ewk_nlo_w"](gen_v_pt)
            elif df['is_nlo_z']:
                all_weights["theory"] = evaluator["qcd_nnlo_z"](gen_v_pt) * evaluator["ewk_nlo_z"](gen_v_pt)
            elif df['is_lo_g']:
                all_weights["theory"] = evaluator["ewk_nlo_g"](gen_v_pt) * evaluator["qcd_nlo_g"](gen_v_pt) * evaluator["qcd_nnlo_g"](gen_v_pt)
            else:
                all_weights["theory"] = np.ones(df.size)

            for name, iw in all_weights.items():
                weight = weight * iw
                if name != 'pileup':
                    weight_nopu = weight_nopu * iw
                if name != 'prefire':
                    weight_nopref = weight_nopref * iw

        # Save per-event values for synchronization
        if cfg.RUN.KINEMATICS.SAVE:
            for event in cfg.RUN.KINEMATICS.EVENTS:
                mask = df['event'] == event
                if not mask.any():
                    continue
                output['kinematics']['event'] += [event]
                output['kinematics']['met'] += [met_pt[mask].flatten()]
                output['kinematics']['met_phi'] += [met_phi[mask].flatten()]
                output['kinematics']['recoil'] += [df['recoil_pt'][mask].flatten()]
                output['kinematics']['recoil_phi'] += [df['recoil_phi'][mask].flatten()]

                output['kinematics']['ak4pt0'] += [ak4[leadak4_index][mask].pt.flatten()]
                output['kinematics']['ak4eta0'] += [ak4[leadak4_index][mask].eta.flatten()]
                output['kinematics']['leadbtag'] += [jet_btag_val[jet_acceptance & (ak4.pt>20)][mask].max()]

                output['kinematics']['nLooseMu'] += [muons.counts[mask]]
                output['kinematics']['nTightMu'] += [muons[is_tight_muon].counts[mask].flatten()]
                output['kinematics']['mupt0'] += [muons[leadmuon_index][mask].pt.flatten()]
                output['kinematics']['mueta0'] += [muons[leadmuon_index][mask].eta.flatten()]

                output['kinematics']['nLooseEl'] += [electrons.counts[mask]]
                output['kinematics']['nTightEl'] += [electrons[is_tight_electron].counts[mask].flatten()]
                output['kinematics']['elpt0'] += [electrons[leadelectron_index][mask].pt.flatten()]
                output['kinematics']['eleta0'] += [electrons[leadelectron_index][mask].eta.flatten()]

                output['kinematics']['nLooseGam'] += [photons.counts[mask]]
                output['kinematics']['nTightGam'] += [photons[is_tight_photon].counts[mask].flatten()]
                output['kinematics']['gpt0'] += [photons[leadphoton_index][mask].pt.flatten()]
                output['kinematics']['geta0'] += [photons[leadphoton_index][mask].eta.flatten()]


        # Sum of all weights to use for normalization
        # TODO: Deal with systematic variations
        output['nevents'][dataset] += df.size
        if not df['is_data']:
            output['sumw'][dataset] +=  df['genEventSumw']
            output['sumw2'][dataset] +=  df['genEventSumw2']
            output['sumw_pileup'][dataset] +=  all_weights['pileup'].sum()

        regions = monojet_regions(cfg)
        for region, cuts in regions.items():
            # Blinding
            if(self._blind and df['is_data'] and region.startswith('sr')):
                continue

            # Cutflow plot for signal and control regions
            if any(x in region for x in ["sr", "cr", "tr"]):
                output['cutflow_' + region]['all']+=df.size
                for icut, cutname in enumerate(cuts):
                    output['cutflow_' + region][cutname] += selection.all(*cuts[:icut+1]).sum()

            mask = selection.all(*cuts)


            # Save the event numbers of events passing this selection
            if cfg.RUN.SAVE.PASSING:
                output['selected_events'][region] += list(df['event'][mask])


            # Multiplicities
            def fill_mult(name, candidates):
                output[name].fill(
                                  dataset=dataset,
                                  region=region,
                                  multiplicity=candidates[mask].counts,
                                  weight=weight[mask]
                                  )

            fill_mult('ak8_mult', ak8)
            fill_mult('ak4_mult', ak4)
            fill_mult('bjet_mult',bjets)
            fill_mult('loose_ele_mult',electrons)
            fill_mult('tight_ele_mult',electrons[is_tight_electron])
            fill_mult('loose_muo_mult',muons)
            fill_mult('tight_muo_mult',muons[is_tight_muon])
            fill_mult('tau_mult',taus)
            fill_mult('photon_mult',photons)
            fill_mult('hlt_single_muon_mult',hlt_single_muons)
            fill_mult('muons_hltmatch_mult',muons_hltmatch)

            def ezfill(name, **kwargs):
                """Helper function to make filling easier."""
                output[name].fill(
                                  dataset=dataset,
                                  region=region,
                                  **kwargs
                                  )
            # Monitor weights
            for wname, wvalue in all_weights.items():
                ezfill("weights", weight_type=wname, weight_value=wvalue[mask])
                ezfill("weights_wide", weight_type=wname, weight_value=wvalue[mask])

            # All ak4
            # This is a workaround to create a weight array of the right dimension
            w_alljets = weight_shape(ak4[mask].eta, weight[mask])
            w_alljets_nopref = weight_shape(ak4[mask].eta, weight_nopref[mask])

            ezfill('ak4_eta',    jeteta=ak4[mask].eta.flatten(), weight=w_alljets)
            ezfill('ak4_phi',    jetphi=ak4[mask].phi.flatten(), weight=w_alljets)
            ezfill('ak4_pt',     jetpt=ak4[mask].pt.flatten(),   weight=w_alljets)

            ezfill('ak4_eta_nopref',    jeteta=ak4[mask].eta.flatten(), weight=w_alljets_nopref)
            ezfill('ak4_phi_nopref',    jetphi=ak4[mask].phi.flatten(), weight=w_alljets_nopref)
            ezfill('ak4_pt_nopref',     jetpt=ak4[mask].pt.flatten(),   weight=w_alljets_nopref)

            # Leading ak4
            w_leadak4 = weight_shape(ak4[leadak4_index].eta[mask], weight[mask])
            ezfill('ak4_eta0',   jeteta=ak4[leadak4_index].eta[mask].flatten(),    weight=w_leadak4)
            ezfill('ak4_phi0',   jetphi=ak4[leadak4_index].phi[mask].flatten(),    weight=w_leadak4)
            ezfill('ak4_pt0',    jetpt=ak4[leadak4_index].pt[mask].flatten(),      weight=w_leadak4)
            ezfill('ak4_ptraw0',    jetpt=ak4[leadak4_index].ptraw[mask].flatten(),      weight=w_leadak4)
            ezfill('ak4_chf0',    frac=ak4[leadak4_index].chf[mask].flatten(),      weight=w_leadak4)
            ezfill('ak4_nhf0',    frac=ak4[leadak4_index].nhf[mask].flatten(),      weight=w_leadak4)
            ezfill('ak4_nconst0',    nconst=ak4[leadak4_index].nconst[mask].flatten(),      weight=w_leadak4)

            ezfill('drelejet',    dr=df['dREleJet'][mask],      weight=weight[mask])
            ezfill('drmuonjet',    dr=df['dRMuonJet'][mask],      weight=weight[mask])
            ezfill('drphotonjet',    dr=df['dRPhotonJet'][mask],  weight=weight[mask])

            # Two-dimensional
            ezfill('ak4_pt0_eta0', jetpt=ak4[leadak4_index].pt[mask].flatten(), jeteta=ak4[leadak4_index].eta[mask].flatten(), weight=w_leadak4)
            ezfill('ak4_pt0_chf0', jetpt=ak4[leadak4_index].pt[mask].flatten(), frac=ak4[leadak4_index].chf[mask].flatten(), weight=w_leadak4)
            ezfill('ak4_pt0_nhf0', jetpt=ak4[leadak4_index].pt[mask].flatten(), frac=ak4[leadak4_index].nhf[mask].flatten(), weight=w_leadak4)
            ezfill('ak4_pt0_nconst0', jetpt=ak4[leadak4_index].pt[mask].flatten(), nconst=ak4[leadak4_index].nconst[mask].flatten(), weight=w_leadak4)

            # AK8 jets
            if region=='inclusive' or region.endswith('v'):
                # All
                w_allak8 = weight_shape(ak8.eta[mask], weight[mask])

                ezfill('ak8_eta',    jeteta=ak8[mask].eta.flatten(), weight=w_allak8)
                ezfill('ak8_phi',    jetphi=ak8[mask].phi.flatten(), weight=w_allak8)
                ezfill('ak8_pt',     jetpt=ak8[mask].pt.flatten(),   weight=w_allak8)
                ezfill('ak8_mass',   mass=ak8[mask].mass.flatten(),  weight=w_allak8)

                # Leading
                w_leadak8 = weight_shape(ak8[leadak8_index].eta[mask], weight[mask])

                ezfill('ak8_eta0',       jeteta=ak8[leadak8_index].eta[mask].flatten(),    weight=w_leadak8)
                ezfill('ak8_phi0',       jetphi=ak8[leadak8_index].phi[mask].flatten(),    weight=w_leadak8)
                ezfill('ak8_pt0',        jetpt=ak8[leadak8_index].pt[mask].flatten(),      weight=w_leadak8 )
                ezfill('ak8_mass0',      mass=ak8[leadak8_index].mass[mask].flatten(),     weight=w_leadak8)
                ezfill('ak8_tau210',     tau21=ak8[leadak8_index].tau21[mask].flatten(),     weight=w_leadak8)
                ezfill('ak8_wvsqcd0',    tagger=ak8[leadak8_index].wvsqcdmd[mask].flatten(),     weight=w_leadak8)
                ezfill('ak8_wvsqcdmd0',  tagger=ak8[leadak8_index].wvsqcd[mask].flatten(),     weight=w_leadak8)
                ezfill('ak8_zvsqcd0',    tagger=ak8[leadak8_index].zvsqcdmd[mask].flatten(),     weight=w_leadak8)
                ezfill('ak8_zvsqcdmd0',  tagger=ak8[leadak8_index].zvsqcd[mask].flatten(),     weight=w_leadak8)

            # B tag discriminator
            btag = getattr(ak4[jet_acceptance], cfg.BTAG.ALGO)
            w_btag = weight_shape(btag[mask], weight[mask])
            ezfill('ak4_btag', btag=btag[mask].flatten(), weight=w_btag )

            # MET
            ezfill('dpfcalo',            dpfcalo=df["dPFCalo"][mask],       weight=weight[mask] )
            ezfill('met',                met=met_pt[mask],            weight=weight[mask] )
            ezfill('met_phi',            phi=met_phi[mask],            weight=weight[mask] )
            ezfill('met_noweight',       met=met_pt[mask],            weight=np.ones(weight[mask].size) )
            ezfill('recoil',             recoil=df["recoil_pt"][mask],      weight=weight[mask] )
            ezfill('recoil_phi',         phi=df["recoil_phi"][mask],      weight=weight[mask] )
            ezfill('recoil_noweight',    recoil=df["recoil_pt"][mask],      weight=np.ones(weight[mask].size) )
            ezfill('dphijm',             dphi=df["minDPhiJetMet"][mask],    weight=weight[mask] )
            ezfill('dphijr',             dphi=df["minDPhiJetRecoil"][mask],    weight=weight[mask] )

            # Muons
            if '_1m_' in region or '_2m_' in region:
                w_allmu = weight_shape(muons.pt[mask], weight[mask])
                ezfill('muon_pt',   pt=muons.pt[mask].flatten(),    weight=w_allmu )
                ezfill('muon_mt',   mt=df['MT_mu'][mask],           weight=weight[mask])
                ezfill('muon_eta',  eta=muons.eta[mask].flatten(),  weight=w_allmu)
                ezfill('muon_phi',  phi=muons.phi[mask].flatten(),  weight=w_allmu)
                ezfill('muon_dxy',  dxy=muons.dxy[mask].flatten(),  weight=w_allmu)
                ezfill('muon_dz',  dz=muons.dz[mask].flatten(),  weight=w_allmu)

                # Leading muon
                w_leadmu = weight_shape(muons[leadmuon_index].pt[mask], weight[mask])
                ezfill('muon_pt0',   pt=muons[leadmuon_index].pt[mask].flatten(),    w_leadmu=w_leadmu )
                ezfill('muon_eta0',  eta=muons[leadmuon_index].eta[mask].flatten(),  w_leadmu=w_leadmu)
                ezfill('muon_phi0',  phi=muons[leadmuon_index].phi[mask].flatten(),  w_leadmu=w_leadmu)
                ezfill('muon_dxy0',  dxy=muons[leadmuon_index].dxy[mask].flatten(),  weight=w_leadmu)
                ezfill('muon_dz0',  dz=muons[leadmuon_index].dz[mask].flatten(),  weight=w_leadmu)

                # HLT Matched muons
                w_muons_hltmatch = weight_shape(muons_hltmatch.pt[mask], weight[mask])
                ezfill('muons_hltmatch_eta',  eta=muons_hltmatch.eta[mask].flatten(),  weight=w_muons_hltmatch)
                ezfill('muons_hltmatch_pt',  pt=muons_hltmatch.pt[mask].flatten(),  weight=w_muons_hltmatch)

            # Dimuon
            if '_2m_' in region:
                w_dimu = weight_shape(dimuons.pt[mask], weight[mask])

                ezfill('dimuon_pt',     pt=dimuons.pt[mask].flatten(),              weight=w_dimu)
                ezfill('dimuon_eta',    eta=dimuons.eta[mask].flatten(),            weight=w_dimu)
                ezfill('dimuon_mass',   dilepton_mass=dimuons.mass[mask].flatten(), weight=w_dimu )

            # Electrons
            if '_1e_' in region or '_2e_' in region:
                w_allel = weight_shape(electrons.pt[mask], weight[mask])
                ezfill('electron_pt',   pt=electrons.pt[mask].flatten(),    weight=w_allel)
                ezfill('electron_mt',   mt=df['MT_el'][mask],               weight=weight[mask])
                ezfill('electron_eta',  eta=electrons.eta[mask].flatten(),  weight=w_allel)
                ezfill('electron_phi',  phi=electrons.phi[mask].flatten(),  weight=w_allel)
                ezfill('electron_dz',  dz=electrons.dz[mask].flatten(),  weight=w_allel)
                ezfill('electron_dxy',  dxy=electrons.dxy[mask].flatten(),  weight=w_allel)

                w_leadel = weight_shape(electrons[leadelectron_index].pt[mask], weight[mask])
                ezfill('electron_pt0',   pt=electrons[leadelectron_index].pt[mask].flatten(),    weight=w_leadel)
                ezfill('electron_eta0',  eta=electrons[leadelectron_index].eta[mask].flatten(),  weight=w_leadel)
                ezfill('electron_phi0',  phi=electrons[leadelectron_index].phi[mask].flatten(),  weight=w_leadel)

            # Dielectron
            if '_2e_' in region:
                w_diel = weight_shape(dielectrons.pt[mask], weight[mask])
                ezfill('dielectron_pt',     pt=dielectrons.pt[mask].flatten(),                  weight=w_diel)
                ezfill('dielectron_eta',    eta=dielectrons.eta[mask].flatten(),                weight=w_diel)
                ezfill('dielectron_mass',   dilepton_mass=dielectrons.mass[mask].flatten(),     weight=w_diel)

            # Photon
            if '_g_' in region:
                w_leading_photon = weight_shape(photons[leadphoton_index].pt[mask],weight[mask]);
                ezfill('photon_pt0',              pt=photons[leadphoton_index].pt[mask].flatten(),    weight=w_leading_photon)
                ezfill('photon_pt0_noweight',     pt=photons[leadphoton_index].pt[mask].flatten(),    weight=np.ones(w_leading_photon.size))
                ezfill('photon_eta0',             eta=photons[leadphoton_index].eta[mask].flatten(),  weight=w_leading_photon)
                ezfill('photon_phi0',             phi=photons[leadphoton_index].phi[mask].flatten(),  weight=w_leading_photon)
                ezfill('photon_pt0_recoil',       pt=photons[leadphoton_index].pt[mask].flatten(), recoil=df['recoil_pt'][mask&(leadphoton_index.counts>0)],  weight=w_leading_photon)

                # w_drphoton_jet = weight_shape(df['dRPhotonJet'][mask], weight[mask])


            # PV
            ezfill('npv', nvtx=df['PV_npvs'][mask], weight=weight[mask])
            ezfill('npvgood', nvtx=df['PV_npvsGood'][mask], weight=weight[mask])

            ezfill('npv_nopu', nvtx=df['PV_npvs'][mask], weight=weight_nopu[mask])
            ezfill('npvgood_nopu', nvtx=df['PV_npvsGood'][mask], weight=weight_nopu[mask])

            ezfill('rho_all', rho=df['fixedGridRhoFastjetAll'][mask], weight=weight[mask])
            ezfill('rho_central', rho=df['fixedGridRhoFastjetCentral'][mask], weight=weight[mask])
            ezfill('rho_all_nopu', rho=df['fixedGridRhoFastjetAll'][mask], weight=weight_nopu[mask])
            ezfill('rho_central_nopu', rho=df['fixedGridRhoFastjetCentral'][mask], weight=weight_nopu[mask])
        return output

    def postprocess(self, accumulator):
        return accumulator
