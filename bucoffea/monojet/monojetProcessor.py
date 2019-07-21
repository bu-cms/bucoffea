import os
import numpy as np
import re
from coffea import hist
import coffea.processor as processor

from dynaconf import settings as cfg

from bucoffea.monojet.definitions import monojet_accumulator, monojet_evaluator, setup_candidates, setup_gen_candidates,monojet_regions
from bucoffea.helpers import min_dphi_jet_met, recoil, mt, weight_shape, bucoffea_path
from bucoffea.helpers.dataset import is_lo_z, is_lo_w, is_data, extract_year


def combine_masks(df, masks):
    """Returns the OR of the masks in the list

    :param df: Data frame
    :type df: LazyDataFrame
    :param masks: Mask names as saved in the df
    :type masks: List
    :return: OR of all masks for each event
    :rtype: array
    """
    # Start with array of False
    decision = np.ones(df.size)==0

    # Flip to true if any is passed
    for t in masks:
        try:
            decision = decision | df[t]
        except KeyError:
            continue
    return decision

class monojetProcessor(processor.ProcessorABC):
    def __init__(self, blind=True):
        self._year=None
        self._blind=blind
        self._accumulator = monojet_accumulator()

    @property
    def accumulator(self):
        return self._accumulator

    def _configure(self,df):
        dataset = df['dataset']
        self._year = extract_year(dataset)

        # Reload config based on year
        cfg.DYNACONF_WORKS="merge_configs"
        cfg.MERGE_ENABLED_FOR_DYNACONF=True
        cfg.SETTINGS_FILE_FOR_DYNACONF = bucoffea_path("config/monojet.yaml")
        cfg.ENV_FOR_DYNACONF = f"era{self._year}"
        cfg.reload()

    def process(self, df):
        if not df.size:
            return self.accumulator.identity()
        self._configure(df)
        dataset = df['dataset']
        df['is_lo_w'] = is_lo_w(dataset)
        df['is_lo_z'] = is_lo_z(dataset)
        df['is_data'] = is_data(dataset)

        if not df['is_data']:
            gen_v_pt = df['LHE_Vpt']

        # Candidates
        # Already pre-filtered!
        # All leptons are at least loose
        # Check out setup_candidates for filtering details
        ak4, ak8, muons, electrons, taus, photons, hlt = setup_candidates(df, cfg)

        # Muons
        is_tight_muon = muons.tightId \
                      & (muons.iso < cfg.MUON.CUTS.TIGHT.ISO) \
                      & (muons.pt>cfg.MUON.CUTS.TIGHT.PT) \
                      & (np.abs(muons.eta)<cfg.MUON.CUTS.TIGHT.ETA)

        dimuons = muons.distincts()
        dimuon_charge = dimuons.i0['charge'] + dimuons.i1['charge']

        df['MT_mu'] = ((muons.counts==1) * mt(muons.pt, muons.phi, df['MET_pt'], df['MET_phi'])).max()

        # Electrons
        is_tight_electron = electrons.tightId \
                            & (electrons.pt > cfg.ELECTRON.CUTS.TIGHT.PT) \
                            & (np.abs(electrons.eta) < cfg.ELECTRON.CUTS.TIGHT.ETA)

        dielectrons = electrons.distincts()
        dielectron_charge = dielectrons.i0['charge'] + dielectrons.i1['charge']

        df['MT_el'] = ((electrons.counts==1) * mt(electrons.pt, electrons.phi, df['MET_pt'], df['MET_phi'])).max()


        # ak4
        jet_acceptance = np.abs(ak4.eta)<2.4

        # B tagged ak4
        btag_cut = cfg.BTAG.CUTS[cfg.BTAG.algo][cfg.BTAG.wp]
        jet_btag_val = getattr(ak4, cfg.BTAG.algo)
        jet_btagged = jet_btag_val > btag_cut
        bjets = ak4[ jet_acceptance \
                     & jet_btagged \
                     & (ak4.pt>20) ]

        # Recoil
        df['recoil_pt'], df['recoil_phi'] = recoil(df['MET_pt'],df['MET_phi'], electrons, muons, photons)
        df["dPFCalo"] = (df['MET_pt'] - df["CaloMET_pt"]) / df["recoil_pt"]
        df["minDPhiJetRecoil"] = min_dphi_jet_met(ak4, df['recoil_phi'], njet=4, ptmin=30)
        df["minDPhiJetMet"] = min_dphi_jet_met(ak4, df['MET_phi'], njet=4, ptmin=30)
        selection = processor.PackedSelection()

        selection.add('inclusive', np.ones(df.size)==1)


        # Triggers
        if cfg.RUN.SYNC: # Synchronization mode
            pass_all = np.ones(df.size)==1
            selection.add('filt_met', pass_all)
            selection.add('trig_met', pass_all)
            selection.add('trig_ele', pass_all)
            selection.add('trig_mu',  pass_all)

        else:
            selection.add('filt_met', df['Flag_METFilters'])
            selection.add('trig_met', combine_masks(df, cfg.TRIGGERS.MET))

            # Trigger overlap
            if df['is_data']:
                if "SinglePhoton" in dataset:
                    trig_ele = combine_masks(df, cfg.TRIGGERS.ELECTRON.SINGLE_BACKUP) & (~combine_masks(df, cfg.TRIGGERS.ELECTRON.SINGLE))
                else:
                    trig_ele = combine_masks(df, cfg.TRIGGERS.ELECTRON.SINGLE)
            else:
                trig_ele = combine_masks(df, cfg.TRIGGERS.ELECTRON.SINGLE_BACKUP) | combine_masks(df, cfg.TRIGGERS.ELECTRON.SINGLE)

            selection.add('trig_ele', trig_ele)
            selection.add('trig_mu', combine_masks(df, cfg.TRIGGERS.MUON.SINGLE))
            selection.add('trig_ht_for_g_eff', combine_masks(df, cfg.TRIGGERS.HT.GAMMAEFF))

        # Trigger objects
        hlt_muons = hlt[hlt.id==13]
        hlt_single_muons = hlt_muons[hlt_muons.filter & 8 == 8]
        hlt_double_muons = hlt_muons[hlt_muons.filter & 16 == 16]

        selection.add('one_hlt_muon', hlt_single_muons.counts>=1)
        selection.add('two_hlt_muons', (hlt_single_muons.counts + 2*hlt_double_muons.counts)>=2)

        # Common selection
        selection.add('veto_ele', electrons.counts==0)
        selection.add('veto_muo', muons.counts==0)
        selection.add('veto_photon', photons.counts==0)
        selection.add('veto_tau', taus.counts==0)
        selection.add('veto_b', bjets.counts==0)
        selection.add('mindphijr',df['minDPhiJetRecoil'] > cfg.SELECTION.SIGNAL.MINDPHIJR)
        selection.add('dpfcalo',np.abs(df['dPFCalo']) < cfg.SELECTION.SIGNAL.DPFCALO)
        selection.add('recoil', df['recoil_pt']>cfg.SELECTION.SIGNAL.RECOIL)

        # AK4 Jet
        leadak4_index=ak4.pt.argmax()
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
        selection.add('mt_el', df['MT_el'] < cfg.SELECTION.CONTROL.SINGLEEL.MT)

        # Photon CR
        selection.add('trig_photon', combine_masks(df, cfg.TRIGGERS.PHOTON.SINGLE))
        leadphoton_index=photons.pt.argmax()

        is_tight_photon = photons.mediumId \
                         & (photons.pt > cfg.PHOTON.CUTS.TIGHT.PT) \
                         & (np.abs(photons.eta) < cfg.PHOTON.CUTS.TIGHT.ETA)

        selection.add('one_photon', photons.counts==1)
        selection.add('at_least_one_tight_photon', is_tight_photon.any())

        # Fill histograms
        output = self.accumulator.identity()

        # Gen
        if not df['is_data']:
            output['genvpt_check'].fill(vpt=gen_v_pt,type="Nano", dataset=dataset)

        # Weights
        evaluator = monojet_evaluator(cfg)
        all_weights = {}
        if df['is_data']:
            weight = np.ones(df.size)
        else:
            weight = df['Generator_weight']

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
            csev_sf_index = 0.5 * photons.barrel + 2.5 * ~photons.barrel
            all_weights["photon_csev"] = evaluator['photon_csev'](csev_sf_index).prod()

            all_weights["pileup"] = evaluator['pileup'](df['Pileup_nTrueInt'])

            if df['is_lo_w']:
                all_weights["theory"] = evaluator["qcd_ew_nlo_w"](gen_v_pt)
            elif df['is_lo_z']:
                all_weights["theory"] = evaluator["qcd_ew_nlo_z"](gen_v_pt)
            else:
                all_weights["theory"] = np.ones(df.size)
            for iw in all_weights.values():
                weight = weight * iw

        # Save per-event values for synchronization
        if cfg.RUN.KINEMATICS.SAVE:
            for event in cfg.RUN.KINEMATICS.EVENTS:
                mask = df['event'] == event
                if not mask.any():
                    continue
                output['kinematics']['event'] += [event]
                output['kinematics']['met'] += [df['MET_pt'][mask]]
                output['kinematics']['met_phi'] += [df['MET_phi'][mask]]
                output['kinematics']['recoil'] += [df['recoil_pt'][mask]]
                output['kinematics']['recoil_phi'] += [df['recoil_phi'][mask]]

                output['kinematics']['ak4pt0'] += [ak4[leadak4_index][mask].pt]
                output['kinematics']['ak4eta0'] += [ak4[leadak4_index][mask].eta]
                output['kinematics']['leadbtag'] += [jet_btag_val[jet_acceptance & (ak4.pt>20)][mask].max()]

                output['kinematics']['nLooseMu'] += [muons.counts[mask]]
                output['kinematics']['nTightMu'] += [muons[is_tight_muon].counts[mask]]
                output['kinematics']['mupt0'] += [muons[leadmuon_index][mask].pt]
                output['kinematics']['mueta0'] += [muons[leadmuon_index][mask].eta]

                output['kinematics']['nLooseEl'] += [electrons.counts[mask]]
                output['kinematics']['nTightEl'] += [electrons[is_tight_electron].counts[mask]]
                output['kinematics']['elpt0'] += [electrons[leadelectron_index][mask].pt]
                output['kinematics']['eleta0'] += [electrons[leadelectron_index][mask].eta]

                output['kinematics']['nLooseGam'] += [photons.counts[mask]]
                output['kinematics']['nTightGam'] += [photons[is_tight_photon].counts[mask]]
                output['kinematics']['gpt0'] += [photons[leadphoton_index][mask].pt]
                output['kinematics']['geta0'] += [photons[leadphoton_index][mask].eta]


        # Sum of all weights to use for normalization
        # TODO: Deal with systematic variations
        if not df['is_data']:
            output['sumw'][dataset] +=  df['genEventSumw']
            output['sumw2'][dataset] +=  df['genEventSumw2']

        regions = monojet_regions()
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

            # All ak4
            # This is a workaround to create a weight array of the right dimension
            w_alljets = weight_shape(ak4[mask].eta, weight[mask])


            ezfill('ak4eta',    jeteta=ak4[mask].eta.flatten(), weight=w_alljets)
            ezfill('ak4pt',     jetpt=ak4[mask].pt.flatten(),   weight=w_alljets)

            # Leading ak4
            leadak4_indices = ak4.pt.argmax()
            w_leadak4 = weight_shape(ak4[leadak4_indices].eta[mask], weight[mask])
            ezfill('ak4eta0',   jeteta=ak4[leadak4_indices].eta[mask].flatten(),    weight=w_leadak4)
            ezfill('ak4pt0',    jetpt=ak4[leadak4_indices].pt[mask].flatten(),      weight=w_leadak4)

            # All ak8
            w_allak8 = weight_shape(ak8.eta[mask], weight[mask])

            ezfill('ak8eta',    jeteta=ak8[mask].eta.flatten(), weight=w_allak8)
            ezfill('ak8pt',     jetpt=ak8[mask].pt.flatten(),   weight=w_allak8)
            ezfill('ak8mass',   mass=ak8[mask].mass.flatten(),  weight=w_allak8)

            # Leading ak8
            leadak8_indices = ak8.pt.argmax()
            w_leadak8 = weight_shape(ak8[leadak8_indices].eta[mask], weight[mask])

            ezfill('ak8eta0',   jeteta=ak8[leadak8_indices].eta[mask].flatten(),    weight=w_leadak8)
            ezfill('ak8pt0',    jetpt=ak8[leadak8_indices].pt[mask].flatten(),      weight=w_leadak8 )
            ezfill('ak8mass0',  mass=ak8[leadak8_indices].mass[mask].flatten(),     weight=w_leadak8)

            # B tag discriminator
            btag = getattr(ak4, cfg.BTAG.ALGO)
            w_btag = weight_shape(btag[mask], weight[mask])
            ezfill('ak4btag', btag=btag[mask].flatten(), weight=w_btag )

            # MET
            ezfill('dpfcalo',   dpfcalo=df["dPFCalo"][mask],    weight=weight[mask] )
            ezfill('met',       met=df["MET_pt"][mask],         weight=weight[mask] )
            ezfill('recoil',    recoil=df["recoil_pt"][mask],   weight=weight[mask] )
            ezfill('dphijm',    dphi=df["minDPhiJetMet"][mask], weight=weight[mask] )

            # Muons
            w_allmu = weight_shape(muons.pt[mask], weight[mask])
            ezfill('muon_pt',   pt=muons.pt[mask].flatten(),    weight=w_allmu )
            ezfill('muon_mt',   mt=df['MT_mu'][mask],           weight=weight[mask])
            ezfill('muon_eta',  eta=muons.eta[mask].flatten(),  weight=w_allmu)
            # Dimuon
            w_dimu = weight_shape(dimuons.pt[mask], weight[mask])

            ezfill('dimuon_pt',     pt=dimuons.pt[mask].flatten(),              weight=w_dimu)
            ezfill('dimuon_eta',    eta=dimuons.eta[mask].flatten(),            weight=w_dimu)
            ezfill('dimuon_mass',   dilepton_mass=dimuons.mass[mask].flatten(), weight=w_dimu )

            # Electrons
            w_allel = weight_shape(electrons.pt[mask], weight[mask])
            ezfill('electron_pt',   pt=electrons.pt[mask].flatten(),    weight=w_allel)
            ezfill('electron_mt',   mt=df['MT_el'][mask],               weight=weight[mask])
            ezfill('electron_eta',  eta=electrons.eta[mask].flatten(),  weight=w_allel)

            # Dielectron
            w_diel = weight_shape(dielectrons.pt[mask], weight[mask])
            ezfill('dielectron_pt',     pt=dielectrons.pt[mask].flatten(),                  weight=w_diel)
            ezfill('dielectron_eta',    eta=dielectrons.eta[mask].flatten(),                weight=w_diel)
            ezfill('dielectron_mass',   dilepton_mass=dielectrons.mass[mask].flatten(),     weight=w_diel)

            # Photon
            w_leading_photon = weight_shape(photons[leadphoton_index].pt[mask],weight[mask]);
            ezfill('photonpt0',     pt=photons[leadphoton_index].pt[mask].flatten(),    weight=w_leading_photon)
            ezfill('photoneta0',    eta=photons[leadphoton_index].eta[mask].flatten(),  weight=w_leading_photon)
            ezfill('photonphi0',    phi=photons[leadphoton_index].phi[mask].flatten(),  weight=w_leading_photon)

        return output

    def postprocess(self, accumulator):
        return accumulator



if __name__ == "__main__":
    main()