import copy
import re

import numpy as np

import coffea.processor as processor

from dynaconf import settings as cfg
from bucoffea.monojet.definitions import (
                                          monojet_accumulator,
                                          setup_candidates,
                                          monojet_regions,
                                          theory_weights_monojet,
                                          pileup_weights,
                                          candidate_weights,
                                          photon_trigger_sf,
                                          photon_impurity_weights,
                                          data_driven_qcd_dataset
                                         )
from bucoffea.helpers import (
                              min_dphi_jet_met,
                              recoil,
                              mt,
                              weight_shape,
                              bucoffea_path,
                              dphi,
                              mask_and,
                              mask_or,
                              evaluator_from_config,
                              candidates_in_hem
                             )
from bucoffea.helpers.weights import (
                              get_veto_weights,
                              diboson_nlo_weights,
                              btag_weights
                             )

from bucoffea.helpers.dataset import (
                                      is_lo_z,
                                      is_lo_znunu,
                                      is_lo_w,
                                      is_lo_g,
                                      is_nlo_z,
                                      is_nlo_w,
                                      is_nlo_g,
                                      has_v_jet,
                                      is_data,
                                      extract_year,
                                      rand_dataset_dict
                                     )
from bucoffea.helpers.gen import (
                                  setup_gen_candidates,
                                  setup_dressed_gen_candidates,
                                  fill_gen_v_info,
                                  get_gen_photon_pt
                                 )

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
        df['is_lo_g'] = is_lo_g(dataset)
        df['is_nlo_z'] = is_nlo_z(dataset)
        df['is_nlo_w'] = is_nlo_w(dataset)
        df['is_nlo_g'] = is_nlo_g(dataset)
        df['has_v_jet'] = has_v_jet(dataset)
        df['has_lhe_v_pt'] = df['is_lo_w'] | df['is_lo_z'] | df['is_nlo_z'] | df['is_nlo_w'] | df['is_lo_g']
        df['is_data'] = is_data(dataset)

        gen_v_pt = None
        if not df['is_data']:
            gen = setup_gen_candidates(df)
        if df['is_lo_w'] or df['is_lo_z'] or df['is_nlo_z'] or df['is_nlo_w']:
            dressed = setup_dressed_gen_candidates(df)
            fill_gen_v_info(df, gen, dressed)
            gen_v_pt = df['gen_v_pt_combined']
        elif df['is_lo_g'] or df['is_nlo_g']:
            gen_v_pt = get_gen_photon_pt(gen)

        # Candidates
        # Already pre-filtered!
        # All leptons are at least loose
        # Check out setup_candidates for filtering details
        met_pt, met_phi, ak4, bjets, ak8, muons, electrons, taus, photons = setup_candidates(df, cfg)
        bjets = bjets[bjets.looseId]

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

        # Photons
        # Angular distance leading photon - leading jet
        phojet_pairs = ak4[:,:1].cross(photons[:,:1])
        df['dRPhotonJet'] = np.hypot(phojet_pairs.i0.eta-phojet_pairs.i1.eta , dphi(phojet_pairs.i0.phi,phojet_pairs.i1.phi)).min()

        # Recoil
        df['recoil_pt'], df['recoil_phi'] = recoil(met_pt,met_phi, electrons, muons, photons)
        df["dPFCaloSR"] = (met_pt - df["CaloMET_pt"]) / met_pt
        df["dPFCalo"] = (met_pt - df["CaloMET_pt"]) / df["recoil_pt"]

        df["minDPhiJetRecoil"] = min_dphi_jet_met(ak4, df['recoil_phi'], njet=4, ptmin=30, etamax=2.4)
        df["minDPhiJetMet"] = min_dphi_jet_met(ak4, met_phi, njet=4, ptmin=30, etamax=2.4)
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

        # B jets are treated using veto weights
        # So accept them in MC, but reject in data
        if df['is_data']:
            selection.add('veto_b', bjets.counts==0)
        else:
            selection.add('veto_b', pass_all)

        selection.add('mindphijr',df['minDPhiJetRecoil'] > cfg.SELECTION.SIGNAL.MINDPHIJR)
        selection.add('mindphijm',df['minDPhiJetMet'] > cfg.SELECTION.SIGNAL.MINDPHIJR)
        selection.add('dpfcalo_sr',np.abs(df['dPFCaloSR']) < cfg.SELECTION.SIGNAL.DPFCALO)
        selection.add('dpfcalo',np.abs(df['dPFCalo']) < cfg.SELECTION.SIGNAL.DPFCALO)
        selection.add('recoil', df['recoil_pt']>cfg.SELECTION.SIGNAL.RECOIL)
        selection.add('met_sr', met_pt>cfg.SELECTION.SIGNAL.RECOIL)


        if df['year'] == 2018:
            selection.add('hemveto',df['hemveto'])
            selection.add('hemveto_metphi', (met_pt>470) | (met_phi>-0.62) | (met_phi<-1.62))
        else:
            selection.add('hemveto',pass_all)
            selection.add('hemveto_metphi', pass_all)
        # AK4 Jet
        leadak4_pt_eta = (ak4.pt.max() > cfg.SELECTION.SIGNAL.leadak4.PT) \
                         & (ak4.abseta[leadak4_index] < cfg.SELECTION.SIGNAL.leadak4.ETA).any()
        selection.add('leadak4_pt_eta', leadak4_pt_eta)

        selection.add('leadak4_id',(ak4.tightId[leadak4_index] \
                                                    & (ak4.chf[leadak4_index] >cfg.SELECTION.SIGNAL.leadak4.CHF) \
                                                    & (ak4.nhf[leadak4_index]<cfg.SELECTION.SIGNAL.leadak4.NHF)).any())

        # AK8 Jet
        leadak8_index=ak8.pt.argmax()
        leadak8_pt_eta = (ak8.pt.max() > cfg.SELECTION.SIGNAL.leadak8.PT) \
                         & (ak8.abseta[leadak8_index] < cfg.SELECTION.SIGNAL.leadak8.ETA).any()
        selection.add('leadak8_pt_eta', leadak8_pt_eta)

        selection.add('leadak8_id',(ak8.tightId[leadak8_index]).any())

        # Mono-V selection
        selection.add('leadak8_tau21', ((ak8.tau2[leadak8_index] / ak8.tau1[leadak8_index]) < cfg.SELECTION.SIGNAL.LEADAK8.TAU21).any())
        selection.add('leadak8_mass', ((ak8.mass[leadak8_index] > cfg.SELECTION.SIGNAL.LEADAK8.MASS.MIN) \
                                    & (ak8.mass[leadak8_index] < cfg.SELECTION.SIGNAL.LEADAK8.MASS.MAX)).any())
        selection.add('leadak8_wvsqcd_loosemd', ((ak8.wvsqcdmd[leadak8_index] > cfg.WTAG.LOOSEMD)
                                    & (ak8.wvsqcdmd[leadak8_index] < cfg.WTAG.TIGHTMD)).any())
        selection.add('leadak8_wvsqcd_tightmd', ((ak8.wvsqcdmd[leadak8_index] > cfg.WTAG.TIGHTMD)).any())
        selection.add('leadak8_wvsqcd_loose', ((ak8.wvsqcd[leadak8_index] > cfg.WTAG.LOOSE)
                                    & (ak8.wvsqcd[leadak8_index] < cfg.WTAG.TIGHT)).any())
        selection.add('leadak8_wvsqcd_medium', ((ak8.wvsqcd[leadak8_index] > cfg.WTAG.MEDIUM)
                                    ).any())
        selection.add('leadak8_wvsqcd_tight', ((ak8.wvsqcd[leadak8_index] > cfg.WTAG.TIGHT)).any())

        selection.add('veto_vtag',
            ~(
            selection.all("leadak8_pt_eta", "leadak8_id", "leadak8_wvsqcd_tight", "leadak8_mass")
            | selection.all("leadak8_pt_eta", "leadak8_id", "leadak8_wvsqcd_loose", "leadak8_mass")
            )
        )
        selection.add('only_one_ak8', ak8.counts==1)

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
        selection.add('at_least_one_tight_el', df['is_tight_electron'].any())

        selection.add('dielectron_mass', ((dielectrons.mass > cfg.SELECTION.CONTROL.DOUBLEEL.MASS.MIN)  \
                                        & (dielectrons.mass < cfg.SELECTION.CONTROL.DOUBLEEL.MASS.MAX)).any())
        selection.add('dielectron_charge', (dielectron_charge==0).any())
        selection.add('two_electrons', electrons.counts==2)

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
        if gen_v_pt is not None:
            output['genvpt_check'].fill(vpt=gen_v_pt,type="Nano", dataset=dataset, weight=df['Generator_weight'])

        if 'LHE_HT' in df:
            output['lhe_ht'].fill(dataset=dataset, ht=df['LHE_HT'])

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

            # B jet veto weights
            bsf_variations = btag_weights(bjets,cfg)
            weights.add("bveto", (1-bsf_variations["central"]).prod())

            weights = pileup_weights(weights, df, evaluator, cfg)
            if not (gen_v_pt is None):
                weights = theory_weights_monojet(weights, df, evaluator, gen_v_pt)

            # Diboson NLO
            diboson_nlo_weights(df, evaluator, gen)
            weights.add('weight_diboson_nlo', df['weight_diboson_nlo'])

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
                output['kinematics']['leadbtag'] += [ak4.pt.max()<0][mask]

                output['kinematics']['nLooseMu'] += [muons.counts[mask]]
                output['kinematics']['nTightMu'] += [muons[df['is_tight_muon']].counts[mask].flatten()]
                output['kinematics']['mupt0'] += [muons[leadmuon_index][mask].pt.flatten()]
                output['kinematics']['mueta0'] += [muons[leadmuon_index][mask].eta.flatten()]
                output['kinematics']['muphi0'] += [muons[leadmuon_index][mask].phi.flatten()]

                output['kinematics']['nLooseEl'] += [electrons.counts[mask]]
                output['kinematics']['nTightEl'] += [electrons[df['is_tight_electron']].counts[mask].flatten()]
                output['kinematics']['elpt0'] += [electrons[leadelectron_index][mask].pt.flatten()]
                output['kinematics']['eleta0'] += [electrons[leadelectron_index][mask].eta.flatten()]

                output['kinematics']['nLooseGam'] += [photons.counts[mask]]
                output['kinematics']['nTightGam'] += [photons[df['is_tight_photon']].counts[mask].flatten()]
                output['kinematics']['gpt0'] += [photons[leadphoton_index][mask].pt.flatten()]
                output['kinematics']['geta0'] += [photons[leadphoton_index][mask].eta.flatten()]


        # Randomized Parameter data sets
        # keep track of the mapping
        rand_datasets = rand_dataset_dict(df.keys(), df['year'])

        # Sum of all weights to use for normalization
        output['nevents'][dataset] += df.size
        if not df['is_data']:
            if len(rand_datasets):
                # For randomized datasets, save the normalization separately per sub-dataset
                # but also integread for the whole dataset, so that we can use both the sub
                # and total datasets for plotting
                for ds, short in rand_datasets.items():
                    dsmask = df[f'GenModel_{ds}']
                    output['nevents'][short] += dsmask.sum()
                    # Split per sub-dataset
                    output['sumw'][short] +=  getattr(df, f'genEventSumw_{ds}', 0)
                    output['sumw2'][short] +=  getattr(df,f'genEventSumw2_{ds}', 0)
                    output['sumw_pileup'][short] +=  weights.partial_weight(include=['pileup'])[dsmask].sum()

                    # Integrated for the whole dataset
                    output['sumw'][dataset] +=  getattr(df, f'genEventSumw_{ds}', 0)
                    output['sumw2'][dataset] +=  getattr(df,f'genEventSumw2_{ds}', 0)
                    output['sumw_pileup'][dataset] +=  weights.partial_weight(include=['pileup'])[dsmask].sum()
            else:
                # For normal datasets, no splitting is necessary
                output['sumw'][dataset] +=  df[f'genEventSumw']
                output['sumw2'][dataset] +=  df[f'genEventSumw2']
                output['sumw_pileup'][dataset] +=  weights.partial_weight(include=['pileup']).sum()
        regions = monojet_regions(cfg)

        # Get veto weights (only for MC)
        if not df['is_data']:
            veto_weights = get_veto_weights(df, cfg, evaluator, electrons, muons, taus, do_variations=True)

        for region, cuts in regions.items():

            if re.match('sr_.*', region):
                recoil_pt = met_pt
                recoil_phi = met_phi
            else:
                recoil_pt = df['recoil_pt']
                recoil_phi = df['recoil_phi']

            exclude = [None]
            region_weights = copy.deepcopy(weights)

            if not df['is_data']:
                ### Trigger weights
                if re.match(r'cr_(\d+)e.*', region):
                    p_pass_data = 1 - (1-evaluator["trigger_electron_eff_data"](electrons.etasc, electrons.pt)).prod()
                    p_pass_mc   = 1 - (1-evaluator["trigger_electron_eff_mc"](electrons.etasc, electrons.pt)).prod()
                    trigger_weight = p_pass_data/p_pass_mc
                    trigger_weight[np.isnan(trigger_weight)] = 1
                    region_weights.add('trigger_ele', trigger_weight)
                elif re.match(r'cr_(\d+)m.*', region) or re.match('sr_.*', region):
                    region_weights.add('trigger_met', evaluator["trigger_met"](recoil_pt))
                elif re.match(r'cr_g.*', region):
                    photon_trigger_sf(region_weights, photons, df)

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
                    region_weights.add("vetoweight", veto_weights.partial_weight(include=["nominal"]))

            if not (df['is_data']):
                genVs = gen[((gen.pdg==23) | (gen.pdg==24) | (gen.pdg==-24)) & (gen.pt>10)]
                leadak8 = ak8[ak8.pt.argmax()]
                leadak8_matched_mask = leadak8.match(genVs, deltaRCut=0.8)
                matched_leadak8 = leadak8[leadak8_matched_mask]
                unmatched_leadak8 = leadak8[~leadak8_matched_mask]
                for wp in ['loose','tight','medium']:
                    if re.match(f'.*_{wp}_v.*', region):
                        if ('nomistag' in region) or wp=='medium': 
                            matched_weights = evaluator[f'wtag_{wp}'](matched_leadak8.pt).prod()
                        #elif re.match(r'cr_g.*', region):
                        #    matched_weights = evaluator[f'wtag_{wp}'](matched_leadak8.pt).prod() \
                        #            * evaluator[f'wtag_mistag_g_{wp}'](unmatched_leadak8.pt).prod()
                        else:
                            matched_weights = evaluator[f'wtag_{wp}'](matched_leadak8.pt).prod() \
                                    * evaluator[f'wtag_mistag_g_{wp}'](unmatched_leadak8.pt).prod()

                        region_weights.add('wtag_{wp}', matched_weights)



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
                if re.match(cfg.RUN.SAVE.TREEREGIONS, region):
                    # General properties
                    output['tree_int64'][region]["event"]                   += processor.column_accumulator(df["event"][mask])
                    output['tree_int64'][region]["run"]                     += processor.column_accumulator(df["run"][mask])
                    output['tree_int64'][region]["lumi"]                    += processor.column_accumulator(df["luminosityBlock"][mask])



                    # Selection bits
                    if region=='inclusive':
                        for name in selection.names:
                            output['tree_bool'][region][name] += processor.column_accumulator(np.bool_(selection.all(*[name])[mask]))
                    else:
                        output['tree_float16'][region]["recoil_pt"]             += processor.column_accumulator(recoil_pt[mask])
                        output['tree_float16'][region]["recoil_phi"]            += processor.column_accumulator(recoil_phi[mask])
                        output['tree_float16'][region]["met_pt"]                += processor.column_accumulator(met_pt[mask])
                        output['tree_float16'][region]["met_phi"]               += processor.column_accumulator(met_phi[mask])
                        output['tree_float16'][region]["met_pt_nojer"]          += processor.column_accumulator(df['MET_pt_nom' if df['year']==2018 else 'METFixEE2017_pt_nom'][mask])
                        output['tree_float16'][region]["met_phi_nojer"]         += processor.column_accumulator(df['MET_phi_nom' if df['year']==2018 else 'METFixEE2017_phi_nom'][mask])
                        output['tree_float16'][region]["leadak4_pt"]            += processor.column_accumulator(ak4[leadak4_index].pt.max()[mask])
                        output['tree_float16'][region]["leadak4_eta"]           += processor.column_accumulator(ak4[leadak4_index].eta.max()[mask])
                        output['tree_float16'][region]["leadak4_phi"]           += processor.column_accumulator(ak4[leadak4_index].phi.max()[mask])

                        output['tree_float16'][region]["mindphijr"]            += processor.column_accumulator(df['minDPhiJetRecoil'][mask])

                        output['tree_float16'][region]["calomet_pt"]            += processor.column_accumulator(df['CaloMET_pt'][mask])


                        # MC quantities
                        if not df['is_data']:
                            if gen_v_pt is not None:
                                output['tree_float16'][region]["gen_v_pt"]              += processor.column_accumulator(gen_v_pt[mask] if gen_v_pt is not None else np.zeros(len(df["event"][mask])))

                            for name, w in region_weights._weights.items():
                                output['tree_float16'][region][f"weight_{name}"] += processor.column_accumulator(np.float16(w[mask]))


                        # Transverse mass single muon
                        if re.match('.*_1m_.*', region):
                            output['tree_float16'][region]["mt"]  += processor.column_accumulator(df['MT_mu'][mask])

                        # Transverse mass single electron
                        if re.match('.*_1e_.*', region):
                            output['tree_float16'][region]["mt"]  += processor.column_accumulator(df['MT_el'][mask])

                        # Photon
                        if re.match('.*_g_.*', region):
                            output['tree_float16'][region]["gam0pt"]  += processor.column_accumulator(photons.pt[leadphoton_index][mask].max())
                            output['tree_float16'][region]["gam0eta"] += processor.column_accumulator(photons.eta[leadphoton_index][mask].max())

                        # Leading muon
                        if re.match('.*_(\d)m_.*', region):
                            output['tree_float16'][region]["mu0pt"]  += processor.column_accumulator(muons.pt[leadmuon_index][mask].max())
                            output['tree_float16'][region]["mu0eta"] += processor.column_accumulator(muons.eta[leadmuon_index][mask].max())
                            output['tree_float16'][region]["mu0tight"] += processor.column_accumulator(muons.tightId[leadmuon_index][mask].max())

                        # Trailing muon
                        if re.match('.*_2m_.*', region):
                            output['tree_float16'][region]["mu1pt"]  += processor.column_accumulator(muons.pt[~leadmuon_index][mask].max())
                            output['tree_float16'][region]["mu1eta"] += processor.column_accumulator(muons.eta[~leadmuon_index][mask].max())
                            output['tree_float16'][region]["mu1tight"] += processor.column_accumulator(muons.tightId[~leadmuon_index][mask].max())

                        # Leading electron
                        if re.match('.*_(\d)e_.*', region):
                            output['tree_float16'][region]["el0pt"]  += processor.column_accumulator(electrons.pt[leadelectron_index][mask].max())
                            output['tree_float16'][region]["el0eta"] += processor.column_accumulator(electrons.eta[leadelectron_index][mask].max())
                            output['tree_float16'][region]["el0tight"] += processor.column_accumulator(electrons.tightId[leadelectron_index][mask].max())

                        # Trailing electron
                        if re.match('.*_2e_.*', region):
                            output['tree_float16'][region]["el1pt"]  += processor.column_accumulator(electrons.pt[~leadelectron_index][mask].max())
                            output['tree_float16'][region]["el1eta"] += processor.column_accumulator(electrons.eta[~leadelectron_index][mask].max())
                            output['tree_float16'][region]["el1tight"] += processor.column_accumulator(electrons.tightId[~leadelectron_index][mask].max())



            if region=='inclusive':
                continue
            # Save the event numbers of events passing this selection
            if cfg.RUN.SAVE.PASSING:
                # Save only every Nth event
                save_mask = mask & ((df['event']%cfg.RUN.SAVE.PRESCALE)== 0)
                output['selected_events'][region] += processor.column_accumulator(df['event'][save_mask].astype(np.uint64))


            # Multiplicities
            def fill_mult(name, candidates):
                output[name].fill(
                                  dataset=dataset,
                                  region=region,
                                  multiplicity=candidates[mask].counts,
                                  weight=region_weights.partial_weight(exclude=exclude)[mask]
                                  )

            fill_mult('ak8_mult', ak8)
            fill_mult('ak4_mult', ak4[(ak4.pt>30)&(ak4.abseta<2.4)])
            fill_mult('bjet_mult',bjets)
            fill_mult('loose_ele_mult',electrons)
            fill_mult('tight_ele_mult',electrons[df['is_tight_electron']])
            fill_mult('loose_muo_mult',muons)
            fill_mult('tight_muo_mult',muons[df['is_tight_muon']])
            fill_mult('tau_mult',taus)
            fill_mult('photon_mult',photons)

            def ezfill(name, **kwargs):
                """Helper function to make filling easier."""
                if not ('dataset' in kwargs):
                    kwargs['dataset'] = dataset
                output[name].fill(
                                  region=region,
                                  **kwargs
                                  )
            # Monitor weights
            for wname, wvalue in region_weights._weights.items():
                ezfill("weights", weight_type=wname, weight_value=wvalue[mask])

            # All ak4
            # This is a workaround to create a weight array of the right dimension
            w_alljets = weight_shape(ak4[mask].eta, region_weights.partial_weight(exclude=exclude)[mask])

            ezfill('ak4_eta',    jeteta=ak4[mask].eta.flatten(), weight=w_alljets)
            ezfill('ak4_phi',    jetphi=ak4[mask].phi.flatten(), weight=w_alljets)
            ezfill('ak4_eta_phi', phi=ak4[mask].phi.flatten(),eta=ak4[mask].eta.flatten(), weight=w_alljets)
            ezfill('ak4_pt',     jetpt=ak4[mask].pt.flatten(),   weight=w_alljets)
            ezfill('ak4_deepcsv', deepcsv=ak4[mask].deepcsv.flatten(),   weight=w_alljets)

            w_bjets = weight_shape(bjets[mask].eta, region_weights.partial_weight(exclude=["bveto"])[mask])
            ezfill('bjet_eta',    jeteta=bjets[mask].eta.flatten(), weight=w_bjets)
            ezfill('bjet_phi',    jetphi=bjets[mask].phi.flatten(), weight=w_bjets)
            ezfill('bjet_pt',     jetpt=bjets[mask].pt.flatten(),   weight=w_bjets)

            # Leading ak4
            w_leadak4 = weight_shape(ak4[leadak4_index].eta[mask], region_weights.partial_weight(exclude=exclude)[mask])
            ezfill('ak4_eta0',   jeteta=ak4[leadak4_index].eta[mask].flatten(),    weight=w_leadak4)
            ezfill('ak4_phi0',   jetphi=ak4[leadak4_index].phi[mask].flatten(),    weight=w_leadak4)
            ezfill('ak4_pt0',    jetpt=ak4[leadak4_index].pt[mask].flatten(),      weight=w_leadak4)
            ezfill('ak4_ptraw0',    jetpt=ak4[leadak4_index].ptraw[mask].flatten(),      weight=w_leadak4)
            ezfill('ak4_chf0',    frac=ak4[leadak4_index].chf[mask].flatten(),      weight=w_leadak4)
            ezfill('ak4_nhf0',    frac=ak4[leadak4_index].nhf[mask].flatten(),      weight=w_leadak4)
            ezfill('ak4_nef0',    frac=ak4[leadak4_index].nef[mask].flatten(),      weight=w_leadak4)

            rw=region_weights.partial_weight(exclude=exclude)
            ezfill('drelejet',    dr=df['dREleJet'][mask],      weight=rw[mask])
            ezfill('drmuonjet',    dr=df['dRMuonJet'][mask],      weight=rw[mask])
            ezfill('drphotonjet',    dr=df['dRPhotonJet'][mask],  weight=rw[mask])

            # AK8 jets
            if region=='inclusive' or region.endswith('v'):
                # All
                w_allak8 = weight_shape(ak8.eta[mask], region_weights.partial_weight(exclude=exclude)[mask])

                ezfill('ak8_eta',    jeteta=ak8[mask].eta.flatten(), weight=w_allak8)
                ezfill('ak8_phi',    jetphi=ak8[mask].phi.flatten(), weight=w_allak8)
                ezfill('ak8_pt',     jetpt=ak8[mask].pt.flatten(),   weight=w_allak8)
                ezfill('ak8_mass',   mass=ak8[mask].mass.flatten(),  weight=w_allak8)

                # Leading
                w_leadak8 = weight_shape(ak8[leadak8_index].eta[mask], region_weights.partial_weight(exclude=exclude)[mask])

                ezfill('ak8_eta0',       jeteta=ak8[leadak8_index].eta[mask].flatten(),    weight=w_leadak8)
                ezfill('ak8_phi0',       jetphi=ak8[leadak8_index].phi[mask].flatten(),    weight=w_leadak8)
                ezfill('ak8_pt0',        jetpt=ak8[leadak8_index].pt[mask].flatten(),      weight=w_leadak8 )
                ezfill('ak8_mass0',      mass=ak8[leadak8_index].mass[mask].flatten(),     weight=w_leadak8)
                ezfill('ak8_tau210',     tau21=ak8[leadak8_index].tau21[mask].flatten(),     weight=w_leadak8)
                ezfill('ak8_wvsqcd0',    tagger=ak8[leadak8_index].wvsqcd[mask].flatten(),     weight=w_leadak8)
                ezfill('ak8_wvsqcdmd0',  tagger=ak8[leadak8_index].wvsqcdmd[mask].flatten(),     weight=w_leadak8)
                ezfill('ak8_zvsqcd0',    tagger=ak8[leadak8_index].zvsqcd[mask].flatten(),     weight=w_leadak8)
                ezfill('ak8_zvsqcdmd0',  tagger=ak8[leadak8_index].zvsqcdmd[mask].flatten(),     weight=w_leadak8)
                ezfill('ak8_tvsqcd0',    tagger=ak8[leadak8_index].tvsqcd[mask].flatten(),     weight=w_leadak8)
                ezfill('ak8_tvsqcdmd0',    tagger=ak8[leadak8_index].tvsqcdmd[mask].flatten(),     weight=w_leadak8)
                ezfill('ak8_wvstqcd0',    tagger=ak8[leadak8_index].wvstqcd[mask].flatten(),     weight=w_leadak8)
                ezfill('ak8_wvstqcdmd0',    tagger=ak8[leadak8_index].wvstqcdmd[mask].flatten(),     weight=w_leadak8)

                # histogram with only gen-matched lead ak8 pt
                if not df['is_data']:
                    w_matchedleadak8 = weight_shape(matched_leadak8.eta[mask], region_weights.partial_weight(exclude=exclude)[mask])
                    ezfill('ak8_Vmatched_pt0', jetpt=matched_leadak8.pt[mask].flatten(),      weight=w_matchedleadak8 )


                # specifically for deepak8 mistag rate measurement
                if cfg.RUN.MONOVMISTAG and 'inclusive_v' in region:
                    ezfill('ak8_passloose_pt0', wppass=ak8[leadak8_index].wvsqcd[mask].max()>cfg.WTAG.LOOSE, jetpt=ak8[leadak8_index].pt[mask].max(),      weight=w_leadak8 )
                    ezfill('ak8_passmedium_pt0', wppass=ak8[leadak8_index].wvsqcd[mask].max()>cfg.WTAG.LOOSE, jetpt=ak8[leadak8_index].pt[mask].max(),      weight=w_leadak8 )
                    ezfill('ak8_passtight_pt0', wppass=ak8[leadak8_index].wvsqcd[mask].max()>cfg.WTAG.TIGHT, jetpt=ak8[leadak8_index].pt[mask].max(),      weight=w_leadak8 )
                    ezfill('ak8_passloosemd_pt0', wppass=ak8[leadak8_index].wvsqcdmd[mask].max()>cfg.WTAG.LOOSEMD, jetpt=ak8[leadak8_index].pt[mask].max(),      weight=w_leadak8 )
                    ezfill('ak8_passtightmd_pt0', wppass=ak8[leadak8_index].wvsqcdmd[mask].max()>cfg.WTAG.TIGHTMD, jetpt=ak8[leadak8_index].pt[mask].max(),      weight=w_leadak8 )
                    ezfill('ak8_passloose_mass0', wppass=ak8[leadak8_index].wvsqcd[mask].max()>cfg.WTAG.LOOSE, mass=ak8[leadak8_index].mass[mask].max(),      weight=w_leadak8 )
                    ezfill('ak8_passmedium_mass0', wppass=ak8[leadak8_index].wvsqcd[mask].max()>cfg.WTAG.LOOSE, mass=ak8[leadak8_index].mass[mask].max(),      weight=w_leadak8 )
                    ezfill('ak8_passtight_mass0', wppass=ak8[leadak8_index].wvsqcd[mask].max()>cfg.WTAG.TIGHT, mass=ak8[leadak8_index].mass[mask].max(),      weight=w_leadak8 )
                    ezfill('ak8_passloosemd_mass0', wppass=ak8[leadak8_index].wvsqcdmd[mask].max()>cfg.WTAG.LOOSEMD, mass=ak8[leadak8_index].mass[mask].max(),      weight=w_leadak8 )
                    ezfill('ak8_passtightmd_mass0', wppass=ak8[leadak8_index].wvsqcdmd[mask].max()>cfg.WTAG.TIGHTMD, mass=ak8[leadak8_index].mass[mask].max(),      weight=w_leadak8 )

            # MET
            rw = region_weights.partial_weight(exclude=exclude)
            ezfill('dpfcalo',            dpfcalo=df["dPFCalo"][mask], weight=rw[mask])
            ezfill('met',                met=met_pt[mask],            weight=rw[mask] )
            ezfill('met_phi',            phi=met_phi[mask],           weight=rw[mask] )
            ezfill('recoil',             recoil=recoil_pt[mask],      weight=rw[mask] )
            ezfill('recoil_phi',         phi=recoil_phi[mask],        weight=rw[mask] )
            ezfill('recoil_nopog',       recoil=recoil_pt[mask],      weight=region_weights.partial_weight(include=['pileup','theory','gen','prefire'])[mask])
            ezfill('recoil_nopref',      recoil=recoil_pt[mask],      weight=region_weights.partial_weight(exclude=['prefire']+exclude)[mask])
            ezfill('recoil_nopu',        recoil=recoil_pt[mask],      weight=region_weights.partial_weight(exclude=['pileup']+exclude)[mask])
            ezfill('recoil_notrg',       recoil=recoil_pt[mask],      weight=region_weights.partial_weight(exclude=['trigger']+exclude)[mask])
            ezfill('ak4_pt0_over_recoil',    ratio=ak4.pt.max()[mask]/recoil_pt[mask],      weight=region_weights.partial_weight(exclude=exclude)[mask])
            ezfill('dphijm',             dphi=df["minDPhiJetMet"][mask],    weight=region_weights.partial_weight(exclude=exclude)[mask] )
            ezfill('dphijr',             dphi=df["minDPhiJetRecoil"][mask],    weight=region_weights.partial_weight(exclude=exclude)[mask] )

            # Diboson NLO
            ezfill(
                    'recoil_nodibosonnlo',
                    recoil=recoil_pt[mask],
                    weight=region_weights.partial_weight(exclude=['weight_diboson_nlo']+exclude)[mask]
                    )

            if not df['is_data']:
                ezfill(
                        'recoil_dibosonnlo_up',
                        recoil=recoil_pt[mask],
                        weight=(region_weights.partial_weight(exclude=exclude)*(1+df['weight_diboson_nlo_rel_unc']))[mask]
                        )
                ezfill(
                        'recoil_dibosonnlo_dn',
                        recoil=recoil_pt[mask],
                        weight=(region_weights.partial_weight(exclude=exclude)*(1-df['weight_diboson_nlo_rel_unc']))[mask]
                        )

            # Randomized parameter samples
            for ds, short in rand_datasets.items():
                dsmask = df[f'GenModel_{ds}']
                ezfill('recoil', recoil=recoil_pt[mask&dsmask],      weight=rw[mask&dsmask], dataset=short )

            if cfg.RUN.BTAG_STUDY:
                ezfill('recoil_hardbveto',   recoil=recoil_pt[mask&(bjets.counts==0)],      weight=region_weights.partial_weight(exclude=exclude+['bveto'])[mask&(bjets.counts==0)])
                if not df['is_data']:
                    rw = region_weights.partial_weight(exclude=exclude+['bveto'])
                    ezfill('recoil_bveto_up',    recoil=recoil_pt[mask],  weight=(rw*(1-bsf_variations['up']).prod())[mask])
                    ezfill('recoil_bveto_down',  recoil=recoil_pt[mask],  weight=(rw*(1-bsf_variations['down']).prod())[mask])

            if cfg.RUN.PHOTON_ID_STUDY and (not df['is_data']) and ('cr_g' in region) and (df['year']!=2016):
                photon_id_sf_nom = evaluator['photon_id_tight_tnp'](np.abs(photons[df['is_tight_photon']].eta))
                photon_id_sf_err = evaluator['photon_id_tight_tnp_error'](np.abs(photons[df['is_tight_photon']].eta))

                rw = region_weights.partial_weight(exclude=exclude+['photon_id_tight'])
                ezfill('recoil_photon_id_up', recoil=recoil_pt[mask], weight=(rw * (photon_id_sf_nom+photon_id_sf_err).prod())[mask])
                ezfill('recoil_photon_id_dn', recoil=recoil_pt[mask], weight=(rw * (photon_id_sf_nom-photon_id_sf_err).prod())[mask])

                photon_id_sf_err_extrap_err = evaluator['photon_id_tight_tnp_extrap_unc_slope'](np.abs(photons[df['is_tight_photon']].eta)) * (photons[df['is_tight_photon']].pt - 150)

                ezfill('recoil_photon_id_extrap_up', recoil=recoil_pt[mask], weight=(rw * (photon_id_sf_nom+photon_id_sf_err_extrap_err)).prod()[mask])
                ezfill('recoil_photon_id_extrap_dn', recoil=recoil_pt[mask], weight=(rw * (photon_id_sf_nom-photon_id_sf_err_extrap_err)).prod()[mask])

            if cfg.RUN.ELE_ID_STUDY and (not df['is_data']) and ('cr_1e' in region or 'cr_2e' in region) and (df['year']!=2016):
                # note that electrons in the gap do not count in this study, the "nominal recoil" distribution is different from the default "recoil" distribution
                mask_electron_nogap = (np.abs(electrons.etasc)<1.4442) | (np.abs(electrons.etasc)>1.566)
                electrons_nogap = electrons[mask_electron_nogap]
                electron_is_tight_electron = df['is_tight_electron'][mask_electron_nogap]
                electrons_nogap_tight = electrons_nogap[ electron_is_tight_electron]
                electrons_nogap_loose = electrons_nogap[~electron_is_tight_electron]
                eletight_id_sf = {
                        "up": evaluator['ele_id_tight'](electrons_nogap_tight.etasc, electrons_nogap_tight.pt) + evaluator['ele_id_tight_error'](electrons_nogap_tight.etasc, electrons_nogap_tight.pt),
                        "dn": evaluator['ele_id_tight'](electrons_nogap_tight.etasc, electrons_nogap_tight.pt) - evaluator['ele_id_tight_error'](electrons_nogap_tight.etasc, electrons_nogap_tight.pt),
                        "nm": evaluator['ele_id_tight'](electrons_nogap_tight.etasc, electrons_nogap_tight.pt)}
                eleloose_id_sf = {
                        "up": evaluator['ele_id_loose'](electrons_nogap_loose.etasc, electrons_nogap_loose.pt) + evaluator['ele_id_loose_error'](electrons_nogap_loose.etasc, electrons_nogap_loose.pt),
                        "dn": evaluator['ele_id_loose'](electrons_nogap_loose.etasc, electrons_nogap_loose.pt) - evaluator['ele_id_loose_error'](electrons_nogap_loose.etasc, electrons_nogap_loose.pt),
                        "nm": evaluator['ele_id_loose'](electrons_nogap_loose.etasc, electrons_nogap_loose.pt)}
                if cfg.SF.DIELE_ID_SF.USE_AVERAGE:
                    tight_dielectrons = electrons_nogap[electron_is_tight_electron].distincts()
                    eletight0_sf = evaluator['ele_id_tight'       ](tight_dielectrons.i0.etasc, tight_dielectrons.i0.pt).prod()
                    eletight0_er = evaluator['ele_id_tight_error' ](tight_dielectrons.i0.etasc, tight_dielectrons.i0.pt).prod()
                    eleloose0_sf = evaluator['ele_id_loose'       ](tight_dielectrons.i0.etasc, tight_dielectrons.i0.pt).prod()
                    eleloose0_er = evaluator['ele_id_loose_error' ](tight_dielectrons.i0.etasc, tight_dielectrons.i0.pt).prod()
                    eletight1_sf = evaluator['ele_id_tight'       ](tight_dielectrons.i1.etasc, tight_dielectrons.i1.pt).prod()
                    eletight1_er = evaluator['ele_id_tight_error' ](tight_dielectrons.i1.etasc, tight_dielectrons.i1.pt).prod()
                    eleloose1_sf = evaluator['ele_id_loose'       ](tight_dielectrons.i1.etasc, tight_dielectrons.i1.pt).prod()
                    eleloose1_er = evaluator['ele_id_loose_error' ](tight_dielectrons.i1.etasc, tight_dielectrons.i1.pt).prod()
                    weights_2e_tight_up = 0.5*((eletight0_sf+eletight0_er)*(eleloose1_sf+eleloose1_er) + (eletight1_sf+eletight1_er)*(eleloose0_sf+eleloose0_er))
                    weights_2e_tight_dn = 0.5*((eletight0_sf-eletight0_er)*(eleloose1_sf-eleloose1_er) + (eletight1_sf-eletight1_er)*(eleloose0_sf-eleloose0_er))
                    weights_2e_tight_nm = 0.5*((eletight0_sf)*(eleloose1_sf) + (eletight1_sf)*(eleloose0_sf))

                    eletight_id_sf["up"] = eletight_id_sf["up"]*(tight_dielectrons.counts != 1) + weights_2e_tight_up*(tight_dielectrons.counts == 1)
                    eletight_id_sf["dn"] = eletight_id_sf["dn"]*(tight_dielectrons.counts != 1) + weights_2e_tight_dn*(tight_dielectrons.counts == 1)
                    eletight_id_sf["nm"] = eletight_id_sf["nm"]*(tight_dielectrons.counts != 1) + weights_2e_tight_nm*(tight_dielectrons.counts == 1)

                rw = region_weights.partial_weight(exclude=exclude+['ele_id_tight','ele_id_loose'])
                ezfill('recoil_ele_id_up', recoil=recoil_pt[mask], weight=(rw * eletight_id_sf["up"].prod() * eleloose_id_sf["up"].prod())[mask])
                ezfill('recoil_ele_id_dn', recoil=recoil_pt[mask], weight=(rw * eletight_id_sf["dn"].prod() * eleloose_id_sf["dn"].prod())[mask])
                ezfill('recoil_ele_id_nm', recoil=recoil_pt[mask], weight=(rw * eletight_id_sf["nm"].prod() * eleloose_id_sf["nm"].prod())[mask])

                ### Electron Reco efficiency up&down variations
                ele_reco_sf = {}
                if df['year']==2017:
                    high_et = electrons.pt>20
                    ele_reco_sf["up"] = (evaluator['ele_reco'](electrons.etasc[high_et], electrons.pt[high_et]) + evaluator['ele_reco_error'](electrons.etasc[high_et], electrons.pt[high_et])).prod()
                    ele_reco_sf["dn"] = (evaluator['ele_reco'](electrons.etasc[high_et], electrons.pt[high_et]) - evaluator['ele_reco_error'](electrons.etasc[high_et], electrons.pt[high_et])).prod()
                    ele_reco_sf["up"] *= (evaluator['ele_reco_pt_lt_20'](electrons.etasc[~high_et], electrons.pt[~high_et]) + evaluator['ele_reco_pt_lt_20_error'](electrons.etasc[~high_et], electrons.pt[~high_et])).prod()
                    ele_reco_sf["dn"] *= (evaluator['ele_reco_pt_lt_20'](electrons.etasc[~high_et], electrons.pt[~high_et]) - evaluator['ele_reco_pt_lt_20_error'](electrons.etasc[~high_et], electrons.pt[~high_et])).prod()
                else:
                    ele_reco_sf["up"] = (evaluator['ele_reco'](electrons.etasc, electrons.pt) + evaluator['ele_reco_error'](electrons.etasc, electrons.pt)).prod()
                    ele_reco_sf["dn"] = (evaluator['ele_reco'](electrons.etasc, electrons.pt) - evaluator['ele_reco_error'](electrons.etasc, electrons.pt)).prod()
                rw = region_weights.partial_weight(exclude=exclude+['ele_reco'])
                ezfill('recoil_ele_reco_up', recoil=recoil_pt[mask], weight=(rw * ele_reco_sf["up"])[mask])
                ezfill('recoil_ele_reco_dn', recoil=recoil_pt[mask], weight=(rw * ele_reco_sf["dn"])[mask])

            if re.match('.*no_veto.*', region) and not df['is_data']:
                for variation in veto_weights._weights.keys():
                    ezfill(
                            "recoil_veto_weight",
                            recoil=recoil_pt[mask],
                            weight=region_weights.partial_weight(exclude=exclude+["vetoweight"])[mask]*veto_weights.partial_weight(include=[variation])[mask],
                            variation=variation
                            )

            # Photon CR data-driven QCD estimate
            if df['is_data'] and re.match("cr_g.*", region) and re.match("(SinglePhoton|EGamma).*", dataset):
                w_imp = photon_impurity_weights(photons[leadphoton_index].pt.max()[mask], df["year"])
                output['recoil'].fill(
                                    dataset=data_driven_qcd_dataset(dataset),
                                    region=region,
                                    recoil=recoil_pt[mask],
                                    weight=region_weights.partial_weight(exclude=exclude)[mask] * w_imp
                                )

            if 'noveto' in region:
                continue

            # Muons
            if '_1m_' in region or '_2m_' in region:
                w_allmu = weight_shape(muons.pt[mask], region_weights.partial_weight(exclude=exclude)[mask])
                ezfill('muon_pt',   pt=muons.pt[mask].flatten(),    weight=w_allmu )
                ezfill('muon_mt',   mt=df['MT_mu'][mask],           weight=region_weights.partial_weight(exclude=exclude)[mask])
                ezfill('muon_eta',  eta=muons.eta[mask].flatten(),  weight=w_allmu)
                ezfill('muon_eta_phi', phi=muons.phi[mask].flatten(),eta=muons.eta[mask].flatten(), weight=w_allmu)
                ezfill('muon_phi',  phi=muons.phi[mask].flatten(),  weight=w_allmu)
                ezfill('muon_dxy',  dxy=muons.dxy[mask].flatten(),  weight=w_allmu)
                ezfill('muon_dz',  dz=muons.dz[mask].flatten(),  weight=w_allmu)

                # Leading muon
                w_leadmu = weight_shape(muons[leadmuon_index].pt[mask], region_weights.partial_weight(exclude=exclude)[mask])
                ezfill('muon_pt0',   pt=muons[leadmuon_index].pt[mask].flatten(),    weight=w_leadmu )
                ezfill('muon_eta0',  eta=muons[leadmuon_index].eta[mask].flatten(),  weight=w_leadmu)
                ezfill('muon_phi0',  phi=muons[leadmuon_index].phi[mask].flatten(),  weight=w_leadmu)
                ezfill('muon_dxy0',  dxy=muons[leadmuon_index].dxy[mask].flatten(),  weight=w_leadmu)
                ezfill('muon_dz0',  dz=muons[leadmuon_index].dz[mask].flatten(),  weight=w_leadmu)

            # Dimuon
            if '_2m_' in region:
                w_dimu = weight_shape(dimuons.pt[mask], region_weights.partial_weight(exclude=exclude)[mask])

                ezfill('dimuon_pt',     pt=dimuons.pt[mask].flatten(),              weight=w_dimu)
                ezfill('dimuon_eta',    eta=dimuons.eta[mask].flatten(),            weight=w_dimu)
                ezfill('dimuon_mass',   dilepton_mass=dimuons.mass[mask].flatten(), weight=w_dimu )
                ezfill('dimuon_dr',   dr=dimuons.i0.p4.delta_r(dimuons.i1.p4)[mask].flatten(), weight=w_dimu )

                ezfill('muon_pt1',   pt=muons[~leadmuon_index].pt[mask].flatten(),    weight=w_leadmu )
                ezfill('muon_eta1',  eta=muons[~leadmuon_index].eta[mask].flatten(),  weight=w_leadmu)
                ezfill('muon_phi1',  phi=muons[~leadmuon_index].phi[mask].flatten(),  weight=w_leadmu)

            # Electrons
            if '_1e_' in region or '_2e_' in region:
                w_allel = weight_shape(electrons.pt[mask], region_weights.partial_weight(exclude=exclude)[mask])
                ezfill('electron_pt',   pt=electrons.pt[mask].flatten(),    weight=w_allel)
                ezfill('electron_mt',   mt=df['MT_el'][mask],               weight=region_weights.partial_weight(exclude=exclude)[mask])
                ezfill('electron_eta',  eta=electrons.eta[mask].flatten(),  weight=w_allel)
                ezfill('electron_phi',  phi=electrons.phi[mask].flatten(),  weight=w_allel)
                ezfill('electron_eta_phi', phi=electrons.phi[mask].flatten(),eta=electrons.eta[mask].flatten(), weight=w_allel)
                ezfill('electron_dz',  dz=electrons.dz[mask].flatten(),  weight=w_allel)
                ezfill('electron_dxy',  dxy=electrons.dxy[mask].flatten(),  weight=w_allel)

                w_leadel = weight_shape(electrons[leadelectron_index].pt[mask], region_weights.partial_weight(exclude=exclude)[mask])
                ezfill('electron_pt0',   pt=electrons[leadelectron_index].pt[mask].flatten(),    weight=w_leadel)
                ezfill('electron_eta0',  eta=electrons[leadelectron_index].eta[mask].flatten(),  weight=w_leadel)
                ezfill('electron_phi0',  phi=electrons[leadelectron_index].phi[mask].flatten(),  weight=w_leadel)

                w_trailel = weight_shape(electrons[~leadelectron_index].pt[mask], region_weights.partial_weight(exclude=exclude)[mask])
                ezfill('electron_tightid1',  id=electrons[~leadelectron_index].tightId[mask].flatten(),  weight=w_trailel)

            # Dielectron
            if '_2e_' in region:
                w_diel = weight_shape(dielectrons.pt[mask], region_weights.partial_weight(exclude=exclude)[mask])
                ezfill('dielectron_pt',     pt=dielectrons.pt[mask].flatten(),                  weight=w_diel)
                ezfill('dielectron_eta',    eta=dielectrons.eta[mask].flatten(),                weight=w_diel)
                ezfill('dielectron_mass',   dilepton_mass=dielectrons.mass[mask].flatten(),     weight=w_diel)
                ezfill('dielectron_dr',   dr=dielectrons.i0.p4.delta_r(dielectrons.i1.p4)[mask].flatten(), weight=w_diel )

                ezfill('electron_pt1',   pt=electrons[~leadelectron_index].pt[mask].flatten(),    weight=w_leadel)
                ezfill('electron_eta1',  eta=electrons[~leadelectron_index].eta[mask].flatten(),  weight=w_leadel)
                ezfill('electron_phi1',  phi=electrons[~leadelectron_index].phi[mask].flatten(),  weight=w_leadel)
            # Photon
            if '_g_' in region:
                w_leading_photon = weight_shape(photons[leadphoton_index].pt[mask],region_weights.partial_weight(exclude=exclude)[mask]);
                ezfill('photon_pt0',              pt=photons[leadphoton_index].pt[mask].flatten(),    weight=w_leading_photon)
                ezfill('photon_eta0',             eta=photons[leadphoton_index].eta[mask].flatten(),  weight=w_leading_photon)
                ezfill('photon_phi0',             phi=photons[leadphoton_index].phi[mask].flatten(),  weight=w_leading_photon)
                ezfill('photon_eta_phi', phi=photons[leadphoton_index].phi[mask].flatten(),eta=photons[leadphoton_index].eta[mask].flatten(), weight=w_leading_photon)

                # w_drphoton_jet = weight_shape(df['dRPhotonJet'][mask], region_weights.partial_weight(exclude=exclude)[mask])

            # PV
            ezfill('npv', nvtx=df['PV_npvs'][mask], weight=region_weights.partial_weight(exclude=exclude)[mask])
            ezfill('npvgood', nvtx=df['PV_npvsGood'][mask], weight=region_weights.partial_weight(exclude=exclude)[mask])

            ezfill('npv_nopu', nvtx=df['PV_npvs'][mask], weight=region_weights.partial_weight(exclude=['pileup']+exclude)[mask])
            ezfill('npvgood_nopu', nvtx=df['PV_npvsGood'][mask], weight=region_weights.partial_weight(exclude=['pileup']+exclude)[mask])

            ezfill('rho_all', rho=df['fixedGridRhoFastjetAll'][mask], weight=region_weights.partial_weight(exclude=exclude)[mask])
            ezfill('rho_central', rho=df['fixedGridRhoFastjetCentral'][mask], weight=region_weights.partial_weight(exclude=exclude)[mask])
            ezfill('rho_all_nopu', rho=df['fixedGridRhoFastjetAll'][mask], weight=region_weights.partial_weight(exclude=['pileup']+exclude)[mask])
            ezfill('rho_central_nopu', rho=df['fixedGridRhoFastjetCentral'][mask], weight=region_weights.partial_weight(exclude=['pileup']+exclude)[mask])
        return output

    def postprocess(self, accumulator):
        return accumulator
