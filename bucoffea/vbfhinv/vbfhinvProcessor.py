import coffea.processor as processor
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
                              weight_shape
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
                                      is_nlo_z
                                      )
from bucoffea.helpers.gen import (
                                  setup_gen_candidates,
                                  setup_dressed_gen_candidates,
                                  setup_lhe_cleaned_genjets,
                                  fill_gen_v_info
                                 )
from bucoffea.monojet.definitions import (
                                          candidate_weights, 
                                          pileup_weights,
                                          setup_candidates, 
                                          theory_weights_vbf
                                          )
from bucoffea.vbfhinv.definitions import (
                                           vbfhinv_accumulator, 
                                           vbfhinv_regions
                                         )

def trigger_selection(selection, df, cfg):
    pass_all = np.zeros(df.size) == 0
    pass_none = ~pass_all
    dataset = df['dataset']
    if cfg.RUN.SYNC: # Synchronization mode
        selection.add('filt_met', pass_all)
        selection.add('trig_met', df['HLT_PFMETNoMu120_PFMHTNoMu120_IDTight'])
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


class vbfhinvProcessor(processor.ProcessorABC):
    def __init__(self, blind=True):
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
            gen_v_pt = df['gen_v_pt_dress']
        elif df['is_lo_g']:
            gen = setup_gen_candidates(df)
            gen_v_pt = gen[(gen.pdg==22) & (gen.status==1)].pt.max()

        # Generator-level leading dijet mass
        if df['has_lhe_v_pt']:
            genjets = setup_lhe_cleaned_genjets(df)
            digenjet = genjets[:,:2].distincts()
            df['mjj_gen'] = digenjet.mass.max()

        # Candidates
        # Already pre-filtered!
        # All leptons are at least loose
        # Check out setup_candidates for filtering details
        met_pt, met_phi, ak4, _, muons, electrons, taus, photons = setup_candidates(df, cfg)

        # Filtering ak4 jets according to pileup ID
        ak4 = ak4[ak4.puid]

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
                            & (electrons.abseta < cfg.ELECTRON.CUTS.TIGHT.ETA)

        dielectrons = electrons.distincts()
        dielectron_charge = dielectrons.i0['charge'] + dielectrons.i1['charge']

        df['MT_el'] = ((electrons.counts==1) * mt(electrons.pt, electrons.phi, met_pt, met_phi)).max()

        # ak4
        leadak4_index=ak4.pt.argmax()

        elejet_pairs = ak4[:,:1].cross(electrons)
        df['dREleJet'] = np.hypot(elejet_pairs.i0.eta-elejet_pairs.i1.eta , dphi(elejet_pairs.i0.phi,elejet_pairs.i1.phi)).min()
        muonjet_pairs = ak4[:,:1].cross(muons)
        df['dRMuonJet'] = np.hypot(muonjet_pairs.i0.eta-muonjet_pairs.i1.eta , dphi(muonjet_pairs.i0.phi,muonjet_pairs.i1.phi)).min()

        # B tagged ak4
        btag_cut = cfg.BTAG.CUTS[cfg.BTAG.algo][cfg.BTAG.wp]
        jet_btag_val = getattr(ak4, cfg.BTAG.algo)
        jet_btagged = jet_btag_val > btag_cut
        bjets = ak4[ (ak4.abseta<2.4) \
                     & jet_btagged \
                     & (ak4.pt>20) ]

        # Recoil
        df['recoil_pt'], df['recoil_phi'] = recoil(met_pt,met_phi, electrons, muons, photons)
        df["dPFCalo"] = (met_pt - df["CaloMET_pt"]) / df["recoil_pt"]
        df["minDPhiJetRecoil"] = min_dphi_jet_met(ak4, df['recoil_phi'], njet=4, ptmin=30, etamax=4.7)
        df["minDPhiJetMet"] = min_dphi_jet_met(ak4, met_phi, njet=4, ptmin=30, etamax=4.7)
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
        selection.add('veto_b', bjets.counts==0)
        selection.add('mindphijr',df['minDPhiJetRecoil'] > cfg.SELECTION.SIGNAL.MINDPHIJR)
        selection.add('dpfcalo',np.abs(df['dPFCalo']) < cfg.SELECTION.SIGNAL.DPFCALO)
        selection.add('recoil', df['recoil_pt']>cfg.SELECTION.SIGNAL.RECOIL)

        if(cfg.MITIGATION.HEM and extract_year(df['dataset']) == 2018 and not cfg.RUN.SYNC):
            selection.add('hemveto', df['hemveto'])
        else:
            selection.add('hemveto', np.ones(df.size)==1)

        # AK4 dijet
        diak4 = ak4[:,:2].distincts()
        leadak4_pt_eta = (diak4.i0.pt > cfg.SELECTION.SIGNAL.LEADAK4.PT) & (np.abs(diak4.i0.eta) < cfg.SELECTION.SIGNAL.LEADAK4.ETA)
        trailak4_pt_eta = (diak4.i1.pt > cfg.SELECTION.SIGNAL.TRAILAK4.PT) & (np.abs(diak4.i1.eta) < cfg.SELECTION.SIGNAL.TRAILAK4.ETA)
        hemisphere = (diak4.i0.eta * diak4.i1.eta < 0).any()
        has_track0 = np.abs(diak4.i0.eta) <= 2.5
        has_track1 = np.abs(diak4.i1.eta) <= 2.5

        leadak4_id = diak4.i0.tightId & (has_track0*((diak4.i0.chf > cfg.SELECTION.SIGNAL.LEADAK4.CHF) & (diak4.i0.nhf < cfg.SELECTION.SIGNAL.LEADAK4.NHF)) + ~has_track0)
        trailak4_id = has_track1*((diak4.i1.chf > cfg.SELECTION.SIGNAL.TRAILAK4.CHF) & (diak4.i1.nhf < cfg.SELECTION.SIGNAL.TRAILAK4.NHF)) + ~has_track1

        df['mjj'] = diak4.mass.max()
        df['dphijj'] = dphi(diak4.i0.phi.min(), diak4.i1.phi.max())
        df['detajj'] = np.abs(diak4.i0.eta - diak4.i1.eta).max()

        selection.add('two_jets', diak4.counts>0)
        selection.add('leadak4_pt_eta', leadak4_pt_eta.any())
        selection.add('trailak4_pt_eta', trailak4_pt_eta.any())
        selection.add('hemisphere', hemisphere)
        selection.add('leadak4_id',leadak4_id.any())
        selection.add('trailak4_id',trailak4_id.any())
        selection.add('mjj', df['mjj'] > cfg.SELECTION.SIGNAL.DIJET.SHAPE_BASED.MASS)
        selection.add('dphijj', df['dphijj'] < cfg.SELECTION.SIGNAL.DIJET.SHAPE_BASED.DPHI)
        selection.add('detajj', df['detajj'] > cfg.SELECTION.SIGNAL.DIJET.SHAPE_BASED.DETA)

        # Divide into three categories for trigger study
        if cfg.RUN.TRIGGER_STUDY:
            two_central_jets = (np.abs(diak4.i0.eta) <= 2.4) & (np.abs(diak4.i1.eta) <= 2.4)
            two_forward_jets = (np.abs(diak4.i0.eta) > 2.4) & (np.abs(diak4.i1.eta) > 2.4)
            one_jet_forward_one_jet_central = (~two_central_jets) & (~two_forward_jets)
            selection.add('two_central_jets', two_central_jets.any())
            selection.add('two_forward_jets', two_forward_jets.any())
            selection.add('one_jet_forward_one_jet_central', one_jet_forward_one_jet_central.any())

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
        selection.add('two_electrons', electrons.counts==2)

        # Single Ele CR
        selection.add('met_el', met_pt > cfg.SELECTION.CONTROL.SINGLEEL.MET)
        selection.add('mt_el', df['MT_el'] < cfg.SELECTION.CONTROL.SINGLEEL.MT)

        # Photon CR
        leadphoton_index=photons.pt.argmax()

        df['is_tight_photon'] = photons.mediumId \
                         & (photons.abseta < cfg.PHOTON.CUTS.TIGHT.ETA)

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

            weights = candidate_weights(weights, df, evaluator, muons, electrons, photons)
            weights = pileup_weights(weights, df, evaluator, cfg)
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
                output['kinematics']['leadbtag'] += [jet_btag_val[(ak4.abseta<2.4) & (ak4.pt>20)][mask].max()]

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
                                  weight=weights.weight()[mask]
                                  )

            fill_mult('ak4_mult', ak4)
            fill_mult('bjet_mult',bjets)
            fill_mult('loose_ele_mult',electrons)
            fill_mult('tight_ele_mult',electrons[df['is_tight_electron']])
            fill_mult('loose_muo_mult',muons)
            fill_mult('tight_muo_mult',muons[df['is_tight_muon']])
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
            for wname, wvalue in weights._weights.items():
                ezfill("weights", weight_type=wname, weight_value=wvalue[mask])
                ezfill("weights_wide", weight_type=wname, weight_value=wvalue[mask])

            # All ak4
            # This is a workaround to create a weight array of the right dimension
            w_alljets = weight_shape(ak4[mask].eta, weights.weight()[mask])
            w_alljets_nopref = weight_shape(ak4[mask].eta, weights.partial_weight(exclude=['prefire'])[mask])

            ezfill('ak4_eta',    jeteta=ak4[mask].eta.flatten(), weight=w_alljets)
            ezfill('ak4_phi',    jetphi=ak4[mask].phi.flatten(), weight=w_alljets)
            ezfill('ak4_pt',     jetpt=ak4[mask].pt.flatten(),   weight=w_alljets)

            ezfill('ak4_eta_nopref',    jeteta=ak4[mask].eta.flatten(), weight=w_alljets_nopref)
            ezfill('ak4_phi_nopref',    jetphi=ak4[mask].phi.flatten(), weight=w_alljets_nopref)
            ezfill('ak4_pt_nopref',     jetpt=ak4[mask].pt.flatten(),   weight=w_alljets_nopref)

            # Leading ak4
            w_diak4 = weight_shape(diak4.pt[mask], weights.weight()[mask])
            ezfill('ak4_eta0',      jeteta=diak4.i0.eta[mask].flatten(),    weight=w_diak4)
            ezfill('ak4_phi0',      jetphi=diak4.i0.phi[mask].flatten(),    weight=w_diak4)
            ezfill('ak4_pt0',       jetpt=diak4.i0.pt[mask].flatten(),      weight=w_diak4)
            ezfill('ak4_ptraw0',    jetpt=diak4.i0.ptraw[mask].flatten(),   weight=w_diak4)
            ezfill('ak4_chf0',      frac=diak4.i0.chf[mask].flatten(),      weight=w_diak4)
            ezfill('ak4_nhf0',      frac=diak4.i0.nhf[mask].flatten(),      weight=w_diak4)
            ezfill('ak4_nconst0',   nconst=diak4.i0.nconst[mask].flatten(), weight=w_diak4)

            # Trailing ak4
            ezfill('ak4_eta1',      jeteta=diak4.i1.eta[mask].flatten(),    weight=w_diak4)
            ezfill('ak4_phi1',      jetphi=diak4.i1.phi[mask].flatten(),    weight=w_diak4)
            ezfill('ak4_pt1',       jetpt=diak4.i1.pt[mask].flatten(),      weight=w_diak4)
            ezfill('ak4_ptraw1',    jetpt=diak4.i1.ptraw[mask].flatten(),   weight=w_diak4)
            ezfill('ak4_chf1',      frac=diak4.i1.chf[mask].flatten(),      weight=w_diak4)
            ezfill('ak4_nhf1',      frac=diak4.i1.nhf[mask].flatten(),      weight=w_diak4)
            ezfill('ak4_nconst1',   nconst=diak4.i1.nconst[mask].flatten(), weight=w_diak4)

            # B tag discriminator
            btag = getattr(ak4, cfg.BTAG.ALGO)
            w_btag = weight_shape(btag[mask], weights.weight()[mask])
            ezfill('ak4_btag', btag=btag[mask].flatten(), weight=w_btag )

            # MET
            ezfill('dpfcalo',            dpfcalo=df["dPFCalo"][mask],       weight=weights.weight()[mask] )
            ezfill('met',                met=met_pt[mask],            weight=weights.weight()[mask] )
            ezfill('met_phi',            phi=met_phi[mask],           weight=weights.weight()[mask] )
            ezfill('recoil',             recoil=df["recoil_pt"][mask],      weight=weights.weight()[mask] )
            ezfill('recoil_phi',         phi=df["recoil_phi"][mask],        weight=weights.weight()[mask] )
            ezfill('dphijm',             dphi=df["minDPhiJetMet"][mask],    weight=weights.weight()[mask] )
            ezfill('dphijr',             dphi=df["minDPhiJetRecoil"][mask], weight=weights.weight()[mask] )

            ezfill('dphijj',             dphi=df["dphijj"][mask],   weight=weights.weight()[mask] )
            ezfill('detajj',             deta=df["detajj"][mask],   weight=weights.weight()[mask] )
            ezfill('mjj',                mjj=df["mjj"][mask],      weight=weights.weight()[mask] )

            # Two dimensional
            ezfill('recoil_mjj',         recoil=df["recoil_pt"][mask], mjj=df["mjj"][mask], weight=weights.weight()[mask] )

            # Muons
            if '_1m_' in region or '_2m_' in region:
                w_allmu = weight_shape(muons.pt[mask], weights.weight()[mask])
                ezfill('muon_pt',   pt=muons.pt[mask].flatten(),    weight=w_allmu )
                ezfill('muon_mt',   mt=df['MT_mu'][mask],           weight=weights.weight()[mask])
                ezfill('muon_eta',  eta=muons.eta[mask].flatten(),  weight=w_allmu)
                ezfill('muon_phi',  phi=muons.phi[mask].flatten(),  weight=w_allmu)

            # Dimuon
            if '_2m_' in region:
                w_dimu = weight_shape(dimuons.pt[mask], weights.weight()[mask])
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
            if '_1e_' in region or '_2e_' in region:
                w_allel = weight_shape(electrons.pt[mask], weights.weight()[mask])
                ezfill('electron_pt',   pt=electrons.pt[mask].flatten(),    weight=w_allel)
                ezfill('electron_mt',   mt=df['MT_el'][mask],               weight=weights.weight()[mask])
                ezfill('electron_eta',  eta=electrons.eta[mask].flatten(),  weight=w_allel)
                ezfill('electron_phi',  phi=electrons.phi[mask].flatten(),  weight=w_allel)

            # Dielectron
            if '_2e_' in region:
                w_diel = weight_shape(dielectrons.pt[mask], weights.weight()[mask])
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
                w_leading_photon = weight_shape(photons[leadphoton_index].pt[mask],weights.weight()[mask]);
                ezfill('photon_pt0',              pt=photons[leadphoton_index].pt[mask].flatten(),    weight=w_leading_photon)
                ezfill('photon_eta0',             eta=photons[leadphoton_index].eta[mask].flatten(),  weight=w_leading_photon)
                ezfill('photon_phi0',             phi=photons[leadphoton_index].phi[mask].flatten(),  weight=w_leading_photon)
                ezfill('photon_pt0_recoil',       pt=photons[leadphoton_index].pt[mask].flatten(), recoil=df['recoil_pt'][mask&(leadphoton_index.counts>0)],  weight=w_leading_photon)
                ezfill('photon_eta_phi',          eta=photons[leadphoton_index].eta[mask].flatten(), phi=photons[leadphoton_index].phi[mask].flatten(),  weight=w_leading_photon)

                # w_drphoton_jet = weight_shape(df['dRPhotonJet'][mask], weights.weight()[mask])


            # PV
            ezfill('npv', nvtx=df['PV_npvs'][mask], weight=weights.weight()[mask])
            ezfill('npvgood', nvtx=df['PV_npvsGood'][mask], weight=weights.weight()[mask])

            ezfill('npv_nopu', nvtx=df['PV_npvs'][mask], weight=weights.partial_weight(exclude=['pileup'])[mask])
            ezfill('npvgood_nopu', nvtx=df['PV_npvsGood'][mask], weight=weights.partial_weight(exclude=['pileup'])[mask])

            ezfill('rho_all', rho=df['fixedGridRhoFastjetAll'][mask], weight=weights.weight()[mask])
            ezfill('rho_central', rho=df['fixedGridRhoFastjetCentral'][mask], weight=weights.weight()[mask])
            ezfill('rho_all_nopu', rho=df['fixedGridRhoFastjetAll'][mask], weight=weights.partial_weight(exclude=['pileup'])[mask])
            ezfill('rho_central_nopu', rho=df['fixedGridRhoFastjetCentral'][mask], weight=weights.partial_weight(exclude=['pileup'])[mask])
        return output

    def postprocess(self, accumulator):
        return accumulator
