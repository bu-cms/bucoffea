import copy
import coffea.processor as processor
import re
import numpy as np
import copy
from dynaconf import settings as cfg

from bucoffea.helpers import (
                              bucoffea_path,
                              dphi,
                              evaluator_from_config,
                              mask_and, 
                              mask_or, 
                              min_dphi_jet_met,
                              mjj,
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
                                          theory_weights_vbf,
                                          photon_trigger_sf,
                                          photon_impurity_weights,
                                          data_driven_qcd_dataset
                                          )
from bucoffea.vbfhinv.definitions import (
                                           vbfhinv_accumulator,
                                           vbfhinv_regions
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




def get_veto_weights(df, evaluator, electrons, muons, taus):
    """
    Calculate veto weights for SR W

    The weights are effectively:

        w = product(1-SF)

    where the product runs overveto-able e, mu, tau.
    """
    veto_weights = processor.Weights(size=df.size, storeIndividual=True)

    for variation in [
                      'nominal',
                    #   'ele_reco_up','ele_reco_dn',
                    #   'ele_id_up','ele_id_dn',
                    #   'muon_id_up','muon_id_dn',
                    #   'muon_iso_up','muon_iso_dn',
                    #   'tau_id_up','tau_id_dn'
                      ]:

        def varied_weight(sfname, *args):
            '''Helper function to easily get the correct weights for a given variation'''

            # For the nominal variation, just pass through
            if 'nominal' in variation:
                return evaluator[sfname](*args)

            # If this variation is unrelated to the SF at hand,
            # pass through as well
            if not (re.sub('_(up|dn)', '', variation) in sfname):
                return evaluator[sfname](*args)

            # Direction of variation
            sgn = 1 if variation.endswith("up") else -1
            return evaluator[sfname](*args) + sgn * evaluator[f"{sfname}_error"](*args)


        ### Electrons
        if extract_year(df['dataset']) == 2017:
            high_et = electrons.pt>20

            # Low pt SFs
            low_pt_args = (electrons.etasc[~high_et], electrons.pt[~high_et])
            ele_reco_sf_low = varied_weight('ele_reco_pt_lt_20', *low_pt_args)
            ele_id_sf_low = varied_weight("ele_id_loose", *low_pt_args)

            # High pt SFs
            high_pt_args = (electrons.etasc[high_et], electrons.pt[high_et])

            ele_reco_sf_high = varied_weight("ele_reco", *high_pt_args)
            ele_id_sf_high = varied_weight("ele_id_loose", *high_pt_args)

            # Combine
            veto_weight_ele = (1 - ele_reco_sf_low*ele_id_sf_low).prod() * (1-ele_reco_sf_high*ele_id_sf_high).prod()
        else:
            # No split for 2018
            args = (electrons.etasc, electrons.pt)
            ele_reco_sf = varied_weight("ele_reco", *args)
            ele_id_sf = varied_weight("ele_id_loose", *args)

            # Combine
            veto_weight_ele = (1 - ele_id_sf*ele_reco_sf).prod()

        ### Muons
        args = (muons.pt, muons.abseta)
        veto_weight_muo = (1 - varied_weight("muon_id_loose", *args)*varied_weight("muon_iso_loose", *args)).prod()

        ### Taus
        # Taus have their variations saves as separate histograms,
        # so our cool trick from above is replaced by the pedestrian way
        if "tau_id" in variation:
            direction = variation.split("_")[-1]
            tau_sf_name = f"tau_id_{direction}"
        else:
            tau_sf_name = "tau_id"
        veto_weight_tau = (1 - evaluator[tau_sf_name](taus.pt)).prod()

        ### Combine
        veto_weights.add(variation, veto_weight_ele * veto_weight_muo * veto_weight_tau)

    return veto_weights



class vbfhinvProcessor(processor.ProcessorABC):
    def __init__(self, blind=True):
        self._year=None
        self._blind=blind
        self._configure()
        self._variations = ['', '_jerup', '_jerdown', '_jesup', '_jesdown']
        self._accumulator = vbfhinv_accumulator(cfg, variations=self._variations)

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
        df['is_lo_w_ewk'] = is_lo_w_ewk(dataset)
        df['is_lo_z_ewk'] = is_lo_z_ewk(dataset)
        df['is_lo_g'] = is_lo_g(dataset)
        df['is_nlo_z'] = is_nlo_z(dataset)
        df['is_nlo_w'] = is_nlo_w(dataset)
        df['has_lhe_v_pt'] = df['is_lo_w'] | df['is_lo_z'] | df['is_nlo_z'] | df['is_nlo_w'] | df['is_lo_g'] | df['is_lo_w_ewk'] | df['is_lo_z_ewk']
        df['is_data'] = is_data(dataset)

        if df['is_data']:
            return self.accumulator.identity()

        gen_v_pt = None
        if df['is_lo_w'] or df['is_lo_z'] or df['is_nlo_z'] or df['is_nlo_w'] or df['is_lo_z_ewk'] or df['is_lo_w_ewk']:
            gen = setup_gen_candidates(df)
            dressed = setup_dressed_gen_candidates(df)
            fill_gen_v_info(df, gen, dressed)
            gen_v_pt = df['gen_v_pt_combined']
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
        vmap, _, muons, electrons, taus, photons = setup_candidates(df, cfg, variations=self._variations)

        # vmap holds information about ak4, met and selection
        # packers for each JES/JER variation 
        # Check out monojet/definitions.py for the object definition

        #################################
        # First process the part which is 
        # unrelated to JES/JER variations
        #################################
        
        # Muons
        df['is_tight_muon'] = muons.tightId \
                      & (muons.iso < cfg.MUON.CUTS.TIGHT.ISO) \
                      & (muons.pt>cfg.MUON.CUTS.TIGHT.PT) \
                      & (muons.abseta<cfg.MUON.CUTS.TIGHT.ETA)

        dimuons = muons.distincts()
        dimuon_charge = dimuons.i0['charge'] + dimuons.i1['charge']
        
        # Electrons
        df['is_tight_electron'] = electrons.tightId \
                            & (electrons.pt > cfg.ELECTRON.CUTS.TIGHT.PT) \
                            & (electrons.absetasc < cfg.ELECTRON.CUTS.TIGHT.ETA)

        dielectrons = electrons.distincts()
        dielectron_charge = dielectrons.i0['charge'] + dielectrons.i1['charge']

        # Selection packer for the nominal (no varition) case
        selection_nom = processor.PackedSelection() 

        # Triggers
        pass_all = np.ones(df.size)==1
        selection_nom.add('inclusive', pass_all)
        selection_nom = trigger_selection(selection_nom, df, cfg)
        
        selection_nom.add('mu_pt_trig_safe', muons.pt.max() > 30)

        # Common selection
        selection_nom.add('veto_ele', electrons.counts==0)
        selection_nom.add('veto_muo', muons.counts==0)
        selection_nom.add('veto_photon', photons.counts==0)
        selection_nom.add('veto_tau', taus.counts==0)

        if(cfg.MITIGATION.HEM and extract_year(df['dataset']) == 2018 and not cfg.RUN.SYNC):
            selection_nom.add('hemveto', df['hemveto'])
        else:
            selection_nom.add('hemveto', np.ones(df.size)==1)
            
        selection_nom.add('one_muon', muons.counts==1)

        # Dimuon CR
        leadmuon_index=muons.pt.argmax()
        selection_nom.add('at_least_one_tight_mu', df['is_tight_muon'].any())
        selection_nom.add('dimuon_mass', ((dimuons.mass > cfg.SELECTION.CONTROL.DOUBLEMU.MASS.MIN) \
                                    & (dimuons.mass < cfg.SELECTION.CONTROL.DOUBLEMU.MASS.MAX)).any())
        selection_nom.add('dimuon_charge', (dimuon_charge==0).any())
        selection_nom.add('two_muons', muons.counts==2)

        # Diele CR
        leadelectron_index=electrons.pt.argmax()
        selection_nom.add('one_electron', electrons.counts==1)
        selection_nom.add('two_electrons', electrons.counts==2)
        selection_nom.add('at_least_one_tight_el', df['is_tight_electron'].any())

        selection_nom.add('dielectron_mass', ((dielectrons.mass > cfg.SELECTION.CONTROL.DOUBLEEL.MASS.MIN)  \
                                        & (dielectrons.mass < cfg.SELECTION.CONTROL.DOUBLEEL.MASS.MAX)).any())
        selection_nom.add('dielectron_charge', (dielectron_charge==0).any())
        selection_nom.add('two_electrons', electrons.counts==2)

        df['is_tight_photon'] = photons.mediumId & photons.barrel

        # Photon CR
        leadphoton_index=photons.pt.argmax()
        selection_nom.add('one_photon', photons.counts==1)
        selection_nom.add('at_least_one_tight_photon', df['is_tight_photon'].any())
        selection_nom.add('photon_pt', photons.pt.max() > cfg.PHOTON.CUTS.TIGHT.PT)
        selection_nom.add('photon_pt_trig', photons.pt.max() > cfg.PHOTON.CUTS.TIGHT.PTTRIG)
     
        vmap.set_selection_packer(var='', sel=selection_nom)

        # Process for each JES/JER variation
        for var in self._variations:
            # Get the correct objects/quantities for each variation
            # For other variations, copy the common selections and
            # add on top of those.
            if var == '':
                selection = vmap.get_selection_packer(var='')
            else:
                selection = copy.deepcopy(selection_nom) 
                vmap.set_selection_packer(var=var, sel=selection)    

            bjets = vmap.get_bjets(var)
            ak4 = vmap.get_ak4(var) 
            diak4 = vmap.get_diak4(var) 
            met = vmap.get_met(var) 

            met_pt = getattr(met, f'pt{var}').flatten()
            met_phi = getattr(met, f'phi{var}').flatten()

            selection.add(f'veto_b{var}', bjets.counts==0)

            # Filtering ak4 jets according to pileup ID
            ak4_puid = getattr(ak4, f'puid{var}')
            bjets_puid = getattr(bjets, f'puid{var}')
            
            ak4 = ak4[ak4_puid]
            bjets = bjets[bjets_puid] 

            ak4_pt = getattr(ak4, f'pt{var}')

            df[f'MT_mu{var}'] = ((muons.counts==1) * mt(muons.pt, muons.phi, met_pt, met_phi)).max()
            selection.add(f'mt_mu{var}', df[f'MT_mu{var}'] < cfg.SELECTION.CONTROL.SINGLEMU.MT)
            
            df[f'MT_el{var}'] = ((electrons.counts==1) * mt(electrons.pt, electrons.phi, met_pt, met_phi)).max()
            selection.add(f'mt_el{var}', df[f'MT_el{var}'] < cfg.SELECTION.CONTROL.SINGLEEL.MT)

            selection.add(f'met_el{var}', met_pt > cfg.SELECTION.CONTROL.SINGLEEL.MET)
           
            # Leading jet
            leadak4_index=ak4_pt.argmax()

            elejet_pairs = ak4[:,:1].cross(electrons)
            df[f'dREleJet{var}'] = np.hypot(elejet_pairs.i0.eta-elejet_pairs.i1.eta , dphi(elejet_pairs.i0.phi,elejet_pairs.i1.phi)).min()
            muonjet_pairs = ak4[:,:1].cross(muons)
            df['dRMuonJet{var}'] = np.hypot(muonjet_pairs.i0.eta-muonjet_pairs.i1.eta , dphi(muonjet_pairs.i0.phi,muonjet_pairs.i1.phi)).min()

            # Recoil
            df[f'recoil_pt{var}'], df[f'recoil_phi{var}'] = recoil(met_pt, met_phi, electrons, muons, photons)
            df[f"dPFCalo{var}"] = (met_pt - df["CaloMET_pt"]) / df[f"recoil_pt{var}"]
            df[f"minDPhiJetRecoil{var}"] = min_dphi_jet_met(ak4, df[f'recoil_phi{var}'], njet=4, ptmin=30, etamax=5.0, var=var)
            df[f"minDPhiJetMet{var}"] = min_dphi_jet_met(ak4, met_phi, njet=4, ptmin=30, etamax=5.0, var=var)
        
            selection.add(f'mindphijr{var}',df[f'minDPhiJetRecoil{var}'] > cfg.SELECTION.SIGNAL.MINDPHIJR)
            selection.add(f'dpfcalo{var}',np.abs(df[f'dPFCalo{var}']) < cfg.SELECTION.SIGNAL.DPFCALO)
            selection.add(f'recoil{var}', df[f'recoil_pt{var}']>cfg.SELECTION.SIGNAL.RECOIL)

            # AK4 dijet
            lead_jet_pt = getattr(diak4.i0, f'pt{var}')
            trail_jet_pt = getattr(diak4.i1, f'pt{var}')

            leadak4_pt_eta = (lead_jet_pt > cfg.SELECTION.SIGNAL.LEADAK4.PT) & (np.abs(diak4.i0.eta) < cfg.SELECTION.SIGNAL.LEADAK4.ETA)
            trailak4_pt_eta = (trail_jet_pt > cfg.SELECTION.SIGNAL.TRAILAK4.PT) & (np.abs(diak4.i1.eta) < cfg.SELECTION.SIGNAL.TRAILAK4.ETA)
            hemisphere = (diak4.i0.eta * diak4.i1.eta < 0).any()
            has_track0 = np.abs(diak4.i0.eta) <= 2.5
            has_track1 = np.abs(diak4.i1.eta) <= 2.5

            leadak4_id = diak4.i0.tightId & (has_track0*((diak4.i0.chf > cfg.SELECTION.SIGNAL.LEADAK4.CHF) & (diak4.i0.nhf < cfg.SELECTION.SIGNAL.LEADAK4.NHF)) + ~has_track0)
            trailak4_id = has_track1*((diak4.i1.chf > cfg.SELECTION.SIGNAL.TRAILAK4.CHF) & (diak4.i1.nhf < cfg.SELECTION.SIGNAL.TRAILAK4.NHF)) + ~has_track1

            df[f'mjj{var}'] = mjj(diak4, var=var)
            df[f'dphijj{var}'] = dphi(diak4.i0.phi.min(), diak4.i1.phi.max())
            df[f'detajj{var}'] = np.abs(diak4.i0.eta - diak4.i1.eta).max()

            selection.add(f'two_jets{var}', diak4.counts>0)
            selection.add(f'leadak4_pt_eta{var}', leadak4_pt_eta.any())
            selection.add(f'trailak4_pt_eta{var}', trailak4_pt_eta.any())
            selection.add(f'hemisphere{var}', hemisphere)
            selection.add(f'leadak4_id{var}',leadak4_id.any())
            selection.add(f'trailak4_id{var}',trailak4_id.any())
            selection.add(f'mjj{var}', df[f'mjj{var}'] > cfg.SELECTION.SIGNAL.DIJET.SHAPE_BASED.MASS)
            selection.add(f'dphijj{var}', df[f'dphijj{var}'] < cfg.SELECTION.SIGNAL.DIJET.SHAPE_BASED.DPHI)
            selection.add(f'detajj{var}', df[f'detajj{var}'] > cfg.SELECTION.SIGNAL.DIJET.SHAPE_BASED.DETA)

            # Calculate ratios: Varied / Nominal
            if var != '':
                df[f'mjj{var}_over_nom'] = df[f'mjj{var}']/df['mjj'] - 1 
                df[f'detajj{var}_over_nom'] = df[f'detajj{var}']/df['detajj'] - 1 
                df[f'dphijj{var}_over_nom'] = df[f'dphijj{var}']/df['dphijj'] - 1 
                df[f'recoil_pt{var}_over_nom'] = df[f'recoil_pt{var}']/df['recoil_pt'] - 1 

                # Get nominal leading and trailing jet pt
                diak4_nom = vmap.get_diak4(var='') 
                lead_jetpt_nom = diak4_nom.i0.pt 
                trail_jetpt_nom = diak4_nom.i1.pt 

                df[f'ak4_pt0{var}_over_nom'] = (lead_jet_pt / lead_jetpt_nom - 1).flatten()
                df[f'ak4_pt1{var}_over_nom'] = (trail_jet_pt / trail_jetpt_nom - 1).flatten()

        # Divide into three categories for trigger study
        if cfg.RUN.TRIGGER_STUDY:
            two_central_jets = (np.abs(diak4.i0.eta) <= 2.4) & (np.abs(diak4.i1.eta) <= 2.4)
            two_forward_jets = (np.abs(diak4.i0.eta) > 2.4) & (np.abs(diak4.i1.eta) > 2.4)
            one_jet_forward_one_jet_central = (~two_central_jets) & (~two_forward_jets)
            selection.add('two_central_jets', two_central_jets.any())
            selection.add('two_forward_jets', two_forward_jets.any())
            selection.add('one_jet_forward_one_jet_central', one_jet_forward_one_jet_central.any())

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

        veto_weights = get_veto_weights(df, evaluator, electrons, muons, taus)
        regions = vbfhinv_regions(cfg, variations=self._variations)

        for region, cuts in regions.items():
            exclude = [None]
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
                    region_weights.add('trigger_met', evaluator["trigger_met"](df['recoil_pt']))
                elif re.match(r'cr_g.*', region):
                    photon_trigger_sf(region_weights, photons, df)

                # Veto weights
                if re.match('.*no_veto.*', region):
                    exclude = [
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


            # This is the default weight for this region
            rweight = region_weights.partial_weight(exclude=exclude)

            # Blinding
            if(self._blind and df['is_data'] and region.startswith('sr')):
                continue

            # Get relevant variation for each region
            if ('up' in region) or ('down' in region):
                var = '_' + region.split('_')[-1]
            else:
                var = ''

            # Get the correct objects/quantities for each variation
            selection = vmap.get_selection_packer(var)
            bjets = vmap.get_bjets(var)
            ak4 = vmap.get_ak4(var) 
            diak4 = vmap.get_diak4(var) 
            met = vmap.get_met(var)

            met_pt_nom = met.pt_nom.flatten() # MET_pt_nom in nanoAOD
            met_pt_jer = met.pt.flatten()     # MET_pt_jer in nanoAOD
            met_pt = getattr(met, f'pt{var}').flatten() # Varied MET pt 

            # Cutflow plot for signal and control regions
            if any(x in region for x in ["sr", "cr", "tr"]):
                output['cutflow_' + region][dataset]['all']+=df.size
                # Get weighted cutflow
                for icut, cutname in enumerate(cuts):
                    output['cutflow_' + region][dataset][cutname] += np.nansum(region_weights.weight()[selection.all(*cuts[:icut+1])] )

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
                                  weight=rweight[mask]
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
            for wname, wvalue in region_weights._weights.items():
                ezfill("weights", weight_type=wname, weight_value=wvalue[mask])
                ezfill("weights_wide", weight_type=wname, weight_value=wvalue[mask])

            # All ak4
            # This is a workaround to create a weight array of the right dimension
            w_alljets = weight_shape(ak4[mask].eta, rweight[mask])
            
            ezfill('ak4_pt',     jetpt=getattr(ak4, f'pt{var}')[mask].flatten(),   weight=w_alljets)

            # Leading ak4
            w_diak4 = weight_shape(diak4.pt[mask], rweight[mask])
            ezfill('ak4_pt0',       jetpt=getattr(diak4.i0, f'pt{var}')[mask].flatten(),      weight=w_diak4)

            # Trailing ak4
            ezfill('ak4_pt1',       jetpt=getattr(diak4.i1, f'pt{var}')[mask].flatten(),      weight=w_diak4)

            # B tag discriminator
            btag = getattr(ak4, cfg.BTAG.ALGO)
            w_btag = weight_shape(btag[mask], rweight[mask])

            # Photon CR data-driven QCD estimate
            if df['is_data'] and re.match("cr_g.*", region) and re.match("(SinglePhoton|EGamma).*", dataset):
                w_imp = photon_impurity_weights(photons[leadphoton_index].pt.max()[mask], df["year"])
                output['mjj'].fill(
                                    dataset=data_driven_qcd_dataset(dataset),
                                    region=region,
                                    mjj=df["mjj"][mask],
                                    weight=rweight[mask] * w_imp
                                )

            # MET
            ezfill('met',                met=met_pt[mask],            weight=rweight[mask] )
            ezfill('met_inc',            met=met_pt,                  weight=rweight )
            ezfill('recoil',             recoil=df[f"recoil_pt{var}"][mask],      weight=rweight[mask] )
            ezfill('mjj',                mjj=df[f"mjj{var}"][mask],      weight=rweight[mask] )

            # Inclusive and masked, JER smeared and nominal (not JER smeared) MET
            if region in ['sr_vbf']:
                ezfill('met_jer',       met=met_pt_jer[mask],       weight=rweight[mask])
                ezfill('met_nom',       met=met_pt_nom[mask],       weight=rweight[mask])
                ezfill('met_jer_inc',   met=met_pt_jer,       weight=rweight)
                ezfill('met_nom_inc',   met=met_pt_nom,       weight=rweight)

            # Muons
            if region in ['cr_1m_vbf', 'cr_2m_vbf']:
                w_allmu = weight_shape(muons.pt[mask], rweight[mask])
                ezfill('muon_pt',   pt=muons.pt[mask].flatten(),    weight=w_allmu)
                ezfill('muon_eta',  eta=muons.eta[mask].flatten(),  weight=w_allmu)
                ezfill('muon_phi',  phi=muons.phi[mask].flatten(),  weight=w_allmu)

            # Dimuon
            if region in ['cr_2m_vbf']:
                w_dimu = weight_shape(dimuons.pt[mask], rweight[mask])
                ezfill('muon_pt0',      pt=dimuons.i0.pt[mask].flatten(),           weight=w_dimu)
                ezfill('muon_pt1',      pt=dimuons.i1.pt[mask].flatten(),           weight=w_dimu)
                ezfill('muon_eta0',     eta=dimuons.i0.eta[mask].flatten(),         weight=w_dimu)
                ezfill('muon_eta1',     eta=dimuons.i1.eta[mask].flatten(),         weight=w_dimu)
                ezfill('muon_phi0',     phi=dimuons.i0.phi[mask].flatten(),         weight=w_dimu)
                ezfill('muon_phi1',     phi=dimuons.i1.phi[mask].flatten(),         weight=w_dimu)
                ezfill('dimuon_pt',     pt=dimuons.pt[mask].flatten(),              weight=w_dimu)
                ezfill('dimuon_eta',    eta=dimuons.eta[mask].flatten(),            weight=w_dimu)
                ezfill('dimuon_mass',   dilepton_mass=dimuons.mass[mask].flatten(), weight=w_dimu)

            # Electrons
            if region in ['cr_1e_vbf', 'cr_2e_vbf']:
                w_allel = weight_shape(electrons.pt[mask], rweight[mask])
                ezfill('electron_pt',   pt=electrons.pt[mask].flatten(),    weight=w_allel)
                ezfill('electron_eta',  eta=electrons.eta[mask].flatten(),  weight=w_allel)
                ezfill('electron_phi',  phi=electrons.phi[mask].flatten(),  weight=w_allel)

            # Dielectron
            if region in ['cr_2e_vbf']:
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
            if region in ['cr_g_vbf']:
                w_leading_photon = weight_shape(photons[leadphoton_index].pt[mask],rweight[mask])
                ezfill('photon_pt0',              pt=photons[leadphoton_index].pt[mask].flatten(),    weight=w_leading_photon)
                ezfill('photon_eta0',             eta=photons[leadphoton_index].eta[mask].flatten(),  weight=w_leading_photon)
                ezfill('photon_phi0',             phi=photons[leadphoton_index].phi[mask].flatten(),  weight=w_leading_photon)
                #ezfill('photon_pt0_recoil',       pt=photons[leadphoton_index].pt[mask].flatten(), recoil=df['recoil_pt'][mask&(leadphoton_index.counts>0)],  weight=w_leading_photon)
                #ezfill('photon_eta_phi',          eta=photons[leadphoton_index].eta[mask].flatten(), phi=photons[leadphoton_index].phi[mask].flatten(),  weight=w_leading_photon)

                # w_drphoton_jet = weight_shape(df['dRPhotonJet'][mask], rweight[mask])

            # Tau
            if 'no_veto' in region:
                w_all_taus = weight_shape(taus.pt[mask], rweight[mask])
                ezfill("tau_pt", pt=taus.pt[mask].flatten(), weight=w_all_taus)

            # Variation / Nominal ratio plots for signal region
            if region.startswith('sr') and var != '':
                ezfill('recoil_varovernom',       ratio=df[f'recoil_pt{var}_over_nom'], weight=rweight)             
                ezfill('mjj_varovernom',          ratio=df[f'mjj{var}_over_nom'],    weight=rweight)             
                ezfill('detajj_varovernom',       ratio=df[f'detajj{var}_over_nom'], weight=rweight)             
                ezfill('dphijj_varovernom',       ratio=df[f'dphijj{var}_over_nom'], weight=rweight)             

            # PV
            if region in ['sr_vbf', 'cr_1m_vbf', 'cr_2m_vbf', 'cr_1e_vbf', 'cr_2e_vbf']:
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
