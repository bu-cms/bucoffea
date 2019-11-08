import os
import re

import coffea.processor as processor
import numpy as np
from awkward import JaggedArray
from coffea import hist
from coffea.analysis_objects import JaggedCandidateArray
from coffea.lumi_tools import LumiMask

from bucoffea.helpers import bucoffea_path,min_dphi_jet_met, object_overlap, weight_shape, mask_and
from bucoffea.helpers.dataset import (extract_year, is_data, is_lo_w, is_lo_z,
                                      is_nlo_w, is_nlo_z)
from bucoffea.helpers.gen import (fill_gen_v_info, find_gen_dilepton, islep,
                                  isnu, setup_dressed_gen_candidates,
                                  setup_gen_candidates)

from dynaconf import settings as cfg

Hist = hist.Hist
Bin = hist.Bin
Cat = hist.Cat


# VID map definition
# 
# 0  1  MinPtCut
# 2  3  PhoSCEtaMultiRangeCut
# 4  5  PhoSingleTowerHadOverEmCut
# 6  7  PhoFull5x5SigmaIEtaIEtaCut
# 8  9  PhoAnyPFIsoWithEACut
# 10 11 PhoAnyPFIsoWithEAAndQuadScalingCut
# 12 13 PhoAnyPFIsoWithEACut

def medium_id_no_sieie(photons):
    # VID bitmask encodes 7 cuts with each two bits.
    # for each cut, the more significant bit indicates
    # whether the medium WP of the cut is passed.
    #
    # To get medium or tighter, require every second bit:
    #   --> 1X-1X-1X-1X-1X-1X-1X
    #   where X = 0 or 1 (we dont care)
    #
    # To ignore the sieie requirement, also ignore bit 7
    #   --> 1X-1X-1X-XX-1X-1X-1X
    # 
    # We are going to check the logical AND of the mask
    # and the bitmap, so set all X to 0
    mask = int('10101000101010',2)
    return (photons.vid & mask) == mask

def medium_id_no_sieie_inv_iso(photons):
    # Same as medium_id_no_sieie, but in first step,
    # ignore isolation bits
    # --> XX-XX-XX-XX-1X-1X-1X
    mask1 = int('00000000101010',2)
    medium_id_no_iso = (photons.vid & mask1) == mask1

    # In second step, require that at least one
    # of the most significant isolation bits
    # fails, which we achieve by using a mask
    # that would *pass* and then requiring that
    # (vid & mask) != mask
    mask2 = int('10101000000000',2)
    inv_iso = (photons.vid & mask2) != mask2

    return medium_id_no_iso & inv_iso


def setup_photons(df):
    # Setup photons

    if extract_year(df['dataset']) == 2016:
        id_branch = 'Photon_cutBased'
    else:
        id_branch = 'Photon_cutBasedBitmap'


    photons = JaggedCandidateArray.candidatesfromcounts(
        df['nPhoton'],
        pt=df['Photon_pt'],
        eta=df['Photon_eta'],
        abseta=np.abs(df['Photon_eta']),
        phi=df['Photon_phi'],
        mass=0*df['Photon_pt'],
        mediumId=(df[id_branch]>=2) & df['Photon_electronVeto'],
        r9=df['Photon_r9'],
        barrel=np.abs(df['Photon_eta']) < 1.479,
        vid=df['Photon_vidNestedWPBitmap'],
        eleveto= df['Photon_electronVeto'],
        sieie= df['Photon_sieie'],
    )

    photons = photons[ (photons.pt > 200) & photons.barrel & photons.eleveto]
    return photons

def setup_jets(df):
    ak4 = JaggedCandidateArray.candidatesfromcounts(
        df['nJet'],
        pt=df[f'Jet_pt'],
        eta=df['Jet_eta'],
        abseta=np.abs(df['Jet_eta']),
        phi=df['Jet_phi'],
        mass=np.zeros_like(df['Jet_pt']),
        looseId=(df['Jet_jetId']&2) == 2, # bitmask: 1 = loose, 2 = tight, 3 = tight + lep veto
        tightId=(df['Jet_jetId']&2) == 2, # bitmask: 1 = loose, 2 = tight, 3 = tight + lep veto
        nhf=df['Jet_neHEF'],
        chf=df['Jet_chHEF'],
        ptraw=df['Jet_pt']*(1-df['Jet_rawFactor']),
        nconst=df['Jet_nConstituents']
    )
    return ak4

class photonPurityProcessor(processor.ProcessorABC):
    def __init__(self):

        # Histogram setup
        dataset_ax = Cat("dataset", "Primary dataset")

        sieie_ax = Bin("sieie", r"sieie", 100,0,0.02)
        pt_ax = Bin("pt", r"pt", 50, 200, 1200)
        cat_ax = Cat("cat", r"cat")

        items = {}
        items[f"sieie"] = Hist(
                               "Counts",
                               dataset_ax,
                               sieie_ax,
                               pt_ax,
                               cat_ax
                               )

        items['sumw'] = processor.defaultdict_accumulator(float)
        items['sumw2'] = processor.defaultdict_accumulator(float)

        self._accumulator = processor.dict_accumulator(items)
        self._configure()

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

    @property
    def accumulator(self):
        return self._accumulator


    def process(self, df):
        self._configure(df)
        output = self.accumulator.identity()
        dataset = df['dataset']

        # Lumi mask
        if is_data(dataset):
            year = extract_year(dataset)
            if year == 2016:
                json = bucoffea_path('data/json/Cert_271036-284044_13TeV_ReReco_07Aug2017_Collisions16_JSON.txt')
            elif year == 2017:
                json = bucoffea_path('data/json/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1.txt')
            elif year == 2018:
                json = bucoffea_path('data/json/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt')
            lumi_mask = LumiMask(json)(df['run'], df['luminosityBlock'])
        else:
            lumi_mask = np.ones(df.size)==1

        # MET filters
        if is_data(dataset):
            filt_met = mask_and(df, cfg.FILTERS.DATA)
        else:
            filt_met = mask_and(df, cfg.FILTERS.MC)


        photons = setup_photons(df)

        ak4 = setup_jets(df)
        ak4 = ak4[
                  object_overlap(ak4, photons) \
                  & ak4.tightId \
                  & (ak4.pt > 100) \
                  & (ak4.abseta < 2.4)
                  ]

        event_mask = filt_met \
                     & lumi_mask \
                     & (ak4.counts > 0) \
                     & df['HLT_Photon200'] \
                     & (df['MET_pt'] < 60)

        # Generator weight
        weights = processor.Weights(size=df.size, storeIndividual=True)

        if is_data(dataset):
            weights.add('gen', np.ones(df.size))
        else:
            weights.add('gen', df['Generator_weight'])

        photon_kinematics = (photons.pt > 200) & (photons.barrel)

        # Medium
        vals = photons[photon_kinematics & photons.mediumId].sieie[event_mask]
        pt = photons[photon_kinematics & photons.mediumId].pt[event_mask]
        output['sieie'].fill(
                             dataset=dataset,
                             cat='medium',
                             sieie=vals.flatten(),
                             pt=pt.flatten(),
                             weights=weight_shape(
                                                  vals,
                                                  weights.weight()[event_mask]
                                                  )
                            )

        # No Sieie
        vals = photons[photon_kinematics & medium_id_no_sieie(photons)].sieie[event_mask]
        pt = photons[photon_kinematics & medium_id_no_sieie(photons)].pt[event_mask]
        output['sieie'].fill(
                             dataset=dataset,
                             cat='medium_nosieie',
                             sieie=vals.flatten(),
                             pt=pt.flatten(),
                             weights=weight_shape(
                                                  vals,
                                                  weights.weight()[event_mask]
                                                  )
                            )

        # No Sieie, inverted isolation
        vals = photons[photon_kinematics & medium_id_no_sieie_inv_iso(photons)].sieie[event_mask]
        pt = photons[photon_kinematics & medium_id_no_sieie_inv_iso(photons)].pt[event_mask]
        output['sieie'].fill(
                             dataset=dataset,
                             cat='medium_nosieie_invertiso',
                             sieie=vals.flatten(),
                             pt=pt.flatten(),
                             weights=weight_shape(
                                                  vals,
                                                  weights.weight()[event_mask]
                                                  )
                            )

        # Keep track of weight sum
        if not is_data(dataset):
            output['sumw'][dataset] +=  df['genEventSumw']
            output['sumw2'][dataset] +=  df['genEventSumw2']
        return output

    def postprocess(self, accumulator):
        return accumulator
