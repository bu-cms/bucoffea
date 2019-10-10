import os
import re

import coffea.processor as processor
import numpy as np
from awkward import JaggedArray
from coffea import hist
from coffea.analysis_objects import JaggedCandidateArray

from bucoffea.helpers.dataset import (extract_year, is_lo_w, is_lo_z, is_nlo_w,
                                      is_nlo_z)
from bucoffea.helpers.gen import find_gen_dilepton, setup_gen_candidates

from bucoffea.helpers import min_dphi_jet_met

Hist = hist.Hist
Bin = hist.Bin
Cat = hist.Cat

def vbf_selection(dilepton, dijet, genjets):
    selection = processor.PackedSelection()

    selection.add(
                  'two_jets',
                  dijet.counts>0
                  )
    selection.add(
                  'leadak4_pt_eta',
                  (dijet.i0.pt.max() > 80) & (np.abs(dijet.i0.eta.max()) < 4.7)
                  )
    selection.add(
                  'trailak4_pt_eta',
                  (dijet.i1.pt.max() > 40) & (np.abs(dijet.i1.eta.max()) < 4.7)
                  )
    selection.add(
                  'hemisphere',
                  (dijet.i0.eta.max()*dijet.i1.eta.max() < 0)
                  )
    selection.add(
                  'mindphijr',
                  min_dphi_jet_met(genjets, dilepton.phi.max(), njet=4, ptmin=30) > 0.5
                  )

    return selection

def monojet_selection(dilepton, genjets):
    selection = processor.PackedSelection()

    selection.add(
                  'at_least_one_jet',
                  genjets.counts>0
                  )
    selection.add(
                  'leadak4_pt_eta',
                  (genjets.pt.max() > 100) & (np.abs(genjets[genjets.pt.argmax()].eta.max()) < 2.4)
                  )
    selection.add(
                  'mindphijr',
                  min_dphi_jet_met(genjets, dilepton.phi.max(), njet=4, ptmin=30) > 0.5
                  )

    return selection

class lheVProcessor(processor.ProcessorABC):
    def __init__(self):

        # Histogram setup
        dataset_ax = Cat("dataset", "Primary dataset")

        vpt_ax = Bin("vpt",r"$p_{T}^{V}$ (GeV)", 50, 0, 2000)
        jpt_ax = Bin("jpt",r"$p_{T}^{j}$ (GeV)", 50, 0, 2000)
        mjj_ax = Bin("mjj",r"$m(jj)$ (GeV)", 75, 0, 7500)

        items = {}
        items["gen_vpt_inclusive"] = Hist("Counts",
                                dataset_ax,
                                vpt_ax)
        items["gen_vpt_monojet"] = Hist("Counts",
                                dataset_ax,
                                jpt_ax,
                                vpt_ax)
        items["gen_vpt_vbf"] = Hist("Counts",
                                dataset_ax,
                                jpt_ax,
                                mjj_ax,
                                vpt_ax)
        items['sumw'] = processor.defaultdict_accumulator(float)
        items['sumw2'] = processor.defaultdict_accumulator(float)

        self._accumulator = processor.dict_accumulator(items)

    @property
    def accumulator(self):
        return self._accumulator


    def process(self, df):
        output = self.accumulator.identity()
        dataset = df['dataset']

        # Dilepton
        gen = setup_gen_candidates(df)

        genjets = JaggedCandidateArray.candidatesfromcounts(
                df['nGenJet'],
                pt=df['GenJet_pt'],
                eta=df['GenJet_eta'],
                phi=df['GenJet_phi'],
                mass=df['GenJet_mass']
            )

        if is_lo_z(dataset) or is_nlo_z(dataset):
            pdgsum = 0
            # gen_v = gen[(gen.status==62) & (gen.pdg==23)]
        elif is_lo_w(dataset) or is_nlo_w(dataset):
            pdgsum = 1
            # gen_v = gen[(gen.status==62) & (gen.pdg==24)]
        gen_dilep = find_gen_dilepton(gen, pdgsum)
        gen_dilep = gen_dilep[gen_dilep.mass.argmax()]

        # Select jets not overlapping with leptons
        genjets = genjets[
            (~genjets.match(gen_dilep.i0,deltaRCut=0.4)) &
            (~genjets.match(gen_dilep.i1,deltaRCut=0.4)) \
            ]

        # Dijet for VBF
        dijet = genjets[:,:2].distincts()

        # Selection
        vbf_sel = vbf_selection(gen_dilep, dijet, genjets)
        monojet_sel = monojet_selection(gen_dilep, genjets)

        nominal = df['Generator_weight']

        output['gen_vpt_inclusive'].fill(
                                dataset=dataset,
                                vpt=gen_dilep.pt.max(),
                                jpt=genjets.pt.max(),
                                weight=nominal
                                )

        mask_vbf = vbf_sel.all(*vbf_sel.names)
        print(mask_vbf.size, mask_vbf)
        print(nominal.size, nominal)
        output['gen_vpt_vbf'].fill(
                                dataset=dataset,
                                vpt=gen_dilep.pt.max()[mask_vbf],
                                jpt=genjets.pt.max()[mask_vbf],
                                mjj = dijet.mass.max()[mask_vbf],
                                weight=nominal[mask_vbf]
                                )

        mask_monojet = monojet_sel.all(*monojet_sel.names)

        print(mask_monojet.size, mask_monojet)
        print(nominal.size, nominal)
        output['gen_vpt_monojet'].fill(
                                dataset=dataset,
                                vpt=gen_dilep.pt.max()[mask_monojet],
                                jpt=genjets.pt.max()[mask_monojet],
                                weight=nominal[mask_monojet]
                                )

        # Keep track of weight sum
        output['sumw'][dataset] +=  df['genEventSumw']
        output['sumw2'][dataset] +=  df['genEventSumw2']
        return output

    def postprocess(self, accumulator):
        return accumulator
