import os
import re

import coffea.processor as processor
import numpy as np
from awkward import JaggedArray
from coffea import hist
from coffea.analysis_objects import JaggedCandidateArray

from bucoffea.helpers.dataset import (extract_year, is_lo_w, is_lo_z, is_nlo_w,
                                      is_nlo_z)
from bucoffea.helpers.gen import find_gen_dilepton, setup_gen_candidates, setup_dressed_gen_candidates, isnu, islep, fill_gen_v_info

from bucoffea.helpers import min_dphi_jet_met

Hist = hist.Hist
Bin = hist.Bin
Cat = hist.Cat

def vbf_selection(vphi, dijet, genjets):
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
                  min_dphi_jet_met(genjets, vphi.max(), njet=4, ptmin=30) > 0.5
                  )

    return selection

def monojet_selection(vphi, genjets):
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
                  min_dphi_jet_met(genjets, vphi.max(), njet=4, ptmin=30) > 0.5
                  )

    return selection


class lheVProcessor(processor.ProcessorABC):
    def __init__(self):

        # Histogram setup
        dataset_ax = Cat("dataset", "Primary dataset")

        vpt_ax = Bin("vpt",r"$p_{T}^{V}$ (GeV)", 50, 0, 2000)
        jpt_ax = Bin("jpt",r"$p_{T}^{j}$ (GeV)", 50, 0, 2000)
        mjj_ax = Bin("mjj",r"$m(jj)$ (GeV)", 75, 0, 7500)
        res_ax = Bin("res",r"pt: dressed / stat1 - 1", 80,-0.2,0.2)

        items = {}
        for tag in ['stat1','dress','lhe']:
            items[f"gen_vpt_inclusive_{tag}"] = Hist("Counts",
                                    dataset_ax,
                                    vpt_ax)
            items[f"gen_vpt_monojet_{tag}"] = Hist("Counts",
                                    dataset_ax,
                                    jpt_ax,
                                    vpt_ax)
            items[f"gen_vpt_vbf_{tag}"] = Hist("Counts",
                                    dataset_ax,
                                    jpt_ax,
                                    mjj_ax,
                                    vpt_ax)
        items["resolution"] = Hist("Counts",
                                dataset_ax,
                                res_ax)
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
        

        genjets_all = JaggedCandidateArray.candidatesfromcounts(
                df['nGenJet'],
                pt=df['GenJet_pt'],
                eta=df['GenJet_eta'],
                abseta=np.abs(df['GenJet_eta']),
                phi=df['GenJet_phi'],
                mass=df['GenJet_mass']
            )
        gen = setup_gen_candidates(df)
        dressed = setup_dressed_gen_candidates(df)
        fill_gen_v_info(df, gen, dressed)

        for tag in ['lhe','dress','stat1']:
            # Select jets not overlapping with leptons
            genjets = genjets_all[
            (~genjets_all.match(dressed,deltaRCut=0.4)) &
            (~genjets_all.match(gen[islep(gen)],deltaRCut=0.4)) \
            ]

            # Dijet for VBF
            dijet = genjets[:,:2].distincts()

            # Selection
            vbf_sel = vbf_selection(df[f'gen_v_phi_{tag}'], dijet, genjets)
            monojet_sel = monojet_selection(df[f'gen_v_phi_{tag}'], genjets)

            nominal = df['Generator_weight']

            output[f'gen_vpt_inclusive_{tag}'].fill(
                                    dataset=dataset,
                                    vpt=df[f'gen_v_pt_{tag}'],
                                    jpt=genjets.pt.max(),
                                    weight=nominal
                                    )

            mask_vbf = vbf_sel.all(*vbf_sel.names)
            output[f'gen_vpt_vbf_{tag}'].fill(
                                    dataset=dataset,
                                    vpt=df[f'gen_v_pt_{tag}'][mask_vbf],
                                    jpt=genjets.pt.max()[mask_vbf],
                                    mjj = dijet.mass.max()[mask_vbf],
                                    weight=nominal[mask_vbf]
                                    )

            mask_monojet = monojet_sel.all(*monojet_sel.names)

            output[f'gen_vpt_monojet_{tag}'].fill(
                                    dataset=dataset,
                                    vpt=df[f'gen_v_pt_{tag}'][mask_monojet],
                                    jpt=genjets.pt.max()[mask_monojet],
                                    weight=nominal[mask_monojet]
                                    )

        output['resolution'].fill(
            dataset=dataset,
            res = df['gen_v_pt_dress'] / df['gen_v_pt_stat1'] - 1
        )

        # Keep track of weight sum
        output['sumw'][dataset] +=  df['genEventSumw']
        output['sumw2'][dataset] +=  df['genEventSumw2']
        return output

    def postprocess(self, accumulator):
        return accumulator
