import os
import re

import coffea.processor as processor
import numpy as np
from awkward import JaggedArray
from coffea import hist

from bucoffea.helpers.dataset import (extract_year, is_lo_w, is_lo_z, is_nlo_w,
                                      is_nlo_z)
from bucoffea.helpers.gen import find_gen_dilepton, setup_gen_candidates

Hist = hist.Hist
Bin = hist.Bin
Cat = hist.Cat


class lheVProcessor(processor.ProcessorABC):
    def __init__(self):

        # Histogram setup
        dataset_ax = Cat("dataset", "Primary dataset")
        type_ax = Cat("type", "Calculation type")
        weight_type_ax = Cat("weight_type", "Weight type")
        weight_index_ax = Bin("weight_index",r"weight index", 50,-0.5,49.5)

        vpt_ax = Bin("vpt",r"$p_{T}^{V}$ (GeV)", 100, 0, 2000)

        items = {}
        items["gen_vpt"] = Hist("Counts",
                                dataset_ax,
                                weight_type_ax,
                                weight_index_ax,
                                type_ax,
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
        if is_lo_z(dataset) or is_nlo_z(dataset):
            pdgsum = 0
            gen_v = gen[(gen.status==62) & (gen.pdg==23)]
        elif is_lo_w(dataset) or is_nlo_w(dataset):
            pdgsum = 1
            gen_v = gen[(gen.status==62) & (gen.pdg==24)]
        gen_dilep = find_gen_dilepton(gen, pdgsum)


        nominal = df['Generator_weight']

        # Fill
        output['gen_vpt'].fill(
                                dataset=dataset,
                                vpt=df['LHE_Vpt'],
                                weight_type='nominal',
                                weight_index=0,
                                weight=nominal,
                                type='nano'
                                )

        output['gen_vpt'].fill(
                                dataset=dataset,
                                vpt=gen_dilep[gen_dilep.mass.argmax()].pt.max(),
                                weight_type='nominal',
                                weight_index=0,
                                weight=nominal,
                                type='dilepton'
                                )
                            
        output['gen_vpt'].fill(
                                dataset=dataset,
                                vpt=gen_v.pt.max(),
                                weight_type='nominal',
                                weight_index=0,
                                weight=nominal,
                                type='genv'
                                )

        # PDF variations
        w_pdf = JaggedArray.fromcounts(df['nLHEPdfWeight'],df['LHEPdfWeight'])
        for i in range(df['nLHEPdfWeight'][0]):
            output['gen_vpt'].fill(
                                dataset=dataset,
                                vpt=df['LHE_Vpt'],
                                weight_type='pdf',
                                weight_index=i,
                                weight=nominal*w_pdf[:,i],
                                type='nano'
                                )
            output['gen_vpt'].fill(
                                dataset=dataset,
                                vpt=gen_dilep.pt.max(),
                                weight_type='pdf',
                                weight_index=i,
                                weight=nominal*w_pdf[:,i],
                                type='dilepton'
                                )

        # Scale variations
        w_scale = JaggedArray.fromcounts(df['nLHEScaleWeight'],df['LHEScaleWeight'])
        for i in range(df['nLHEScaleWeight'][0]):
            output['gen_vpt'].fill(
                                dataset=dataset,
                                vpt=df['LHE_Vpt'],
                                weight_type='scale',
                                weight_index=i,
                                weight=nominal*w_scale[:,i],
                                type='nano'
                                )
            output['gen_vpt'].fill(
                                dataset=dataset,
                                vpt=gen_dilep.pt.max(),
                                weight_type='scale',
                                weight_index=i,
                                weight=nominal*w_pdf[:,i],
                                type='dilepton'
                                )

        # Keep track of weight sum
        output['sumw'][dataset] +=  df['genEventSumw']
        output['sumw2'][dataset] +=  df['genEventSumw2']
        return output

    def postprocess(self, accumulator):
        return accumulator
