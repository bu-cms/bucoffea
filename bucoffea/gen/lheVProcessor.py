import os
import numpy as np
import re
from coffea import hist
import coffea.processor as processor
from awkward import JaggedArray

Hist = hist.Hist
Bin = hist.Bin
Cat = hist.Cat


class lheVProcessor(processor.ProcessorABC):
    def __init__(self):

        # Histogram setup
        dataset_ax = Cat("dataset", "Primary dataset")
        weight_type_ax = Cat("weight_type", "Weight type")
        weight_index_ax = Bin("weight_index",r"weight index", 50,-0.5,49.5)

        vpt_ax = Bin("vpt",r"$p_{T}^{V}$ (GeV)", 100, 0, 2000)

        items = {}
        items["gen_vpt"] = Hist("Counts",
                                dataset_ax,
                                weight_type_ax,
                                weight_index_ax,
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

        nominal = df['Generator_weight']

        output['gen_vpt'].fill(
                                dataset=dataset,
                                vpt=df['LHE_Vpt'],
                                weight_type='nominal',
                                weight_index=0,
                                weight=nominal,
                                )

        w_pdf = JaggedArray.fromcounts(df['nLHEPdfWeight'],df['LHEPdfWeight'])
        for i in range(df['nLHEPdfWeight'][0]):
            output['gen_vpt'].fill(
                                dataset=dataset,
                                vpt=df['LHE_Vpt'],
                                weight_type='pdf',
                                weight_index=i,
                                weight=nominal*w_pdf[:,i]
                                )

        w_scale = JaggedArray.fromcounts(df['nLHEScaleWeight'],df['LHEScaleWeight'])
        for i in range(df['nLHEScaleWeight'][0]):
            output['gen_vpt'].fill(
                                dataset=dataset,
                                vpt=df['LHE_Vpt'],
                                weight_type='scale',
                                weight_index=i,
                                weight=nominal*w_scale[:,i]
                                )

        # Keep track of weight sum
        output['sumw'][dataset] +=  df['genEventSumw']
        output['sumw2'][dataset] +=  df['genEventSumw2']
        return output

    def postprocess(self, accumulator):
        return accumulator



if __name__ == "__main__":
    main()