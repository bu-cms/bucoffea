import os
import numpy as np
import re
from coffea import hist
import coffea.processor as processor


#  export LHAPDF_DATA_PATH=/cvmfs/sft.cern.ch/lcg/releases/MCGenerators/lhapdf/6.2.1-7149a/x86_64-centos7-gcc8-opt/share/LHAPDF/
def setup_pdf_objs(all_pdfs, cache={}):
    import lhapdf
    lhapdf.setVerbosity(1)
    os.environ['LHAPDF_DATA_PATH'] = '/cvmfs/sft.cern.ch/lcg/releases/MCGenerators/lhapdf/6.2.1-7149a/x86_64-centos7-gcc8-opt/share/LHAPDF/'

    paths = [
    '/cvmfs/sft.cern.ch/lcg/external/lhapdfsets/current/',
    ]
    lhapdf.setPaths(paths)

    pdf_objs = {}
    for p in all_pdfs:
        if p in cache:
            pdf_objs[p] = cache[p]
        else:
            pdf_objs[p] = lhapdf.mkPDF(p)
            cache[p] = pdf_objs[p]
    return pdf_objs

Hist = hist.Hist
Bin = hist.Bin
Cat = hist.Cat




class pdfWeightProcessor(processor.ProcessorABC):
    def __init__(self):

        # Histogram setup
        dataset_ax = Cat("dataset", "Primary dataset")
        pdf_ax = Cat("pdf", "pdf")
        vpt_ax = Bin("vpt",r"$p_{T}^{V}$ (GeV)", 100, 0, 2000)
        sign_ax = Bin("sign",r"Sign", 2, -2, 2)

        items = {}
        items["gen_vpt"] = Hist("Counts", dataset_ax, pdf_ax, vpt_ax)
        items["gen_weight_sign"] = Hist("Counts", dataset_ax, sign_ax)
        items['sumw'] = processor.defaultdict_accumulator(float)
        items['sumw2'] = processor.defaultdict_accumulator(float)
        self._accumulator = processor.dict_accumulator(items)

        self._all_pdfs = [303600,263000,262000,306000]

    @property
    def accumulator(self):
        return self._accumulator


    def process(self, df):
        output = self.accumulator.identity()
        dataset = df['dataset']


        base_pdf = 306000

        pdf_weights = {
            'none' : np.ones(df.size)
        }

        pdf_objs = setup_pdf_objs(self._all_pdfs)
        for ipdf in self._all_pdfs:
            w = np.ones(df.size)
            for i in range(df.size):
                id1 = df['Generator_id1'][i]
                id2 = df['Generator_id2'][i]
                x1 = df['Generator_x1'][i]
                x2 = df['Generator_x2'][i]
                q2 = df['Generator_scalePDF'][i]

                w[i] = pdf_objs[ipdf].xfxQ2(id1, x1, q2) * pdf_objs[ipdf].xfxQ2(id2, x2, q2)
                # print(id1, id2, x1, x2, q, w[i])
            pdf_weights[ipdf] = w


        output["gen_weight_sign"].fill(
            sign=np.sign(df['Generator_weight']),
            dataset=dataset
        )
        for ipdf in pdf_weights.keys():
            weight = df['Generator_weight']
            if ipdf !='none':
                reweight = pdf_weights[ipdf]/pdf_weights[base_pdf]
                reweight[abs(reweight)>10] = 1
                weight = weight * reweight
            output['gen_vpt'].fill(
                                vpt = df['LHE_Vpt'],
                                dataset=dataset,
                                pdf=str(ipdf),
                                weight=weight
                                )
        output['sumw'][dataset] +=  df['genEventSumw']
        output['sumw2'][dataset] +=  df['genEventSumw2']
        return output

    def postprocess(self, accumulator):
        return accumulator