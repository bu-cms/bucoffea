import os
import numpy as np
import re
from coffea import hist
import coffea.processor as processor


#  export LHAPDF_DATA_PATH=/cvmfs/sft.cern.ch/lcg/releases/MCGenerators/lhapdf/6.2.1-7149a/x86_64-centos7-gcc8-opt/share/LHAPDF/
def setup_pdf_objs(all_pdfs):
    import lhapdf
    lhapdf.setVerbosity(1)
    os.environ['LHAPDF_DATA_PATH'] = '/cvmfs/sft.cern.ch/lcg/releases/MCGenerators/lhapdf/6.2.1-7149a/x86_64-centos7-gcc8-opt/share/LHAPDF/'

    paths = [
    '/cvmfs/sft.cern.ch/lcg/external/lhapdfsets/current/',
    ]
    lhapdf.setPaths(paths)

    pdf_objs = {}
    for p in all_pdfs:
        pdf_objs[p] = lhapdf.mkPDF(p)
    return pdf_objs

Hist = hist.Hist
Bin = hist.Bin
Cat = hist.Cat




class pdfWeightProcessor(processor.ProcessorABC):
    def __init__(self):
        dataset_ax = Cat("dataset", "Primary dataset")

        # Histogram setup
        pdf_ax = Cat("pdf", "pdf")
        vpt_ax = Bin("vpt",r"$p_{T}^{V}$ (GeV)", 100, 0, 2000)
        items = {}
        items["gen_vpt"] = Hist("Counts", dataset_ax, vpt_ax, vpt_ax)
        items['sumw'] = processor.defaultdict_accumulator(float)
        items['sumw2'] = processor.defaultdict_accumulator(float)
        self._accumulator = processor.dict_accumulator(items)

    @property
    def accumulator(self):
        return self._accumulator


    def process(self, df):
        output = self.accumulator.identity()
        dataset = df['dataset']

        all_pdfs = [303600,263000]
        pdf_objs = setup_pdf_objs(all_pdfs)
        base_pdf = 303600

        pdf_weights = {}

        for ipdf in all_pdfs:
            w = np.ones(df.size)
            for i in range(df.size):
                id1 = df['Generator_id1'][i]
                id2 = df['Generator_id2'][i]
                x1 = df['Generator_x1'][i]
                x2 = df['Generator_x2'][i]
                q2 = df['Generator_scalePDF'][i]

                w[i] = pdf_objs[ipdf].xfxQ2(id1, x1, q2) * pdf_objs[ipdf].xfxQ2(id2, x2, q2)
            pdf_weights[ipdf] = w


        for ipdf in pdf_weights.keys():
            output['gen_vpt'].fill(
                                vpt = df['LHE_Vpt'],
                                dataset=dataset,
                                pdf=str(ipdf),
                                weight=pdf_weights[ipdf]/pdf_weights[base_pdf]
                                )
        output['sumw'][dataset] +=  df['genEventSumw']
        output['sumw2'][dataset] +=  df['genEventSumw2']
        return output

    def postprocess(self, accumulator):
        return accumulator



if __name__ == "__main__":
    main()