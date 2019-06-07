from coffea import hist
from coffea.analysis_objects import JaggedCandidateArray
import coffea.processor as processor
from awkward import JaggedArray
import numpy as np

import lz4.frame as lz4f
import cloudpickle

class monojetProcessor(processor.ProcessorABC):
    def __init__(self, year="2018"):
        dataset_axis = hist.Cat("dataset", "Primary dataset")
        met_axis = hist.Bin("met", r"$p_{T}^{miss}$ (GeV)", 600, 0.25, 1000)
        jet_pt_axis = hist.Bin("jetpt", r"$p_{T}", 600, 0.25, 1000)

        self._accumulator = processor.dict_accumulator({
            "sr_met" : hist.Hist("Counts", dataset_axis, met_axis),
            "sr_jet0_pt" : hist.Hist("Counts", dataset_axis, jet_pt_axis),
        })

    @property
    def accumulator(self):
        return self._accumulator

    def process(self, df):
        output = self.accumulator.identity()
        dataset = df["dataset"]
        output['sr_met'].fill(dataset=dataset,
                            met=df["MET_pt"])
        output['sr_jet0_pt'].fill(dataset=dataset,
                            jetpt=df["Jet_pt"])
        return output

    def postprocess(self, accumulator):
        return accumulator


def main():
    fileset = {
        'NonthDM' : [
            # 'data/F78663A9-8E7F-B74E-8F6F-9C9A61A27AE5.root'
            "root://cms-xrd-global.cern.ch///store/mc/RunIISummer16NanoAODv4/NonthDMMonoJet_MX-1500_l1-2p_l2-0p04_13TeV-madgraph/NANOAODSIM/PUMoriond17_Nano14Dec2018_102X_mcRun2_asymptotic_v6-v1/260000/F78663A9-8E7F-B74E-8F6F-9C9A61A27AE5.root"

        ]
    }

    output = processor.run_uproot_job(fileset,
                                  treename='Events',
                                  processor_instance=monojetProcessor(),
                                  executor=processor.futures_executor,
                                  executor_args={'workers': 4, 'function_args': {'flatten': True}},
                                  chunksize=500000,
                                 )

    with lz4f.open("hists.cpkl.lz4", mode="wb", compression_level=5) as fout:
        cloudpickle.dump(output, fout)

if __name__ == "__main__":
    main()