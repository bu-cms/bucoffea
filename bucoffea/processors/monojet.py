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

        ],
        "Znunu_ht600to800" : [
            "/store/mc/RunIISummer16NanoAODv4/ZJetsToNuNu_HT-600To800_13TeV-madgraph/NANOAODSIM/PUMoriond17_Nano14Dec2018_102X_mcRun2_asymptotic_v6-v1/280000/F4921B81-C2E3-6546-9C00-D908A264FFD8.root",
            # "/store/mc/RunIISummer16NanoAODv4/ZJetsToNuNu_HT-600To800_13TeV-madgraph/NANOAODSIM/PUMoriond17_Nano14Dec2018_102X_mcRun2_asymptotic_v6-v1/280000/67F57995-381B-B14B-939B-23B6103564C9.root",
            # "/store/mc/RunIISummer16NanoAODv4/ZJetsToNuNu_HT-600To800_13TeV-madgraph/NANOAODSIM/PUMoriond17_Nano14Dec2018_102X_mcRun2_asymptotic_v6-v1/280000/6496BE8A-9967-7641-98FE-730C54E69A7E.root",
            # "/store/mc/RunIISummer16NanoAODv4/ZJetsToNuNu_HT-600To800_13TeV-madgraph/NANOAODSIM/PUMoriond17_Nano14Dec2018_102X_mcRun2_asymptotic_v6-v1/280000/4F02C217-D5AF-0443-A407-CCEEF8388273.root",
            # "/store/mc/RunIISummer16NanoAODv4/ZJetsToNuNu_HT-600To800_13TeV-madgraph/NANOAODSIM/PUMoriond17_Nano14Dec2018_102X_mcRun2_asymptotic_v6-v1/280000/365C896D-252D-8D4D-B350-48458EBA23D3.root",
            # "/store/mc/RunIISummer16NanoAODv4/ZJetsToNuNu_HT-600To800_13TeV-madgraph/NANOAODSIM/PUMoriond17_Nano14Dec2018_102X_mcRun2_asymptotic_v6-v1/280000/06FFDDCB-394D-DB40-B59E-DC0DAAAA4C5F.root",
        ],
        "Znunu_ht200to400" : [
            "/store/mc/RunIISummer16NanoAODv4/ZJetsToNuNu_HT-200To400_13TeV-madgraph/NANOAODSIM/PUMoriond17_Nano14Dec2018_102X_mcRun2_asymptotic_v6_ext1-v1/40000/FFD69E5A-A941-2D41-A629-9D62D9E8BE9A.root",
            # "/store/mc/RunIISummer16NanoAODv4/ZJetsToNuNu_HT-200To400_13TeV-madgraph/NANOAODSIM/PUMoriond17_Nano14Dec2018_102X_mcRun2_asymptotic_v6_ext1-v1/40000/FE9DACDF-9EDF-3349-8C70-0F9BA0E61FC5.root",
            # "/store/mc/RunIISummer16NanoAODv4/ZJetsToNuNu_HT-200To400_13TeV-madgraph/NANOAODSIM/PUMoriond17_Nano14Dec2018_102X_mcRun2_asymptotic_v6_ext1-v1/40000/F9D59028-EB30-A144-BEB4-0BDA742A0B14.root",
            # "/store/mc/RunIISummer16NanoAODv4/ZJetsToNuNu_HT-200To400_13TeV-madgraph/NANOAODSIM/PUMoriond17_Nano14Dec2018_102X_mcRun2_asymptotic_v6_ext1-v1/40000/D53E6946-5023-4141-B97A-4A25F29B78E5.root",
            # "/store/mc/RunIISummer16NanoAODv4/ZJetsToNuNu_HT-200To400_13TeV-madgraph/NANOAODSIM/PUMoriond17_Nano14Dec2018_102X_mcRun2_asymptotic_v6_ext1-v1/40000/BFA74576-4D2D-8543-B94D-271EA05515AB.root",
            # "/store/mc/RunIISummer16NanoAODv4/ZJetsToNuNu_HT-200To400_13TeV-madgraph/NANOAODSIM/PUMoriond17_Nano14Dec2018_102X_mcRun2_asymptotic_v6_ext1-v1/40000/B59FF6B2-0116-504C-A118-9D953AA9DE69.root",
            # "/store/mc/RunIISummer16NanoAODv4/ZJetsToNuNu_HT-200To400_13TeV-madgraph/NANOAODSIM/PUMoriond17_Nano14Dec2018_102X_mcRun2_asymptotic_v6_ext1-v1/40000/AE997E4F-8745-2B42-874C-24BEB7A1A5DD.root",
            # "/store/mc/RunIISummer16NanoAODv4/ZJetsToNuNu_HT-200To400_13TeV-madgraph/NANOAODSIM/PUMoriond17_Nano14Dec2018_102X_mcRun2_asymptotic_v6_ext1-v1/40000/8876546B-7CB6-044F-BE53-1938AC50FE78.root",
            # "/store/mc/RunIISummer16NanoAODv4/ZJetsToNuNu_HT-200To400_13TeV-madgraph/NANOAODSIM/PUMoriond17_Nano14Dec2018_102X_mcRun2_asymptotic_v6_ext1-v1/40000/85EDD994-32B6-3C41-A3C8-540EBE7B7083.root",
            # "/store/mc/RunIISummer16NanoAODv4/ZJetsToNuNu_HT-200To400_13TeV-madgraph/NANOAODSIM/PUMoriond17_Nano14Dec2018_102X_mcRun2_asymptotic_v6_ext1-v1/40000/65EB082B-22A8-7B44-A931-6713EC695892.root",
            # "/store/mc/RunIISummer16NanoAODv4/ZJetsToNuNu_HT-200To400_13TeV-madgraph/NANOAODSIM/PUMoriond17_Nano14Dec2018_102X_mcRun2_asymptotic_v6_ext1-v1/40000/58FCF07C-A3AA-4847-9D63-C1E82F2C25F4.root",
            # "/store/mc/RunIISummer16NanoAODv4/ZJetsToNuNu_HT-200To400_13TeV-madgraph/NANOAODSIM/PUMoriond17_Nano14Dec2018_102X_mcRun2_asymptotic_v6_ext1-v1/40000/58EAC30F-CE51-7748-9916-487CDACC01AA.root",
            # "/store/mc/RunIISummer16NanoAODv4/ZJetsToNuNu_HT-200To400_13TeV-madgraph/NANOAODSIM/PUMoriond17_Nano14Dec2018_102X_mcRun2_asymptotic_v6_ext1-v1/40000/54A7E4EB-2388-2248-9E1C-E50D5A1EF99D.root",
            # "/store/mc/RunIISummer16NanoAODv4/ZJetsToNuNu_HT-200To400_13TeV-madgraph/NANOAODSIM/PUMoriond17_Nano14Dec2018_102X_mcRun2_asymptotic_v6_ext1-v1/40000/43BDF497-350C-554B-A244-321183DAE9D0.root",
            # "/store/mc/RunIISummer16NanoAODv4/ZJetsToNuNu_HT-200To400_13TeV-madgraph/NANOAODSIM/PUMoriond17_Nano14Dec2018_102X_mcRun2_asymptotic_v6_ext1-v1/40000/0D1DB570-DDCF-BE4A-B523-A736ADD4F702.root",
            # "/store/mc/RunIISummer16NanoAODv4/ZJetsToNuNu_HT-200To400_13TeV-madgraph/NANOAODSIM/PUMoriond17_Nano14Dec2018_102X_mcRun2_asymptotic_v6_ext1-v1/40000/0A6D78C7-5E74-3142-AD42-FC3FB90937F3.root",
        ]
    }

    for dataset, filelist in fileset.items():
        newlist = []
        for file in filelist:
            if file.startswith("/store/"):
                newlist.append("root://cms-xrd-global.cern.ch//" + file)
            else: newlist.append(file)
        fileset[dataset] = newlist

    output = processor.run_uproot_job(fileset,
                                  treename='Events',
                                  processor_instance=monojetProcessor(),
                                  executor=processor.futures_executor,
                                  executor_args={'workers': 1, 'function_args': {'flatten': True}},
                                  chunksize=500000,
                                 )

    with lz4f.open("hists.cpkl.lz4", mode="wb", compression_level=5) as fout:
        cloudpickle.dump(output, fout)

if __name__ == "__main__":
    main()