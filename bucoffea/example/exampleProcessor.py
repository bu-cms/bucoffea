#!/usr/bin/env python

from coffea import hist
import coffea.processor as processor
from coffea.analysis_objects import JaggedCandidateArray
import os

from matplotlib import pyplot as plt
import matplotlib
import numpy as np
# This just tells matplotlib not to open any
# interactive windows.
matplotlib.use('Agg')

class exampleProcessor(processor.ProcessorABC):
    """Dummy processor used to demonstrate the processor principle"""
    def __init__(self):
        dataset_axis = hist.Cat("dataset", "Primary dataset")
        met_axis = hist.Bin("met", r"$p_{T}^{miss}$ (GeV)", 100, 0, 1000)
        jet_pt_axis = hist.Bin("jetpt", r"$p_{T}$", 100, 0, 1000)
        new_axis = hist.Bin("new_variable", r"Leading jet $p_{T}$ + $p_{T}^{miss}$ (GeV) ", 100, 0, 1000)

        self._accumulator = processor.dict_accumulator({
            "met" : hist.Hist("Counts", dataset_axis, met_axis),
            "jet_pt" : hist.Hist("Counts", dataset_axis, jet_pt_axis),
            "jet_pt_met100" : hist.Hist("Counts", dataset_axis, jet_pt_axis),
            "new_variable" : hist.Hist("Counts", dataset_axis, new_axis),
        })

    @property
    def accumulator(self):
        return self._accumulator

    def process(self, df):
        """
        Processing function. This is where the actual analysis happens.
        """
        output = self.accumulator.identity()
        # We can access the data frame as usual
        # The dataset is written into the data frame
        # outside of this function
        dataset = df["dataset"]

        # JaggedCandidateArray lets us combine
        # multiple branches into an array of
        # particle candidates to make our life easier
        jets = JaggedCandidateArray.candidatesfromcounts(
            df['nJet'],
            pt=df[f'Jet_pt'],
            eta=df['Jet_eta'],
            phi=df['Jet_phi'],
            mass=np.zeros_like(df['Jet_pt']),
        )

        # Fill the histograms
        output['met'].fill(dataset=dataset,
                            met=df["MET_pt"].flatten())
        output['jet_pt'].fill(dataset=dataset,
                            jetpt=jets.pt.flatten())

        # We can also do arbitrary transformations
        # E.g.: Sum of MET and the leading jet PTs
        new_variable = df["MET_pt"] + jets.pt.max()
        output['new_variable'].fill(dataset=dataset,
                            new_variable=new_variable)

        # To apply selections, simply mask
        # Let's see events with MET > 100
        mask = df["MET_pt"] > 100

        # And plot the leading jet pt for these events
        output['jet_pt_met100'].fill(dataset=dataset,
                            jetpt=jets.pt[mask].max().flatten())

        return output

    def postprocess(self, accumulator):
        return accumulator


def main():
    # Inputs are defined in a dictionary
    # dataset : list of files
    fileset = {
        "ADD" : [
            "root://cms-xrd-global.cern.ch///store/mc/RunIIFall17NanoAODv7/ADDMonoJet_MD_8_d_7_TuneCUETP8M1_13TeV_pythia8/NANOAODSIM/PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/230000/0202BA5A-A595-3C46-B1EB-EDF35E2BD642.root"
        ],
        # "Zjet" : [
        #     "root://cms-xrd-global.cern.ch///store/mc/RunIIFall17NanoAODv7/ZJetsToNuNu_HT-400To600_13TeV-madgraph/NANOAODSIM/PU2017_12Apr2018_Nano02Apr2020_new_pmx_102X_mc2017_realistic_v8-v1/100000/50B7CBF1-0B29-C846-AD85-85259D438DA7.root"
        # ]
    }

    # Run the processor
    output = processor.run_uproot_job(fileset,
                                  treename='Events',
                                  processor_instance=exampleProcessor(),
                                  executor=processor.futures_executor,
                                  executor_args={'workers': 1, 'flatten': True},
                                  chunksize=500000,
                                 )

    # Make a few plots
    outdir = "./tmp_plots"
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    for name in ["met", "jet_pt", "new_variable","jet_pt_met100"]:
        histogram = output[name]
        fig, ax = plt.subplots()
        hist.plot1d(histogram, ax=ax, overlay="dataset")
        ax.set_yscale('log')
        ax.set_ylim(0.1,1e5)

        fig.savefig(os.path.join(outdir, "{}.pdf".format(name)))


if __name__ == "__main__":
    main()