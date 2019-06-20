from coffea import hist
import coffea.processor as processor
import os

import matplotlib
# This just tells matplotlib not to open any
# interactive windows.
matplotlib.use('Agg')

class exampleProcessor(processor.ProcessorABC):
    """Dummy processor used to demonstrate the processor principle"""
    def __init__(self):
        dataset_axis = hist.Cat("dataset", "Primary dataset")
        met_axis = hist.Bin("met", r"$p_{T}^{miss}$ (GeV)", 600, 0.25, 1000)
        jet_pt_axis = hist.Bin("jetpt", r"$p_{T}$", 600, 0.25, 1000)
        new_axis = hist.Bin("new_variable", r"Leading jet $p_{T}$ + $p_{T}^{miss}$ (GeV) ", 600, 0.25, 1000)

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

        # And fill the histograms
        output['met'].fill(dataset=dataset,
                            met=df["MET_pt"].flatten())
        output['jet_pt'].fill(dataset=dataset,
                            jetpt=df["Jet_pt"].flatten())

        # We can also do arbitrary transformations
        # E.g.: Sum of MET and the leading jet PTs
        new_variable = df["MET_pt"] + df["Jet_pt"].max()
        output['new_variable'].fill(dataset=dataset,
                            new_variable=new_variable)

        # To apply selections, simply mask
        # Let's see events with MET > 100
        mask = df["MET_pt"] > 100

        # And plot the leading jet pt for these events
        output['jet_pt_met100'].fill(dataset=dataset,
                            jetpt=df["Jet_pt"][mask].max().flatten())

        return output

    def postprocess(self, accumulator):
        return accumulator


def main():
    # Inputs are defined in a dictionary
    # dataset : list of files
    fileset = {
        'NonthDM' : [
            "root://cms-xrd-global.cern.ch///store/mc/RunIISummer16NanoAODv4/NonthDMMonoJet_MX-1500_l1-2p_l2-0p04_13TeV-madgraph/NANOAODSIM/PUMoriond17_Nano14Dec2018_102X_mcRun2_asymptotic_v6-v1/260000/F78663A9-8E7F-B74E-8F6F-9C9A61A27AE5.root"

        ],
        "Znunu_ht600to800" : [
            "root://cms-xrd-global.cern.ch///store/mc/RunIISummer16NanoAODv4/ZJetsToNuNu_HT-600To800_13TeV-madgraph/NANOAODSIM/PUMoriond17_Nano14Dec2018_102X_mcRun2_asymptotic_v6-v1/280000/F4921B81-C2E3-6546-9C00-D908A264FFD8.root",
        ]
    }

    # Run the processor
    output = processor.run_uproot_job(fileset,
                                  treename='Events',
                                  processor_instance=exampleProcessor(),
                                  executor=processor.futures_executor,
                                  executor_args={'workers': 1, 'function_args': {'flatten': False}},
                                  chunksize=500000,
                                 )

    # Make a few plots
    outdir = "./tmp_plots"
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    for name in ["met", "jet_pt", "new_variable","jet_pt_met100"]:
        histogram = output[name]
        fig, ax, _ = hist.plot1d(histogram,overlay="dataset")
        ax.set_yscale('log')
        ax.set_ylim(0.1,1e5)

        fig.savefig(os.path.join(outdir, "{}.pdf".format(name)))


if __name__ == "__main__":
    main()