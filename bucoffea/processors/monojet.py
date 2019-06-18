from coffea import hist
from coffea.analysis_objects import JaggedCandidateArray
import coffea.processor as processor
from awkward import JaggedArray
import numpy as np

import lz4.frame as lz4f
import cloudpickle

import os
pjoin = os.path.join
def setup_candidates(df):
    muons = JaggedCandidateArray.candidatesfromcounts(
        df['nMuon'].flatten(),
        pt=df['Muon_pt'].flatten(),
        eta=df['Muon_eta'].flatten(),
        phi=df['Muon_phi'].flatten(),
        mass=df['Muon_mass'].flatten(),
        charge=df['Muon_charge'].flatten(),
        mediumId=df['Muon_mediumId'].flatten(),
        tightId=df['Muon_tightId'].flatten()
    )

    electrons = JaggedCandidateArray.candidatesfromcounts(
        df['nElectron'].flatten(),
        pt=df['Electron_pt'].flatten(),
        eta=df['Electron_eta'].flatten(),
        phi=df['Electron_phi'].flatten(),
        mass=df['Electron_mass'].flatten(),
        charge=df['Electron_charge'].flatten(),
        looseId=(df['Electron_cutBased_Sum16']>=1).flatten(),
        tightId=(df['Electron_cutBased_Sum16']==4).flatten()
    )
    jets = JaggedCandidateArray.candidatesfromcounts(
        df['nJet'].flatten(),
        pt=df['Jet_pt'].flatten(),
        eta=df['Jet_eta'].flatten(),
        phi=df['Jet_phi'].flatten(),
        mass=df['Jet_mass'].flatten(),
        nef=df['Jet_neEmEF'].flatten(),
        nhf=df['Jet_neHEF'].flatten(),
        chf=df['Jet_chHEF'].flatten(),
        cef=df['Jet_chEmEF'].flatten(),
    )
    return jets, muons, electrons

class monojetProcessor(processor.ProcessorABC):
    def __init__(self, year="2018"):
        dataset_axis = hist.Cat("dataset", "Primary dataset")
        met_axis = hist.Bin("met", r"$p_{T}^{miss}$ (GeV)", 100, 0, 1000)
        jet_pt_axis = hist.Bin("jetpt", r"$p_{T}$ (GeV)", 100, 0, 1000)
        jet_eta_axis = hist.Bin("jeteta", r"$\eta$ (GeV)", 50, -5, 5)
        dpfcalo_axis = hist.Bin("dpfcalo", r"$1-Calo/PF$", 20, -1, 1)

        self._accumulator = processor.dict_accumulator({
            "met" : hist.Hist("Counts", dataset_axis, met_axis),
            "jet0pt" : hist.Hist("Counts", dataset_axis, jet_pt_axis),
            "jet0eta" : hist.Hist("Counts", dataset_axis, jet_eta_axis),
            "jetpt" : hist.Hist("Counts", dataset_axis, jet_pt_axis),
            "jeteta" : hist.Hist("Counts", dataset_axis, jet_eta_axis),
            "dpfcalo" : hist.Hist("Counts", dataset_axis, dpfcalo_axis),
        })

    @property
    def accumulator(self):
        return self._accumulator

    def process(self, df):

        selection = processor.PackedSelection()

        # Lepton candidates
        jets, muons, electrons = setup_candidates(df)
        tight_muons = muons[muons.mediumId & (muons.pt>10)]
        loose_muons = muons[muons.mediumId & (muons.pt>20)]
        loose_electrons = electrons[electrons.looseId & (electrons.pt>10)]
        tight_electrons = electrons[electrons.tightId & (electrons.pt>20)]

        # Jets
        jet_acceptance = (jets.eta<2.4)&(jets.eta>-2.4)
        jet_fractions = (jets.chf>0.1)&(jets.nhf<0.8)
        goodjets = jets[jet_fractions & jet_acceptance]

        # MET
        df["dPFCalo"] = 1 - df["CaloMET_pt"] / df["MET_pt"]

        # Selection
        selection.add("met_250", df["MET_pt"]>250)
        selection.add("leadjet_100", goodjets.pt.max()>100)
        selection.add("no_loose_leptons", (loose_electrons.counts==0) & (loose_muons.counts==0))
        selection.add("pf_calo_0p4",np.abs(df["dPFCalo"]) < 0.4)
        # print(jets[jets.pt>100])

        output = self.accumulator.identity()

        dataset = "inclusive"

        selection_masks = {
            "sr" : selection.all(*selection.names),
            "inclusive" : selection.all()
        }

        for seltag, mask in selection_masks.items():
            dataset = seltag
            # All jets
            output['jeteta'].fill(dataset=dataset,
                                    jeteta=jets[mask].eta.flatten())
            output['jetpt'].fill(dataset=dataset,
                                    jetpt=jets[mask].pt.flatten())
            
            # Leading jet (has to be in acceptance)
            output['jet0eta'].fill(dataset=dataset,
                                    jeteta=goodjets[mask].eta[goodjets[mask].pt.argmax()].flatten())
            output['jet0pt'].fill(dataset=dataset,
                                    jetpt=goodjets[mask].pt.max().flatten())
            

            # MET
            output['dpfcalo'].fill(dataset=dataset,
                                    dpfcalo=df["dPFCalo"][mask])
            output['met'].fill(dataset=dataset,
                                    met=df["MET_pt"][mask])
        return output

    def postprocess(self, accumulator):
        return accumulator


def main():
    fileset = {
        "Znunu_ht200to400" : [
            "./data/FFD69E5A-A941-2D41-A629-9D62D9E8BE9A.root"
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
                                  executor_args={'workers': 4, 'function_args': {'flatten': False}},
                                  chunksize=500000,
                                 )

    outdir = "out"
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    for name in ['jet0pt', 'jet0eta','jetpt','jeteta','met','dpfcalo']:
        fig, ax, _ = hist.plot1d(output[name],overlay="dataset",overflow='all')
        fig.suptitle(name)
        # ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_ylim(0.1, 1e8)
        fig.savefig(pjoin(outdir, "{}.pdf".format(name)))

    with lz4f.open("hists.cpkl.lz4", mode="wb", compression_level=5) as fout:
        cloudpickle.dump(output, fout)

if __name__ == "__main__":
    main()