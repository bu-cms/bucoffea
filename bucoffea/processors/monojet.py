from coffea import hist
from coffea.analysis_objects import JaggedCandidateArray
import coffea.processor as processor
from awkward import JaggedArray
import numpy as np

import lz4.frame as lz4f
import cloudpickle
import kaptan
import os
pjoin = os.path.join

config = kaptan.Kaptan(handler="yaml")
with open("config.yaml","r") as cfgfile:
    config.import_config(cfgfile.read())

def setup_candidates(df):
    muons = JaggedCandidateArray.candidatesfromcounts(
        df['nMuon'].flatten(),
        pt=df['Muon_pt'].flatten(),
        eta=df['Muon_eta'].flatten(),
        phi=df['Muon_phi'].flatten(),
        mass=df['Muon_mass'].flatten(),
        charge=df['Muon_charge'].flatten(),
        mediumId=df['Muon_mediumId'].flatten(),
        iso=df["Muon_pfRelIso04_all"].flatten(),
        tightId=df['Muon_tightId'].flatten()
    )

    electrons = JaggedCandidateArray.candidatesfromcounts(
        df['nElectron'].flatten(),
        pt=df['Electron_pt'].flatten(),
        eta=df['Electron_eta'].flatten(),
        phi=df['Electron_phi'].flatten(),
        mass=df['Electron_mass'].flatten(),
        charge=df['Electron_charge'].flatten(),
        looseId=(df['Electron_cutBased']>=1).flatten(),
        tightId=(df['Electron_cutBased']==4).flatten()
    )
    taus = JaggedCandidateArray.candidatesfromcounts(
        df['nTau'].flatten(),
        pt=df['Tau_pt'].flatten(),
        eta=df['Tau_eta'].flatten(),
        phi=df['Tau_phi'].flatten(),
        mass=df['Tau_mass'].flatten(),
        decaymode=df['Tau_idDecayModeNewDMs'].flatten(),
        iso=df['Tau_idMVAnewDM2017v2'].flatten(),
    )
    jets = JaggedCandidateArray.candidatesfromcounts(
        df['nJet'].flatten(),
        pt=df['Jet_pt'].flatten(),
        eta=df['Jet_eta'].flatten(),
        phi=df['Jet_phi'].flatten(),
        mass=df['Jet_mass'].flatten(),
        tightId=df['Jet_jetId'].flatten(),
        csvv2=df["Jet_btagCSVV2"].flatten(),
        deepcsv=df['Jet_btagDeepB'].flatten(),
        # nef=df['Jet_neEmEF'].flatten(),
        nhf=df['Jet_neHEF'].flatten(),
        chf=df['Jet_chHEF'].flatten(),
        clean=df['Jet_cleanmask'].flatten()
        # cef=df['Jet_chEmEF'].flatten(),
    )
    return jets, muons, electrons, taus

def met_triggers():
    return ["HLT_PFMET170_NotCleaned",
    "HLT_PFMET170_NoiseCleaned",
    "HLT_PFMET170_HBHECleaned",
    "HLT_PFMET170_JetIdCleaned",
    "HLT_PFMET170_BeamHaloCleaned",
    "HLT_PFMET170_HBHE_BeamHaloCleaned"]

class monojetProcessor(processor.ProcessorABC):
    def __init__(self, year="2018"):
        self.year=year
        dataset_axis = hist.Cat("dataset", "Primary dataset")
        met_axis = hist.Bin("met", r"$p_{T}^{miss}$ (GeV)", 100, 0, 1000)
        jet_pt_axis = hist.Bin("jetpt", r"$p_{T}$ (GeV)", 100, 0, 1000)
        jet_eta_axis = hist.Bin("jeteta", r"$\eta$ (GeV)", 50, -5, 5)
        dpfcalo_axis = hist.Bin("dpfcalo", r"$1-Calo/PF$", 20, -1, 1)
        btag_axis = hist.Bin("btag", r"B tag discriminator", 20, 0, 1)
        multiplicity_axis = hist.Bin("multiplicity", r"multiplicity", 10, -0.5, 9.5)

        self._accumulator = processor.dict_accumulator({
            "met" : hist.Hist("Counts", dataset_axis, met_axis),
            "jet0pt" : hist.Hist("Counts", dataset_axis, jet_pt_axis),
            "jet0eta" : hist.Hist("Counts", dataset_axis, jet_eta_axis),
            "jetpt" : hist.Hist("Counts", dataset_axis, jet_pt_axis),
            "jeteta" : hist.Hist("Counts", dataset_axis, jet_eta_axis),
            "dpfcalo" : hist.Hist("Counts", dataset_axis, dpfcalo_axis),
            "jetbtag" : hist.Hist("Counts", dataset_axis, btag_axis),
            "jet_mult" : hist.Hist("Jets", dataset_axis, multiplicity_axis),
            "bjet_mult" : hist.Hist("B Jets", dataset_axis, multiplicity_axis),
            "loose_ele_mult" : hist.Hist("Loose electrons", dataset_axis, multiplicity_axis),
            "tight_ele_mult" : hist.Hist("Tight electrons", dataset_axis, multiplicity_axis),
            "loose_muo_mult" : hist.Hist("Loose muons", dataset_axis, multiplicity_axis),
            "tight_muo_mult" : hist.Hist("Tight muons", dataset_axis, multiplicity_axis),
            "veto_tau_mult" : hist.Hist("Veto taus", dataset_axis, multiplicity_axis),
        })

    @property
    def accumulator(self):
        return self._accumulator

    def process(self, df):

        selection = processor.PackedSelection()

        # Lepton candidates
        jets, muons, electrons, taus = setup_candidates(df)
        loose_muons = muons[muons.mediumId & (muons.pt>10) & (muons.iso < 0.25)]
        tight_muons = muons[muons.mediumId & (muons.pt>20) & (muons.iso < 0.15)]
        loose_electrons = electrons[electrons.looseId & (electrons.pt>10)]
        tight_electrons = electrons[electrons.tightId & (electrons.pt>20)]

        # Jets
        clean_jets = jets[jets.clean==1]
        jet_acceptance = (clean_jets.eta<2.4)&(clean_jets.eta>-2.4)
        jet_fractions = (clean_jets.chf>0.1)&(clean_jets.nhf<0.8)

        # B jets
        btag_algo = config.get(f"{self.year}.btagalgo")
        btag_wp = config.get(f"{self.year}.btagwp")
        btag_cut = config.get(f"{self.year}.{btag_algo}.wp.{btag_wp}")
        
        jet_btagged = getattr(clean_jets, btag_algo) > btag_cut
        bjets = clean_jets[jet_acceptance & jet_btagged]
        goodjets = clean_jets[jet_fractions & jet_acceptance & jet_btagged==0]

        # Taus
        veto_taus = taus[(taus.decaymode)&(taus.pt > 18)&((taus.iso&2)==2)]

        # MET
        df["dPFCalo"] = 1 - df["CaloMET_pt"] / df["MET_pt"]

        # Selection
        selection.add("filt_met", df["Flag_METFilters"])
        selection.add("trig_met", df["HLT_PFMET170_NotCleaned"])
        selection.add("met_250", df["MET_pt"]>250)
        selection.add("leadjet_100", (goodjets.counts>0) & (goodjets.pt.max()>100))
        selection.add("no_loose_leptons", (loose_electrons.counts==0) & (loose_muons.counts==0))
        selection.add("pf_calo_0p4",np.abs(df["dPFCalo"]) < 0.4)
        selection.add("bveto",bjets.counts==0)
        selection.add("tauveto",veto_taus.counts==0)
        # print(jets[jets.pt>100])

        output = self.accumulator.identity()

        dataset = "inclusive"

        selection_masks = {
            "sr" : selection.all(*selection.names),
            "inclusive" : selection.all()
        }

        for seltag, mask in selection_masks.items():
            dataset = seltag

            # Multiplicities
            output['jet_mult'].fill(dataset=dataset, multiplicity=clean_jets[mask].counts.flatten())
            output['bjet_mult'].fill(dataset=dataset, multiplicity=bjets[mask].counts.flatten())
            output['loose_ele_mult'].fill(dataset=dataset, multiplicity=loose_electrons[mask].counts.flatten())
            output['tight_ele_mult'].fill(dataset=dataset, multiplicity=tight_electrons[mask].counts.flatten())
            output['loose_muo_mult'].fill(dataset=dataset, multiplicity=loose_muons[mask].counts.flatten())
            output['tight_muo_mult'].fill(dataset=dataset, multiplicity=tight_muons[mask].counts.flatten())
            output['veto_tau_mult'].fill(dataset=dataset, multiplicity=veto_taus[mask].counts.flatten())




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
            
            # B tag discriminator
            output['jetbtag'].fill(dataset=dataset,
                                    btag=getattr(clean_jets[mask&jet_acceptance], btag_algo).flatten())

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
        # "Znunu_ht200to400" : [
        #     "./data/FFD69E5A-A941-2D41-A629-9D62D9E8BE9A.root"
        # ],
        "NonthDM" : [
            "./data/24EE25F5-FB54-E911-AB96-40F2E9C6B000.root"
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
                                  processor_instance=monojetProcessor(2017),
                                  executor=processor.futures_executor,
                                  executor_args={'workers': 4, 'function_args': {'flatten': False}},
                                  chunksize=500000,
                                 )

    outdir = "out"
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    for name in ['jet0pt', 'jet0eta','jetpt','jeteta','met','dpfcalo','jetbtag'] + list(filter(lambda x: "mult" in x, output.keys())):
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