from coffea import hist
from coffea.analysis_objects import JaggedCandidateArray
import coffea.processor as processor
from awkward import JaggedArray
import numpy as np

import lz4.frame as lz4f
import cloudpickle
from copy import deepcopy
import os
pjoin = os.path.join
from collections import defaultdict
os.environ["ENV_FOR_DYNACONF"] = "era2016"
os.environ["SETTINGS_FILE_FOR_DYNACONF"] = os.path.abspath("config.yaml")
from dynaconf import settings as cfg

def setup_candidates(df):
    muons = JaggedCandidateArray.candidatesfromcounts(
        df['nMuon'],
        pt=df['Muon_pt'],
        eta=df['Muon_eta'],
        phi=df['Muon_phi'],
        mass=df['Muon_mass'],
        charge=df['Muon_charge'],
        mediumId=df['Muon_mediumId'],
        iso=df["Muon_pfRelIso04_all"],
        tightId=df['Muon_tightId']
    )

    electrons = JaggedCandidateArray.candidatesfromcounts(
        df['nElectron'],
        pt=df['Electron_pt'],
        eta=df['Electron_eta'],
        phi=df['Electron_phi'],
        mass=df['Electron_mass'],
        charge=df['Electron_charge'],
        looseId=(df['Electron_cutBased']>=1),
        tightId=(df['Electron_cutBased']==4)
    )
    taus = JaggedCandidateArray.candidatesfromcounts(
        df['nTau'],
        pt=df['Tau_pt'],
        eta=df['Tau_eta'],
        phi=df['Tau_phi'],
        mass=df['Tau_mass'],
        decaymode=df['Tau_idDecayModeNewDMs'],
        clean=df['Tau_cleanmask'],
        iso=df['Tau_idMVAnewDM2017v2'],
    )

    taus = taus[ (taus.clean==1) \
                         & (taus.decaymode) \
                         & (taus.pt > cfg.TAU.CUTS.PT)\
                         & (np.abs(taus.eta) < cfg.TAU.CUTS.ETA) \
                         & ((taus.iso&2)==2)]

    photons = JaggedCandidateArray.candidatesfromcounts(
        df['nPhoton'],
        pt=df['Photon_pt'],
        eta=df['Photon_eta'],
        phi=df['Photon_phi'],
        mass=df['Photon_mass'],
        id=(df['Photon_cutBased']==1) & (df['Photon_electronVeto']==1),
        clean=df['Photon_cleanmask'],
    )
    photons = photons[(photons.clean==1) \
              & photons.id \
              & (photons.pt > cfg.PHOTON.CUTS.pt) \
              & (np.abs(photons.eta) < cfg.PHOTON.CUTS.eta)]
    jets = JaggedCandidateArray.candidatesfromcounts(
        df['nJet'],
        pt=df['Jet_pt'],
        eta=df['Jet_eta'],
        phi=df['Jet_phi'],
        mass=df['Jet_mass'],

        # Jet ID bit mask:
        # Bit 0 = Loose
        # Bit 1 = Tight
        tightId=(df['Jet_jetId']&2) == 2,
        csvv2=df["Jet_btagCSVV2"],
        deepcsv=df['Jet_btagDeepB'],
        # nef=df['Jet_neEmEF'],
        nhf=df['Jet_neHEF'],
        chf=df['Jet_chHEF'],
        clean=df['Jet_cleanmask']
        # cef=df['Jet_chEmEF'],
    )
    jets = jets[jets.clean==1]
    return jets, muons, electrons, taus, photons

def define_dphi_jet_met(jets, met_phi, njet=4, ptmin=30):
    """Calculate minimal delta phi between jets and met

    :param jets: Jet candidates to use, must be sorted by pT
    :type jets: JaggedCandidateArray
    :param met_phi: MET phi values, one per event
    :type met_phi: array
    :param njet: Number of leading jets to consider, defaults to 4
    :type njet: int, optional
    """

    # Use the first njet jets with pT > ptmin
    jets=jets[jets.pt>ptmin]
    jets = jets[:,:njet]

    dphi = np.abs((jets.phi - met_phi + np.pi) % (2*np.pi)  - np.pi)

    return dphi.min()

def monojet_accumulator():
    dataset_ax = hist.Cat("dataset", "Primary dataset")
    region_ax = hist.Cat("region", "Selection region")
    met_ax = hist.Bin("met", r"$p_{T}^{miss}$ (GeV)", 100, 0, 1000)
    jet_pt_ax = hist.Bin("jetpt", r"$p_{T}$ (GeV)", 100, 0, 1000)
    jet_eta_ax = hist.Bin("jeteta", r"$\eta$ (GeV)", 50, -5, 5)
    dpfcalo_ax = hist.Bin("dpfcalo", r"$1-Calo/PF$", 20, -1, 1)
    btag_ax = hist.Bin("btag", r"B tag discriminator", 20, 0, 1)
    multiplicity_ax = hist.Bin("multiplicity", r"multiplicity", 10, -0.5, 9.5)
    dphi_ax = hist.Bin("dphi", r"$\Delta\phi$", 50, 0, 3*np.pi)

    Hist = hist.Hist
    items = {}
    items["met"] = Hist("Counts", dataset_ax, region_ax, met_ax)
    items["jet0pt"] = Hist("Counts", dataset_ax, region_ax, jet_pt_ax)
    items["jet0eta"] = Hist("Counts", dataset_ax, region_ax, jet_eta_ax)
    items["jetpt"] = Hist("Counts", dataset_ax, region_ax, jet_pt_ax)
    items["jeteta"] = Hist("Counts", dataset_ax, region_ax, jet_eta_ax)
    items["dpfcalo"] = Hist("Counts", dataset_ax, region_ax, dpfcalo_ax)
    items["jetbtag"] = Hist("Counts", dataset_ax, region_ax, btag_ax)
    items["jet_mult"] = Hist("Jets", dataset_ax, region_ax, multiplicity_ax)
    items["bjet_mult"] = Hist("B Jets", dataset_ax, region_ax, multiplicity_ax)
    items["loose_ele_mult"] = Hist("Loose electrons", dataset_ax, region_ax, multiplicity_ax)
    items["tight_ele_mult"] = Hist("Tight electrons", dataset_ax, region_ax, multiplicity_ax)
    items["loose_muo_mult"] = Hist("Loose muons", dataset_ax, region_ax, multiplicity_ax)
    items["tight_muo_mult"] = Hist("Tight muons", dataset_ax, region_ax, multiplicity_ax)
    items["tau_mult"] = Hist("Taus", dataset_ax, region_ax, multiplicity_ax)
    items["photon_mult"] = Hist("Photons", dataset_ax, region_ax, multiplicity_ax)
    items["dphijm"] = Hist("min(4 leading jets, MET)", dataset_ax, region_ax, dphi_ax)

    items["cutflow_sr_j"] = processor.defaultdict_accumulator(int)
    items["cutflow_sr_v"] = processor.defaultdict_accumulator(int)

    return  processor.dict_accumulator(items)

class monojetProcessor(processor.ProcessorABC):
    def __init__(self, year="2018"):
        self.year=year
        self._accumulator = monojet_accumulator()

    @property
    def accumulator(self):
        return self._accumulator

    def process(self, df):

        # Lepton candidates
        jets, muons, electrons, taus, photons = setup_candidates(df)
        loose_muons = muons[muons.mediumId & (muons.pt>10) & (muons.iso < 0.25)]
        tight_muons = muons[muons.mediumId & (muons.pt>20) & (muons.iso < 0.15)]
        loose_electrons = electrons[electrons.looseId & (electrons.pt>10)]
        tight_electrons = electrons[electrons.tightId & (electrons.pt>20)]

        # Jets
        jet_acceptance = (jets.eta<2.4)&(jets.eta>-2.4)
        jet_fractions = (jets.chf>0.1)&(jets.nhf<0.8)


        btag_cut = cfg.BTAG.CUTS[cfg.BTAG.algo][cfg.BTAG.wp]

        jet_btagged = getattr(jets, cfg.BTAG.algo) > btag_cut
        bjets = jets[(jets.clean==1) & jet_acceptance & jet_btagged]

        # MET
        df["dPFCalo"] = 1 - df["CaloMET_pt"] / df["MET_pt"]
        df["minDPhiJetMet"] = define_dphi_jet_met(jets[jets.clean==1], df['MET_phi'], njet=4, ptmin=30)

        selection = processor.PackedSelection()

        selection.add('inclusive', np.ones(df.size)==1)
        selection.add('filt_met', df['Flag_METFilters'])
        selection.add('trig_met', df[cfg.TRIGGERS.MET])
        selection.add('veto_ele', loose_electrons.counts==0)
        selection.add('veto_muo', loose_muons.counts==0)
        selection.add('veto_photon', photons.counts==0)
        selection.add('veto_tau',taus.counts==0)
        selection.add('veto_b',bjets.counts==0)

        leadjet_index=jets.pt.argmax()
        selection.add('leadjet_pt_eta', (jets.pt.max() > cfg.SELECTION.SIGNAL.LEADJET.PT) \
                                                         & (np.abs(jets.eta[leadjet_index]) < cfg.SELECTION.SIGNAL.LEADJET.ETA).any())
        selection.add('leadjet_id',(jets.tightId[leadjet_index] \
                                                    & (jets.chf[leadjet_index] >cfg.SELECTION.SIGNAL.LEADJET.CHF) \
                                                    & (jets.nhf[leadjet_index]<cfg.SELECTION.SIGNAL.LEADJET.NHF)).any())
        selection.add('dphijm',df['minDPhiJetMet'] > cfg.SELECTION.SIGNAL.MINDPHIJM)
        selection.add('dpfcalo',np.abs(df['dPFCalo']) < cfg.SELECTION.SIGNAL.DPFCALO)
        selection.add('met_signal', df['MET_pt']>cfg.SELECTION.SIGNAL.MET)
        selection.add('tau21', np.ones(df.size)==0)
        selection.add('veto_vtag', np.ones(df.size)==1)


        common_cuts = [
            'filt_met',
            'trig_met',
            'veto_ele',
            'veto_muo',
            'veto_photon',
            'veto_tau',
            'veto_b',
            'leadjet_pt_eta',
            'leadjet_id',
            'dphijm',
            'dpfcalo',
            'met_signal'
        ]
        regions = {}
        regions['inclusive'] = ['inclusive']
        regions['sr_v'] = common_cuts + ['tau21']
        regions['sr_j'] = common_cuts + ['veto_vtag']

        output = self.accumulator.identity()

        for region, cuts in regions.items():

            # Cutflow plot for signal and control regions
            if any(x in region for x in ["sr", "cr"]):
                output['cutflow_' + region]['all']+=df.size
                for icut, cutname in enumerate(cuts):
                    output['cutflow_' + region][cutname] += selection.all(*cuts[:icut+1]).sum()

            dataset = df['dataset']

            mask = selection.all(*cuts)
            print(region, mask)
            # Multiplicities
            def fill_mult(name, candidates):
                output[name].fill(
                                  dataset=dataset,
                                  region=region,
                                  multiplicity=candidates[mask].counts
                                  )

            fill_mult('jet_mult', jets)
            fill_mult('bjet_mult',bjets)
            fill_mult('loose_ele_mult',loose_electrons)
            fill_mult('tight_ele_mult',tight_electrons)
            fill_mult('loose_muo_mult',loose_muons)
            fill_mult('tight_muo_mult',tight_muons)
            fill_mult('tau_mult',taus)
            fill_mult('photon_mult',photons)

            # All jets
            output['jeteta'].fill(
                                  dataset=dataset,
                                  region=region,
                                  jeteta=jets[mask].eta.flatten()
                                  )
            output['jetpt'].fill(
                                 dataset=dataset,
                                 region=region,
                                 jetpt=jets[mask].pt.flatten()
                                 )

            # Leading jet
            leadjet_indices = jets.pt.argmax()
            output['jet0eta'].fill(
                                   dataset=dataset,
                                   region=region,
                                   jeteta=jets[leadjet_indices].eta[mask].flatten()
                                   )
            output['jet0pt'].fill(
                                  dataset=dataset,
                                  region=region,
                                  jetpt=jets[leadjet_indices].pt[mask].flatten()
                                  )

            # B tag discriminator
            output['jetbtag'].fill(
                                   dataset=dataset,
                                   region=region,
                                   btag=getattr(jets[mask&jet_acceptance], cfg.BTAG.algo).flatten()
                                   )

            # MET
            output['dpfcalo'].fill(
                                   dataset=dataset,
                                   region=region,
                                   dpfcalo=df["dPFCalo"][mask]
                                   )
            output['met'].fill(
                               dataset=dataset,
                               region=region,
                               met=df["MET_pt"][mask]
                                )
            output['dphijm'].fill(
                                  dataset=dataset,
                                  region=region,
                                  dphi=df["minDPhiJetMet"][mask]
                                  )
        return output

    def postprocess(self, accumulator):
        return accumulator


def main():
    fileset = {
        # "Znunu_ht200to400" : [
        #     "./data/FFD69E5A-A941-2D41-A629-9D62D9E8BE9A.root"
        # ],
        # "NonthDM" : [
        #     "./data/24EE25F5-FB54-E911-AB96-40F2E9C6B000.root"
        # ]
        "TTbarDM" : [
            "./data/A13AF968-8A88-754A-BE73-7264241D71D5.root"
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
                                  executor_args={'workers': 4, 'function_args': {'flatten': True}},
                                  chunksize=500000,
                                 )
    with lz4f.open("hists.cpkl.lz4", mode="wb", compression_level=5) as fout:
        cloudpickle.dump(output, fout)

    # Debugging / testing output
    debug_plot_output(output)
    debug_print_cutflows(output)

def debug_plot_output(output):
    """Dump all histograms as PDF."""
    outdir = "out"
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    for name in output.keys():
        if name.startswith("_"):
            continue
        if name.startswith("cutflow"):
            continue
        fig, ax, _ = hist.plot1d(
            output[name]. project('dataset'),
            overlay='region',
            overflow='all',
            )
        fig.suptitle(name)
        # ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_ylim(0.1, 1e8)
        fig.savefig(pjoin(outdir, "{}.pdf".format(name)))



def debug_print_cutflows(output):
    """Pretty-print cutflow data to the terminal."""
    import tabulate
    for cutflow_name in [ x for x in output.keys() if x.startswith("cutflow")]:
        table = []
        print("----")
        print(cutflow_name)
        print("----")
        for cut, count in sorted(output[cutflow_name].items(), key=lambda x:x[1], reverse=True):
            table.append([cut, count])
        print(tabulate.tabulate(table, headers=["Cut", "Passing events"]))



if __name__ == "__main__":
    main()