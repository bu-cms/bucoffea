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

from bucoffea.monojet.definitions import monojet_accumulator, setup_candidates, monojet_regions
from bucoffea.helpers import min_dphi_jet_met, recoil, mt


class monojetProcessor(processor.ProcessorABC):
    def __init__(self, year="2018"):
        self.year=year
        self._accumulator = monojet_accumulator()

    @property
    def accumulator(self):
        return self._accumulator

    def process(self, df):

        # Candidates
        # Already pre-filtered!
        # All leptons are at least loose
        # Check out setup_candidates for filtering details
        ak4, ak8, muons, electrons, taus, photons = setup_candidates(df, cfg)

        # Muons
        is_tight_mu = muons.tightId \
                      & (muons.iso < cfg.MUON.CUTS.TIGHT.ISO) \
                      & (muons.pt>cfg.MUON.CUTS.TIGHT.PT) \
                      & (np.abs(muons.eta)<cfg.MUON.CUTS.TIGHT.ETA)

        dimuons = muons.distincts()
        dimuon_charge = dimuons.i0['charge'] + dimuons.i1['charge']

        df['MT_mu'] = ((muons.counts==1) * mt(muons.pt, muons.phi, df['MET_pt'], df['MET_phi'])).max()

        # Electrons
        is_tight_electron = electrons.tightId \
                            & (electrons.pt > cfg.ELECTRON.CUTS.TIGHT.PT) \
                            & (np.abs(electrons.eta) < cfg.ELECTRON.CUTS.TIGHT.ETA)

        dielectrons = electrons.distincts()
        dielectron_charge = dielectrons.i0['charge'] + dielectrons.i1['charge']

        df['MT_el'] = ((electrons.counts==1) * mt(electrons.pt, electrons.phi, df['MET_pt'], df['MET_phi'])).max()


        # ak4
        jet_acceptance = np.abs(ak4.eta)<2.4

        # B ak4
        btag_cut = cfg.BTAG.CUTS[cfg.BTAG.algo][cfg.BTAG.wp]
        jet_btagged = getattr(ak4, cfg.BTAG.algo) > btag_cut
        bjets = ak4[(ak4.clean==1) & jet_acceptance & jet_btagged]

        # MET
        df["dPFCalo"] = 1 - df["CaloMET_pt"] / df["MET_pt"]
        df["minDPhiJetMet"] = min_dphi_jet_met(ak4[ak4.clean==1], df['MET_phi'], njet=4, ptmin=30)
        df['recoil_pt'], df['recoil_phi'] = recoil(df['MET_pt'],df['MET_phi'], electrons, muons)

        selection = processor.PackedSelection()

        # Common selection
        selection.add('inclusive', np.ones(df.size)==1)
        selection.add('filt_met', df['Flag_METFilters'])
        selection.add('trig_met', df[cfg.TRIGGERS.MET])
        selection.add('veto_ele', electrons.counts==0)
        selection.add('veto_muo', muons.counts==0)
        selection.add('veto_photon', photons.counts==0)
        selection.add('veto_tau',taus.counts==0)
        selection.add('veto_b',bjets.counts==0)
        selection.add('dphijm',df['minDPhiJetMet'] > cfg.SELECTION.SIGNAL.MINDPHIJM)
        selection.add('dpfcalo',np.abs(df['dPFCalo']) < cfg.SELECTION.SIGNAL.DPFCALO)
        selection.add('recoil', df['recoil_pt']>cfg.SELECTION.SIGNAL.RECOIL)

        # AK4 Jet
        leadak4_index=ak4.pt.argmax()
        leadak4_pt_eta = (ak4.pt.max() > cfg.SELECTION.SIGNAL.leadak4.PT) \
                         & (np.abs(ak4.eta[leadak4_index]) < cfg.SELECTION.SIGNAL.leadak4.ETA).any()
        selection.add('leadak4_pt_eta', leadak4_pt_eta)
        
        selection.add('leadak4_id',(ak4.tightId[leadak4_index] \
                                                    & (ak4.chf[leadak4_index] >cfg.SELECTION.SIGNAL.leadak4.CHF) \
                                                    & (ak4.nhf[leadak4_index]<cfg.SELECTION.SIGNAL.leadak4.NHF)).any())

        # AK8 Jet
        leadak8_index=ak8.pt.argmax()
        leadak8_pt_eta = (ak8.pt.max() > cfg.SELECTION.SIGNAL.leadak8.PT) \
                         & (np.abs(ak8.eta[leadak8_index]) < cfg.SELECTION.SIGNAL.leadak8.ETA).any()
        selection.add('leadak8_pt_eta', leadak8_pt_eta)

        selection.add('leadak8_id',(ak8.tightId[leadak8_index]).any())

        # Mono-V selection
        selection.add('leadak8_tau21', ((ak8.tau2[leadak8_index] / ak8.tau1[leadak8_index]) < cfg.SELECTION.SIGNAL.LEADAK8.TAU21).any())
        selection.add('leadak8_mass', ((ak8.mass[leadak8_index] > cfg.SELECTION.SIGNAL.LEADAK8.MASS.MIN) \
                                    & (ak8.mass[leadak8_index] < cfg.SELECTION.SIGNAL.LEADAK8.MASS.MAX)).any())
        
        selection.add('veto_vtag', ~selection.all("leadak8_pt_eta", "leadak8_id", "leadak8_tau21", "leadak8_mass"))

        # Dimuon CR
        selection.add('at_least_one_tight_mu', is_tight_mu.any())
        selection.add('dimuon_mass', ((dimuons.mass > cfg.SELECTION.CONTROL.DOUBLEMU.MASS.MIN) \
                                    & (dimuons.mass < cfg.SELECTION.CONTROL.DOUBLEMU.MASS.MAX)).any())
        selection.add('dimuon_charge', (dimuon_charge==0).any())
        selection.add('two_muons', muons.counts==2)

        # Single muon CR
        selection.add('one_muon', muons.counts==1)
        selection.add('mt_mu', df['MT_mu'] < cfg.SELECTION.CONTROL.SINGLEMU.MT)

        # Diele CR
        selection.add('one_electron', electrons.counts==1)
        selection.add('two_electrons', electrons.counts==2)
        selection.add('at_least_one_tight_el', is_tight_electron.any())

        selection.add('dielectron_mass', ((dielectrons.mass > cfg.SELECTION.CONTROL.DOUBLEEL.MASS.MIN)  \
                                        & (dielectrons.mass < cfg.SELECTION.CONTROL.DOUBLEEL.MASS.MAX)).any())
        selection.add('dielectron_charge', (dielectron_charge==0).any())
        selection.add('two_electrons', electrons.counts==2)

        # Single Ele CR
        selection.add('mt_el', df['MT_el'] < cfg.SELECTION.CONTROL.SINGLEEL.MT)


        # Fill histograms
        weight = df['Generator_weight']
        regions = monojet_regions()
        output = self.accumulator.identity()
        for region, cuts in regions.items():

            # Cutflow plot for signal and control regions
            if any(x in region for x in ["sr", "cr"]):
                output['cutflow_' + region]['all']+=df.size
                for icut, cutname in enumerate(cuts):
                    output['cutflow_' + region][cutname] += selection.all(*cuts[:icut+1]).sum()

            dataset = df['dataset']

            mask = selection.all(*cuts)

            # Multiplicities
            def fill_mult(name, candidates):
                output[name].fill(
                                  dataset=dataset,
                                  region=region,
                                  multiplicity=candidates[mask].counts,
                                  weight=weight
                                  )

            fill_mult('ak8_mult', ak8)
            fill_mult('ak4_mult', ak4)
            fill_mult('bjet_mult',bjets)
            fill_mult('loose_ele_mult',electrons)
            fill_mult('tight_ele_mult',electrons[is_tight_electron])
            fill_mult('loose_muo_mult',muons)
            fill_mult('tight_muo_mult',muons[is_tight_mu])
            fill_mult('tau_mult',taus)
            fill_mult('photon_mult',photons)

            # All ak4
            output['ak4eta'].fill(
                                  dataset=dataset,
                                  region=region,
                                  jeteta=ak4[mask].eta.flatten(),
                                  weight=weight
                                  )
            output['ak4pt'].fill(
                                 dataset=dataset,
                                 region=region,
                                 jetpt=ak4[mask].pt.flatten(),
                                 weight=weight
                                 )

            # Leading ak4
            leadak4_indices = ak4.pt.argmax()
            output['ak4eta0'].fill(
                                   dataset=dataset,
                                   region=region,
                                   jeteta=ak4[leadak4_indices].eta[mask].flatten(),
                                   weight=weight
                                   )
            output['ak4pt0'].fill(
                                  dataset=dataset,
                                  region=region,
                                  jetpt=ak4[leadak4_indices].pt[mask].flatten(),
                                  weight=weight
                                  )

            # All ak8
            output['ak8eta'].fill(
                                  dataset=dataset,
                                  region=region,
                                  jeteta=ak8[mask].eta.flatten(),
                                  weight=weight
                                  )
            output['ak8pt'].fill(
                                 dataset=dataset,
                                 region=region,
                                 jetpt=ak8[mask].pt.flatten(),
                                 weight=weight
                                 )
            output['ak8mass'].fill(
                                 dataset=dataset,
                                 region=region,
                                 mass=ak8[mask].mass.flatten(),
                                 weight=weight
                                 )

            # Leading ak8
            leadak8_indices = ak8.pt.argmax()
            output['ak8eta0'].fill(
                                   dataset=dataset,
                                   region=region,
                                   jeteta=ak8[leadak8_indices].eta[mask].flatten(),
                                   weight=weight
                                   )
            output['ak8pt0'].fill(
                                  dataset=dataset,
                                  region=region,
                                  jetpt=ak8[leadak8_indices].pt[mask].flatten(),
                                  weight=weight
                                  )
            output['ak8mass0'].fill(
                                 dataset=dataset,
                                 region=region,
                                 mass=ak8[leadak8_indices].mass[mask].flatten(),
                                 weight=weight
                                 )

            # B tag discriminator
            output['ak4btag'].fill(
                                   dataset=dataset,
                                   region=region,
                                   btag=getattr(ak4[mask&jet_acceptance], cfg.BTAG.algo).flatten(),
                                   weight=weight
                                   )

            # MET
            output['dpfcalo'].fill(
                                   dataset=dataset,
                                   region=region,
                                   dpfcalo=df["dPFCalo"][mask],
                                   weight=weight
                                   )
            output['met'].fill(
                               dataset=dataset,
                               region=region,
                               met=df["MET_pt"][mask],
                               weight=weight
                                )
            output['recoil'].fill(
                               dataset=dataset,
                               region=region,
                               recoil=df["recoil_pt"][mask],
                               weight=weight
                                )
            output['dphijm'].fill(
                                  dataset=dataset,
                                  region=region,
                                  dphi=df["minDPhiJetMet"][mask],
                                  weight=weight
                                  )

            # Muons
            output['muon_pt'].fill(
                dataset=dataset,
                region=region,
                pt=muons.pt[mask].flatten(),
                weight=weight
            )
            output['muon_mt'].fill(
                dataset=dataset,
                region=region,
                mt=df['MT_mu'][mask],
                weight=weight
            )
            output['muon_eta'].fill(
                dataset=dataset,
                region=region,
                eta=muons.eta[mask].flatten(),
                weight=weight
            )
            # Dimuon
            output['dimuon_pt'].fill(
                dataset=dataset,
                region=region,
                pt=dimuons.pt[mask].flatten(),
                weight=weight
            )
            output['dimuon_eta'].fill(
                dataset=dataset,
                region=region,
                eta=dimuons.eta[mask].flatten(),
                weight=weight
            )
            output['dimuon_mass'].fill(
                dataset=dataset,
                region=region,
                dilepton_mass=dimuons.mass[mask].flatten(),
                weight=weight
            )

            # Electrons
            output['electron_pt'].fill(
                dataset=dataset,
                region=region,
                pt=electrons.pt[mask].flatten(),
                weight=weight
            )
            output['electron_mt'].fill(
                dataset=dataset,
                region=region,
                mt=df['MT_el'][mask],
                weight=weight
            )
            output['electron_eta'].fill(
                dataset=dataset,
                region=region,
                eta=electrons.eta[mask].flatten(),
                weight=weight
            )
            # Dielectron
            output['dielectron_pt'].fill(
                dataset=dataset,
                region=region,
                pt=dielectrons.pt[mask].flatten(),
                weight=weight
            )
            output['dielectron_eta'].fill(
                dataset=dataset,
                region=region,
                eta=dielectrons.eta[mask].flatten(),
                weight=weight
            )
            output['dielectron_mass'].fill(
                dataset=dataset,
                region=region,
                dilepton_mass=dielectrons.mass[mask].flatten(),
                weight=weight
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
        "ttbardm_mmed10000_mchi1_nanoaodv5" : [
            "./data/ttbardm_mmed10000_mchi1_nanoaodv5/64888F08-B888-ED40-84A5-F321A4BEAC27.root",
            "./data/ttbardm_mmed10000_mchi1_nanoaodv5/DBCA1773-BF28-E54F-8046-5EB316C2F725.root"
        ],
        # 'dy_zpt200_m50_mlm_2016_nanov4' : [
            # "./data/dy_zpt200_m50_mlm_2016_nanov4.root"
        # ]
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
    from matplotlib import pyplot as plt
    outdir = "out"
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    for name in output.keys():
        if name.startswith("_"):
            continue
        if name.startswith("cutflow"):
            continue
        if np.sum(output[name].values().values()) == 0:
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
        plt.close(fig)



def debug_print_cutflows(output):
    """Pretty-print cutflow data to the terminal."""
    import tabulate
    for cutflow_name in [ x for x in output.keys() if x.startswith("cutflow")]:
        if not len(output[cutflow_name]):
            continue
        table = []
        print("----")
        print(cutflow_name)
        print("----")
        for cut, count in sorted(output[cutflow_name].items(), key=lambda x:x[1], reverse=True):
            table.append([cut, count])
        print(tabulate.tabulate(table, headers=["Cut", "Passing events"]))



if __name__ == "__main__":
    main()