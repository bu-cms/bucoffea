import os
import numpy as np

from coffea import hist
import coffea.processor as processor
from coffea.analysis_objects import JaggedCandidateArray
from coffea.util import load, save
from collections import defaultdict


from dynaconf import settings as cfg


from bucoffea.monojet.definitions import monojet_accumulator, monojet_evaluator, setup_candidates, setup_gen_candidates,monojet_regions
from bucoffea.helpers import min_dphi_jet_met, recoil, mt, weight_shape, bucoffea_path
from bucoffea.helpers.gen import find_gen_dilepton

def is_w_dataset(dataset):
    """Dummy implementation"""
    return "wjet" in dataset
def is_z_dataset(dataset)
    """Dummy implementation"""
    return "zjet" in dataset
class monojetProcessor(processor.ProcessorABC):
    def __init__(self, year="2017",blind=True):
        self._year=year
        self._blind=blind
        self._accumulator = monojet_accumulator()
        cfg.SETTINGS_FILE_FOR_DYNACONF = bucoffea_path("config/monojet.yaml")
        cfg.ENV_FOR_DYNACONF = f"era{self._year}"
        cfg.reload()
        self._evaluator = monojet_evaluator(cfg)
    @property
    def accumulator(self):
        return self._accumulator

    def process(self, df):
        # Gen info
        gen = setup_gen_candidates(df)
        if(is_z_dataset(df['dataset'])):
            gen_v = find_gen_dilepton(gen, pdgsum=0)
        elif(is_w_dataset(df['dataset'])):
            gen_v = find_gen_dilepton(gen, pdgsum=1)
        else:
            gen_v = np.zeros(df.size)
        # Candidates
        # Already pre-filtered!
        # All leptons are at least loose
        # Check out setup_candidates for filtering details
        ak4, ak8, muons, electrons, taus, photons = setup_candidates(df, cfg)

        # Muons
        is_tight_muon = muons.tightId \
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

        # B tagged ak4
        btag_cut = cfg.BTAG.CUTS[cfg.BTAG.algo][cfg.BTAG.wp]
        jet_btagged = getattr(ak4, cfg.BTAG.algo) > btag_cut
        bjets = ak4[ jet_acceptance \
                     & jet_btagged \
                     & (ak4.pt>20) ]

        # Recoil
        df['recoil_pt'], df['recoil_phi'] = recoil(df['MET_pt'],df['MET_phi'], electrons, muons, photons)
        df["dPFCalo"] = (df['MET_pt'] - df["CaloMET_pt"]) / df["recoil_pt"]
        df["minDPhiJetRecoil"] = min_dphi_jet_met(ak4, df['recoil_phi'], njet=4, ptmin=30)
        df["minDPhiJetMet"] = min_dphi_jet_met(ak4, df['MET_phi'], njet=4, ptmin=30)
        selection = processor.PackedSelection()

        # Common selection
        selection.add('inclusive', np.ones(df.size)==1)
        selection.add('filt_met', df['Flag_METFilters'])
        selection.add('trig_met', df[cfg.TRIGGERS.MET])
        selection.add('veto_ele', electrons.counts==0)
        selection.add('veto_muo', muons.counts==0)
        selection.add('veto_photon', photons.counts==0)
        selection.add('veto_tau', taus.counts==0)
        selection.add('veto_b', bjets.counts==0)
        selection.add('mindphijr',df['minDPhiJetRecoil'] > cfg.SELECTION.SIGNAL.MINDPHIJR)
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
        selection.add('at_least_one_tight_mu', is_tight_muon.any())
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

        # Photon CR
        is_tight_photon = photons.mediumId \
                         & (photons.pt > cfg.PHOTON.CUTS.TIGHT.PT) \
                         & (np.abs(photons.eta) < cfg.PHOTON.CUTS.TIGHT.ETA)

        selection.add('one_photon', photons.counts==1)
        selection.add('at_least_one_tight_photon', is_tight_photon.any())

        # Fill histograms
        output = self.accumulator.identity()
        dataset = df['dataset']
        isdata = dataset.startswith("data_")

        all_weights = {}
        if isdata:
            weight = np.ones(df.size)
        else:
            weight = df['Generator_weight']

            all_weights["muon_id_tight"] = self._evaluator['muon_id_tight'](muons[is_tight_muon].pt, muons[is_tight_muon].eta).prod()
            all_weights["muon_id_loose"] = self._evaluator['muon_id_loose']( muons[~is_tight_muon].pt, muons[~is_tight_muon].eta).prod()

            all_weights["ele_id_tight"] = self._evaluator['ele_id_tight'](electrons[is_tight_electron].eta, electrons[is_tight_electron].pt).prod()
            all_weights["ele_id_loose"] = self._evaluator['ele_id_loose'](electrons[~is_tight_electron].eta, electrons[~is_tight_electron].pt).prod()

            all_weights["photon_id_tight"] = self._evaluator['photon_id_tight'](photons[is_tight_photon].eta, photons[is_tight_photon].pt).prod()
            all_weights["pileup"] = self._evaluator['pileup'](df['Pileup_nTrueInt'])

            for iw in all_weights.values():
                weight = weight * iw




        # Sum of all weights to use for normalization
        # TODO: Deal with systematic variations
        output['sumw'][dataset] += weight.sum()
        output['sumw2'][dataset] += (weight**2).sum()

        regions = monojet_regions()
        for region, cuts in regions.items():
            # Blinding
            if(self._blind and df['dataset'].startswith('data') and region.startswith('sr')):
                continue

            # Cutflow plot for signal and control regions
            if any(x in region for x in ["sr", "cr"]):
                output['cutflow_' + region]['all']+=df.size
                for icut, cutname in enumerate(cuts):
                    output['cutflow_' + region][cutname] += selection.all(*cuts[:icut+1]).sum()

            mask = selection.all(*cuts)


            # Save the event numbers of events passing this selection
            if cfg.RUN.SAVEEVENTS:
                output['selected_events'][region] += list(df['event'][mask])

            # Multiplicities
            def fill_mult(name, candidates):
                output[name].fill(
                                  dataset=dataset,
                                  region=region,
                                  multiplicity=candidates[mask].counts,
                                  weight=weight[mask]
                                  )

            fill_mult('ak8_mult', ak8)
            fill_mult('ak4_mult', ak4)
            fill_mult('bjet_mult',bjets)
            fill_mult('loose_ele_mult',electrons)
            fill_mult('tight_ele_mult',electrons[is_tight_electron])
            fill_mult('loose_muo_mult',muons)
            fill_mult('tight_muo_mult',muons[is_tight_muon])
            fill_mult('tau_mult',taus)
            fill_mult('photon_mult',photons)


            def ezfill(name, **kwargs):
                """Helper function to make filling easier."""
                output[name].fill(
                                  dataset=dataset,
                                  region=region,
                                  **kwargs
                                  )
            # Monitor weights
            for wname, wvalue in all_weights.items():
                ezfill("weights", weight_type=wname, weight_value=wvalue[mask])

            # All ak4
            # This is a workaround to create a weight array of the right dimension
            w_alljets = weight_shape(ak4[mask].eta, weight[mask])


            ezfill('ak4eta',    jeteta=ak4[mask].eta.flatten(), weight=w_alljets)
            ezfill('ak4pt',     jetpt=ak4[mask].pt.flatten(),   weight=w_alljets)

            # Leading ak4
            leadak4_indices = ak4.pt.argmax()
            w_leadak4 = weight_shape(ak4[leadak4_indices].eta[mask], weight[mask])
            ezfill('ak4eta0',   jeteta=ak4[leadak4_indices].eta[mask].flatten(),    weight=w_leadak4)
            ezfill('ak4pt0',    jetpt=ak4[leadak4_indices].pt[mask].flatten(),      weight=w_leadak4)

            # All ak8
            w_allak8 = weight_shape(ak8.eta[mask], weight[mask])

            ezfill('ak8eta',    jeteta=ak8[mask].eta.flatten(), weight=w_allak8)
            ezfill('ak8pt',     jetpt=ak8[mask].pt.flatten(),   weight=w_allak8)
            ezfill('ak8mass',   mass=ak8[mask].mass.flatten(),  weight=w_allak8)

            # Leading ak8
            leadak8_indices = ak8.pt.argmax()
            w_leadak8 = weight_shape(ak8[leadak8_indices].eta[mask], weight[mask])

            ezfill('ak8eta0',   jeteta=ak8[leadak8_indices].eta[mask].flatten(),    weight=w_leadak8)
            ezfill('ak8pt0',    jetpt=ak8[leadak8_indices].pt[mask].flatten(),      weight=w_leadak8 )
            ezfill('ak8mass0',  mass=ak8[leadak8_indices].mass[mask].flatten(),     weight=w_leadak8)

            # B tag discriminator
            btag = getattr(ak4, cfg.BTAG.ALGO)
            w_btag = weight_shape(btag[mask], weight[mask])
            ezfill('ak4btag', btag=btag[mask].flatten(), weight=w_btag )

            # MET
            ezfill('dpfcalo',   dpfcalo=df["dPFCalo"][mask],    weight=weight[mask] )
            ezfill('met',       met=df["MET_pt"][mask],         weight=weight[mask] )
            ezfill('recoil',    recoil=df["recoil_pt"][mask],   weight=weight[mask] )
            ezfill('dphijm',    dphi=df["minDPhiJetMet"][mask], weight=weight[mask] )

            # Muons
            w_allmu = weight_shape(muons.pt[mask], weight[mask])
            ezfill('muon_pt',   pt=muons.pt[mask].flatten(),    weight=w_allmu )
            ezfill('muon_mt',   mt=df['MT_mu'][mask],           weight=weight[mask])
            ezfill('muon_eta',  eta=muons.eta[mask].flatten(),  weight=w_allmu)
            # Dimuon
            w_dimu = weight_shape(dimuons.pt[mask], weight[mask])

            ezfill('dimuon_pt',     pt=dimuons.pt[mask].flatten(),              weight=w_dimu)
            ezfill('dimuon_eta',    eta=dimuons.eta[mask].flatten(),            weight=w_dimu)
            ezfill('dimuon_mass',   dilepton_mass=dimuons.mass[mask].flatten(), weight=w_dimu )

            # Electrons
            w_allel = weight_shape(electrons.pt[mask], weight[mask])
            ezfill('electron_pt',   pt=electrons.pt[mask].flatten(),    weight=w_allel)
            ezfill('electron_mt',   mt=df['MT_el'][mask],               weight=weight[mask])
            ezfill('electron_eta',  eta=electrons.eta[mask].flatten(),  weight=w_allel)

            # Dielectron
            w_diel = weight_shape(dielectrons.pt[mask], weight[mask])
            ezfill('dielectron_pt',     pt=dielectrons.pt[mask].flatten(),                  weight=w_diel)
            ezfill('dielectron_eta',    eta=dielectrons.eta[mask].flatten(),                weight=w_diel)
            ezfill('dielectron_mass',   dilepton_mass=dielectrons.mass[mask].flatten(),     weight=w_diel)
        return output

    def postprocess(self, accumulator):
        return accumulator

def extract_year(dataset):
    for x in [6,7,8]:
        if f"201{x}" in dataset:
            return 2010+x
    raise RuntimeError("Could not determine dataset year")


def main():
    fileset = {
        # "Znunu_ht200to400" : [
        #     "./data/FFD69E5A-A941-2D41-A629-9D62D9E8BE9A.root"
        # ],
        # "NonthDM" : [
        #     "./data/24EE25F5-FB54-E911-AB96-40F2E9C6B000.root"
        # ]
        # "ttbardm_mmed10000_mchi1_nanoaodv5" : [
        #     "./data/ttbardm_mmed10000_mchi1_nanoaodv5/64888F08-B888-ED40-84A5-F321A4BEAC27.root",
        #     "./data/ttbardm_mmed10000_mchi1_nanoaodv5/DBCA1773-BF28-E54F-8046-5EB316C2F725.root"
        # ],
        # "a2HDM" : ["./data/a2HDM_monoz_mH900_ma200_2017_v5.root"]
        # "monozvec18" : ["./data/monozll_vec_mmed_1500_mxd_1_2017_v5.root"],
        # "wz_p8_2018" : ["./data/wz_p8_2018_v5.root"],
        # "wz_p8_2017" : ["./data/tt_amc_2017_v5.root"],
        # "wz_p8_2018" : ["./data/tt_amc_2018_v5.root"],
        # "monozvec17" : ["./data/monozll_vec_mmed_1500_mxd_1_2018_v5"]
        # "data_met_run2016c_v4" : [
        #     "./data/data_met_run2016c_v4.root"
        # ]
        # 'dy_zpt200_m50_mlm_2016_nanov4' : [
            # "./data/dy_zpt200_m50_mlm_2016_nanov4.root"
        # ],
        # "tt_amc_2017_v5" : ["./data/tt_amc_2017_v5.root"],
        # "tt_amc_2018_v5" : ["./data/tt_amc_2018_v5.root"],
        "vector_monoj_mmed1500_mdm300_2017_v5" :[
        "./data/vector_monoj_mmed1500_mdm300_2017_v5.root"
        ]
    }

    years = list(set(map(extract_year, fileset.keys())))
    assert(len(years)==1)

    for dataset, filelist in fileset.items():
        newlist = []
        for file in filelist:
            if file.startswith("/store/"):
                newlist.append("root://cms-xrd-global.cern.ch//" + file)
            else: newlist.append(file)
        fileset[dataset] = newlist

    output = processor.run_uproot_job(fileset,
                                  treename='Events',
                                  processor_instance=monojetProcessor(years[0]),
                                  executor=processor.futures_executor,
                                  executor_args={'workers': 4, 'function_args': {'flatten': True}},
                                  chunksize=500000,
                                 )
    save(output, "monojet.coffea")
    # Debugging / testing output
    # debug_plot_output(output)
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
        if any([x in name for x in ['sumw','cutflow']]):
            continue
        if np.sum(output[name].values().values()) == 0:
            continue
        fig, ax, _ = hist.plot1d(
            output[name]. project('dataset').project("region","inclusive"),
            # overlay='region',
            overflow='all',
            )
        fig.suptitle(name)
        # ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_ylim(0.1, 1e8)
        fig.savefig(os.path.join(outdir, "{}.pdf".format(name)))
        plt.close(fig)



def debug_print_cutflows(output):
    """Pretty-print cutflow data to the terminal."""
    import tabulate
    for cutflow_name in [ x for x in output.keys() if x.startswith("cutflow") and x.endswith("j")]:
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