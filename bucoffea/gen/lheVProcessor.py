import copy

import coffea.processor as processor
import numpy as np
from bucoffea.helpers import dphi, min_dphi_jet_met
from bucoffea.helpers.dataset import (is_lo_g, is_lo_g_ewk, is_lo_w, is_lo_z,
                                      is_nlo_g, is_nlo_g_ewk, is_nlo_w,
                                      is_nlo_z)
from bucoffea.helpers.gen import (fill_gen_v_info,
                                  setup_dressed_gen_candidates,
                                  setup_gen_candidates, setup_lhe_candidates,
                                  setup_lhe_cleaned_gen_fatjets,
                                  setup_lhe_cleaned_genjets)
from coffea import hist

Hist = hist.Hist
Bin = hist.Bin
Cat = hist.Cat

def selection(vphi, genjets, fatjets,lhe_charged_leps):
    selection = processor.PackedSelection()

    selection.add(
                  'at_least_one_jet',
                  genjets.counts>0
                  )
    selection.add(
                  'leadak4_pt_eta',
                  (genjets.pt.max() > 100) & (np.abs(genjets[genjets.pt.argmax()].eta.max()) < 2.4)
                  )
    selection.add(
                  'mindphijr',
                  min_dphi_jet_met(genjets, vphi, njet=4, ptmin=30) > 0.5
                  )

    selection.add(
                  'at_least_one_fat_jet',
                  fatjets.counts>0
                  )
    selection.add(
                  'leadak8_pt_eta',
                  (fatjets.pt.max() > 250) & (np.abs(fatjets[fatjets.pt.argmax()].eta.max()) < 2.4)
                  )
    msd = fatjets[fatjets.pt.argmax()].mass.max()
    selection.add(
                  'leadak8_mass',
                  (msd <65) & (msd<120)
                  )
    selection.add(
                  'no_forward_leps',
                  ~(np.abs(lhe_charged_leps.eta)>2.5).any()
                  )

    return selection



regions = {
    'monojet' : [
                 'at_least_one_jet',
                 'leadak4_pt_eta',
                 'mindphijr'
                 ],
    'monov' : [
                'at_least_one_fat_jet',
                'leadak8_pt_eta',
                'leadak8_mass',
                'mindphijr'
                ],
}

regions['monov_nomass'] = copy.deepcopy(regions['monov'])
regions['monov_nomass'].remove('leadak8_mass')

tmp = {}
for region in regions.keys():
    requirements = copy.deepcopy(regions[region])
    tmp[region+"_central"] = requirements
    tmp[region+"_central"].append("no_forward_leps")
regions.update(tmp)



class lheVProcessor(processor.ProcessorABC):
    def __init__(self):

        # Histogram setup
        dataset_ax = Cat("dataset", "Primary dataset")
        region_ax = Cat("region", "Region")

        vpt_ax = Bin("vpt",r"$p_{T}^{V}$ (GeV)", 100, 0, 2000)
        res_ax = Bin("res",r"pt: dressed / stat1 - 1", 80,-0.2,0.2)

        items = {}
        for tag in ['stat1','dress','lhe','combined']:
                items[f"gen_vpt_{tag}"] = Hist("Counts",
                                        dataset_ax,
                                        region_ax,
                                        vpt_ax)

        items["resolution"] = Hist("Counts",
                                dataset_ax,
                                res_ax)
        items['sumw'] = processor.defaultdict_accumulator(float)
        items['sumw2'] = processor.defaultdict_accumulator(float)

        self._accumulator = processor.dict_accumulator(items)

    @property
    def accumulator(self):
        return self._accumulator


    def process(self, df):
        output = self.accumulator.identity()
        dataset = df['dataset']

        genjets = setup_lhe_cleaned_genjets(df)
        genjets_ak8 = setup_lhe_cleaned_gen_fatjets(df)
        lhe = setup_lhe_candidates(df)

        lhe_charged_leps = lhe[
                                (np.abs(lhe.pdg)==11)
                                | (np.abs(lhe.pdg)==13)
                                | (np.abs(lhe.pdg)==15)
                                ]

        # Dilepton
        gen = setup_gen_candidates(df)
        tags = ['stat1','lhe']
        if is_lo_w(dataset) or is_nlo_w(dataset) or is_lo_z(dataset) or is_nlo_z(dataset):
            dressed = setup_dressed_gen_candidates(df)
            fill_gen_v_info(df, gen, dressed)
            tags.extend(['dress','combined'])
        elif is_lo_g(dataset) or is_nlo_g(dataset) or is_lo_g_ewk(dataset) or is_nlo_g_ewk(dataset):
            photons = gen[(gen.status==1)&(gen.pdg==22)]
            df['gen_v_pt_stat1'] = photons.pt.max()
            df['gen_v_phi_stat1'] = photons[photons.pt.argmax()].phi.max()
            df['gen_v_pt_lhe'] = df['LHE_Vpt']
            df['gen_v_phi_lhe'] = np.zeros(df.size)

        dijet = genjets[:,:2].distincts()
        mjj = dijet.mass.max()
        nominal = df['Generator_weight']
        for tag in tags:
            sel = selection(df[f'gen_v_phi_{tag}'], genjets, genjets_ak8, lhe_charged_leps)
            for region, requirements in regions.items():
                mask = sel.all(*requirements)
                output[f'gen_vpt_{tag}'].fill(
                                        dataset=dataset,
                                        vpt=df[f'gen_v_pt_{tag}'][mask],
                                        region=region,
                                        weight=nominal[mask]
                                        )


        # Keep track of weight sum
        output['sumw'][dataset] +=  df['genEventSumw']
        output['sumw2'][dataset] +=  df['genEventSumw2']
        return output

    def postprocess(self, accumulator):
        return accumulator
