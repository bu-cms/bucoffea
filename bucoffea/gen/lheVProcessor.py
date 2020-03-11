import coffea.processor as processor
import numpy as np
from coffea import hist

from bucoffea.helpers import min_dphi_jet_met, dphi
from bucoffea.helpers.dataset import (is_lo_g, is_lo_g_ewk, is_lo_w, is_lo_z,
                                      is_nlo_g,is_nlo_g_ewk, is_nlo_w, is_nlo_z,)
from bucoffea.helpers.gen import (fill_gen_v_info,
                                  setup_dressed_gen_candidates,
                                  setup_gen_candidates,
                                  setup_lhe_cleaned_genjets,
                                  setup_lhe_parton_photon_pairs)

Hist = hist.Hist
Bin = hist.Bin
Cat = hist.Cat

def vbf_selection(vphi, dijet, genjets):
    selection = processor.PackedSelection()
    selection.add(
                  'two_jets',
                  dijet.counts>0
                  )
    selection.add(
                  'leadak4_pt_eta',
                  (dijet.i0.pt.max() > 80) & (np.abs(dijet.i0.eta.max()) < 5.0)
                  )
    selection.add(
                  'trailak4_pt_eta',
                  (dijet.i1.pt.max() > 40) & (np.abs(dijet.i1.eta.max()) < 5.0)
                  )
    selection.add(
                  'hemisphere',
                  (dijet.i0.eta.max()*dijet.i1.eta.max() < 0)
                  )
    selection.add(
                  'mindphijr',
                  min_dphi_jet_met(genjets, vphi, njet=4, ptmin=30, etamax=5.0) > 0.5
                  )
    selection.add(
                  'detajj',
                  np.abs(dijet.i0.eta-dijet.i1.eta).max() > 1
                  )
    selection.add(
                  'dphijj',
                  dphi(dijet.i0.phi,dijet.i1.phi).min() < 1.5
                  )

    return selection

def monojet_selection(vphi, genjets):
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

    return selection


class lheVProcessor(processor.ProcessorABC):
    def __init__(self):

        # Histogram setup
        dataset_ax = Cat("dataset", "Primary dataset")

        vpt_ax = Bin("vpt",r"$p_{T}^{V}$ (GeV)", 50, 0, 2000)
        jpt_ax = Bin("jpt",r"$p_{T}^{j}$ (GeV)", 50, 0, 2000)
        mjj_ax = Bin("mjj",r"$m(jj)$ (GeV)", 75, 0, 7500)
        res_ax = Bin("res",r"pt: dressed / stat1 - 1", 80,-0.2,0.2)
        dr_ax = Bin("dr", r"$\Delta R$", 100, 0, 10)

        items = {}
        for tag in ['stat1','dress','lhe','combined']:
            items[f"gen_vpt_inclusive_{tag}"] = Hist("Counts",
                                    dataset_ax,
                                    vpt_ax)
            items[f"gen_vpt_monojet_{tag}"] = Hist("Counts",
                                    dataset_ax,
                                    jpt_ax,
                                    vpt_ax)
            items[f"gen_vpt_vbf_{tag}"] = Hist("Counts",
                                    dataset_ax,
                                    jpt_ax,
                                    mjj_ax,
                                    vpt_ax)
            items[f"gen_vpt_vbf_{tag}_withDRreq"] = Hist("Counts",
                                            dataset_ax,
                                            jpt_ax,
                                            mjj_ax,
                                            vpt_ax)
            items[f'lhe_mindr_g_parton_{tag}'] = Hist("Counts",
                                                dataset_ax,
                                                dr_ax)
            items[f'lhe_mindr_g_parton_{tag}_noDRreq'] = Hist("Counts",
                                    dataset_ax,
                                    dr_ax)

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

            # Get LHE level photon + parton pairs
            # Calculate minimum deltaR between them in each event
            pairs = setup_lhe_parton_photon_pairs(df)
            min_dr = pairs.i0.p4.delta_r(pairs.i1.p4).min()
            df['lhe_mindr_g_parton'] = min_dr

        dijet = genjets[:,:2].distincts()
        mjj = dijet.mass.max()
        for tag in tags:
            # Dijet for VBF

            # Selection
            vbf_sel = vbf_selection(df[f'gen_v_phi_{tag}'], dijet, genjets)
            monojet_sel = monojet_selection(df[f'gen_v_phi_{tag}'], genjets)

            nominal = df['Generator_weight']

            output[f'gen_vpt_inclusive_{tag}'].fill(
                                    dataset=dataset,
                                    vpt=df[f'gen_v_pt_{tag}'],
                                    weight=nominal
                                    )
                                    
            mask_vbf = vbf_sel.all(*vbf_sel.names)
            output[f'gen_vpt_vbf_{tag}'].fill(
                                    dataset=dataset,
                                    vpt=df[f'gen_v_pt_{tag}'][mask_vbf],
                                    jpt=genjets.pt.max()[mask_vbf],
                                    mjj = mjj[mask_vbf],
                                    weight=nominal[mask_vbf]
                                    )

            # Fill the histogram with minimum deltaR between photons
            # and partons at LHE level
            if is_lo_g(dataset) or is_nlo_g(dataset) or is_lo_g_ewk(dataset) or is_nlo_g_ewk(dataset):
                # Add new deltaR requirement:
                # deltaR > 0.4 for every event
                dr_mask = df['lhe_mindr_g_parton'] > 0.4
                full_mask_vbf = mask_vbf*dr_mask

                # Fill the histogram with the deltaR requirement
                output[f'lhe_mindr_g_parton_{tag}'].fill(
                                            dataset=dataset,
                                            dr=df['lhe_mindr_g_parton'][full_mask_vbf],
                                            weight=nominal[full_mask_vbf]
                                            )
                # Fill the histogram without the deltaR requirement
                output[f'lhe_mindr_g_parton_{tag}_noDRreq'].fill(
                                                dataset=dataset,
                                                dr=df['lhe_mindr_g_parton'][mask_vbf],
                                                weight=nominal[mask_vbf]
                                                )
                
                # Fill V-pt and mjj with the deltaR > 0.4 requirement
                output[f'gen_vpt_vbf_{tag}_withDRreq'].fill(
                                                        dataset=dataset,
                                                        vpt=df[f'gen_v_pt_{tag}'][full_mask_vbf],
                                                        jpt=genjets.pt.max()[full_mask_vbf],
                                                        mjj=mjj[full_mask_vbf],
                                                        weight=nominal[full_mask_vbf]
                                                    )
                                    
            mask_monojet = monojet_sel.all(*monojet_sel.names)

            output[f'gen_vpt_monojet_{tag}'].fill(
                                    dataset=dataset,
                                    vpt=df[f'gen_v_pt_{tag}'][mask_monojet],
                                    jpt=genjets.pt.max()[mask_monojet],
                                    weight=nominal[mask_monojet]
                                    )


        # Keep track of weight sum
        output['sumw'][dataset] +=  df['genEventSumw']
        output['sumw2'][dataset] +=  df['genEventSumw2']

        return output

    def postprocess(self, accumulator):
        return accumulator
