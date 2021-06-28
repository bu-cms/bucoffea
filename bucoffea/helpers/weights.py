import re

import coffea.processor as processor
import numpy as np

from coffea.btag_tools.btagscalefactor import BTagScaleFactor
from bucoffea.helpers.dataset import extract_year
from bucoffea.helpers.gen import get_gen_photon_pt
from bucoffea.helpers.paths import bucoffea_path

def gen_check_for_leptons(leptons, veto_weights, tau=False):
    '''
    Return the veto weights after checking if the leptons in the event are matched to a proper GEN-level lepton.
    For gen-matching on taus, use tau=True, otherwise use tau=False.
    '''
    # For muons and electrons, we require the gen particle flavor to be non-zero
    if not tau:
        gen_match_ok = (leptons.counts == 0) | ((leptons.genpartflav != 0).all() )
    # For taus, we require the gen particle flavor to be exactly 5
    else:
        gen_match_ok = (leptons.counts == 0) | ((leptons.genpartflav == 5).all() )

    # If an event does not pass lepton gen-matching, assign a weight of 0
    new_veto_weights = np.where(gen_match_ok, veto_weights, 0.)

    return new_veto_weights


def get_veto_weights(df, cfg, evaluator, electrons, muons, taus, do_variations=False):
    """
    Calculate veto weights for SR W

    The weights are effectively:

        w = product(1-SF)

    where the product runs overveto-able e, mu, tau.
    """
    veto_weights = processor.Weights(size=df.size, storeIndividual=True)

    variations = ["nominal"]
    if do_variations:
        variations.extend([
                      'ele_reco_up','ele_reco_dn',
                      'ele_id_up','ele_id_dn',
                      'muon_id_up','muon_id_dn',
                      'muon_iso_up','muon_iso_dn',
                      'tau_id_up','tau_id_dn'
                      ])

    for variation in variations:
        def varied_weight(sfname, *args):
            '''Helper function to easily get the correct weights for a given variation'''

            # For the nominal variation, just pass through
            if 'nominal' in variation:
                return evaluator[sfname](*args)

            # If this variation is unrelated to the SF at hand,
            # pass through as well
            if not (re.sub('_(up|dn)', '', variation) in sfname):
                return evaluator[sfname](*args)

            # Direction of variation
            sgn = 1 if variation.endswith("up") else -1
            return evaluator[sfname](*args) + sgn * evaluator[f"{sfname}_error"](*args)


        ### Electrons (For UL: Both 2017 and 2018 have their SFs split by electron pt)
        if extract_year(df['dataset']) == 2017 or cfg.RUN.ULEGACYV8:
            high_et = electrons.pt>20

            # Low pt SFs
            low_pt_args = (electrons.etasc[~high_et], electrons.pt[~high_et])
            ele_reco_sf_low = varied_weight('ele_reco_pt_lt_20', *low_pt_args)
            ele_id_sf_low = varied_weight("ele_id_loose", *low_pt_args)

            # High pt SFs
            high_pt_args = (electrons.etasc[high_et], electrons.pt[high_et])

            ele_reco_sf_high = varied_weight("ele_reco", *high_pt_args)
            ele_id_sf_high = varied_weight("ele_id_loose", *high_pt_args)

            # Combine
            veto_weight_ele = (1 - ele_reco_sf_low*ele_id_sf_low).prod() * (1-ele_reco_sf_high*ele_id_sf_high).prod()
        else:
            # No split for 2018
            args = (electrons.etasc, electrons.pt)
            ele_reco_sf = varied_weight("ele_reco", *args)
            ele_id_sf = varied_weight("ele_id_loose", *args)

            # Combine
            veto_weight_ele = (1 - ele_id_sf*ele_reco_sf).prod()

        # Gen-checking for electrons
        if cfg.ELECTRON.GENCHECK:
            veto_weight_ele = gen_check_for_leptons(electrons, veto_weight_ele)
        
        ### Muons (eta, pT order different for EOY and UL)
        if cfg.RUN.ULEGACYV8:
            args = (muons.abseta, muons.pt)
        else:
            args = (muons.pt, muons.abseta)
        veto_weight_muo = (1 - varied_weight("muon_id_loose", *args)*varied_weight("muon_iso_loose", *args)).prod()

        # Gen-checking for muons
        if cfg.MUON.GENCHECK:
            veto_weight_muo = gen_check_for_leptons(muons, veto_weight_muo)

        ### Taus
        # Taus have their variations saves as separate histograms,
        # so our cool trick from above is replaced by the pedestrian way
        if "tau_id" in variation:
            direction = variation.split("_")[-1]
            tau_sf_name = f"tau_id_{direction}"
        else:
            tau_sf_name = "tau_id"
        veto_weight_tau = (1 - evaluator[tau_sf_name](taus.pt)).prod()

        # If event has at least one tau and it does NOT match to a gen-level tau 
        # with Tau_genPartFlav == 5, assign a weight of 0 and discard that event
        # Right now we're only doing this for VBF (specified in config)
        veto_weight_tau = gen_check_for_leptons(taus, veto_weight_tau, tau=True)

        ### Combine
        total = veto_weight_ele * veto_weight_muo * veto_weight_tau

        # Cap weights just in case
        total[np.abs(total)>5] = 1
        veto_weights.add(variation, total)

    return veto_weights

def diboson_nlo_weights(df, evaluator, gen):

    if re.match('WW_(PSweights_)?201(6|7|8)',df['dataset']):
        gen_wp = gen[(gen.status==62) & (gen.pdg==24)]
        pt = gen_wp.pt.max()
        w = evaluator['total_nlo_ww_ptwp'](pt)
        dw = evaluator['total_nlo_ww_ptwp_error'](pt)

    elif re.match('WZ_(PSweights_)?201(6|7|8)',df['dataset']):
        gen_wp = gen[(gen.status==62) & (gen.pdg==24)]
        gen_wm = gen[(gen.status==62) & (gen.pdg==-24)]

        pt1 = gen_wp.pt.max()
        pt2 = gen_wm.pt.max()

        w1 = evaluator['total_nlo_wpz_ptwp'](pt1)
        w2 = evaluator['total_nlo_wmz_ptwm'](pt2)
        dw1 = evaluator['total_nlo_wpz_ptwp_error'](pt1)
        dw2 = evaluator['total_nlo_wmz_ptwm_error'](pt2)
        pt = np.where(pt1>0,pt1,pt2)
        w = np.where(
            pt1 > 0,
            w1,
            w2
        )
        dw = np.where(
            pt1 > 0,
            dw1,
            dw2
        )
    elif re.match('ZZ_(PSweights_)?201(6|7|8)',df['dataset']):
        gen_z = gen[(gen.status==62) & (np.abs(gen.pdg)==23)]
        pt = gen_z.pt.max()
        w = evaluator['total_nlo_zz_ptz'](pt)
        dw = evaluator['total_nlo_zz_ptz_error'](pt)
    elif re.match("WQQGamma_5f_NLO_FXFX-amcatnlo_201(6|7|8)", df['dataset']):
        pt = get_gen_photon_pt(gen)
        w = evaluator['total_nlo_wg_ptg'](pt)
        dw = evaluator['total_nlo_wg_ptg_error'](pt)
    elif re.match("ZQQGamma_5f_NLO_FXFX-amcatnlo_201(6|7|8)", df['dataset']):
        pt = get_gen_photon_pt(gen)
        w = evaluator['total_nlo_zg_ptg'](pt)
        dw = evaluator['total_nlo_zg_ptg_error'](pt)
    else:
        w = np.ones(df.size)
        pt = w
        dw = np.zeros(df.size)

    w[pt<0] = 1
    dw[pt<0] = 0
    df['weight_diboson_nlo'] = w
    df['weight_diboson_nlo_rel_unc'] = dw/w

def btag_weights(bjets, cfg):
    # Evaluate weight variations
    weight_variations = {}

    # Only calculate for DeepCSV
    if cfg.BTAG.ALGO != "deepcsv":
        weight_variations["central"] = bjets.pt.ones_like()
        weight_variations["up"] = bjets.pt.ones_like()
        weight_variations["down"] = bjets.pt.ones_like()
        return weight_variations

    # Heavy lifting done by coffea implementation
    bsf = BTagScaleFactor(
                          filename=bucoffea_path(cfg.SF.DEEPCSV.FILE),
                          workingpoint=cfg.BTAG.WP.upper(),
                          methods='comb,comb,incl' # Comb for b and c flavors, incl for light
                          )


    for variation in ["central","up","down"]:
        # Use unsmeared jet pt while calculating the b-weights
        weights = bsf.eval(
                        systematic=variation,
                        flavor=bjets.hadflav,
                        abseta=bjets.abseta,
                        pt=bjets.pt / bjets.jercorr,
                        ignore_missing=True)

        # Cap the weights just in case
        weights[np.abs(weights)>5] = 1

        weight_variations[variation] = weights

    return weight_variations

def get_varied_ele_sf(electrons, df, evaluator):
    '''Electron ID and RECO scale factors with variations.'''
    # Ignore the electrons in the gap region
    mask_electron_nogap = (np.abs(electrons.etasc)<1.4442) | (np.abs(electrons.etasc)>1.566)
    electrons_nogap = electrons[mask_electron_nogap]
    electron_is_tight_electron = df['is_tight_electron'][mask_electron_nogap]

    electrons_nogap_tight = electrons_nogap[ electron_is_tight_electron]
    electrons_nogap_loose = electrons_nogap[~electron_is_tight_electron]
    eletight_id_sf = {
        "up":   evaluator['ele_id_tight'](electrons_nogap_tight.etasc, electrons_nogap_tight.pt) + evaluator['ele_id_tight_error'](electrons_nogap_tight.etasc, electrons_nogap_tight.pt),
        "down": evaluator['ele_id_tight'](electrons_nogap_tight.etasc, electrons_nogap_tight.pt) - evaluator['ele_id_tight_error'](electrons_nogap_tight.etasc, electrons_nogap_tight.pt),
        "nom":  evaluator['ele_id_tight'](electrons_nogap_tight.etasc, electrons_nogap_tight.pt)
    }
    eleloose_id_sf = {
        "up":   evaluator['ele_id_loose'](electrons_nogap_loose.etasc, electrons_nogap_loose.pt) + evaluator['ele_id_loose_error'](electrons_nogap_loose.etasc, electrons_nogap_loose.pt),
        "down": evaluator['ele_id_loose'](electrons_nogap_loose.etasc, electrons_nogap_loose.pt) - evaluator['ele_id_loose_error'](electrons_nogap_loose.etasc, electrons_nogap_loose.pt),
        "nom":  evaluator['ele_id_loose'](electrons_nogap_loose.etasc, electrons_nogap_loose.pt)
    }

    high_et = electrons_nogap.pt>20
    ele_reco_sf = evaluator['ele_reco'](electrons_nogap.etasc[high_et], electrons_nogap.pt[high_et]).prod() * \
        evaluator['ele_reco_pt_lt_20'](electrons_nogap.etasc[~high_et], electrons_nogap.pt[~high_et]).prod()

    ele_reco_sf_up = (evaluator['ele_reco'](electrons_nogap.etasc[high_et], electrons_nogap.pt[high_et]) + \
        evaluator['ele_reco_error'](electrons_nogap.etasc[high_et], electrons_nogap.pt[high_et])).prod() * \
        (evaluator['ele_reco_pt_lt_20'](electrons_nogap.etasc[~high_et], electrons_nogap.pt[~high_et]) + \
        evaluator['ele_reco_pt_lt_20_error'](electrons_nogap.etasc[~high_et], electrons_nogap.pt[~high_et])).prod()

    ele_reco_sf_down = (evaluator['ele_reco'](electrons_nogap.etasc[high_et], electrons_nogap.pt[high_et]) - \
        evaluator['ele_reco_error'](electrons_nogap.etasc[high_et], electrons_nogap.pt[high_et])).prod() * \
        (evaluator['ele_reco_pt_lt_20'](electrons_nogap.etasc[~high_et], electrons_nogap.pt[~high_et]) - \
        evaluator['ele_reco_pt_lt_20_error'](electrons_nogap.etasc[~high_et], electrons_nogap.pt[~high_et])).prod()

    ele_reco_sf =  {
        "up" :  ele_reco_sf_up,
        "down" :  ele_reco_sf_down,
        "nom" :  ele_reco_sf,
    }

    return eleloose_id_sf, eletight_id_sf, ele_reco_sf
