#!/usr/bin/env python

import numpy as np
from awkward import JaggedArray
from coffea.analysis_objects import JaggedCandidateArray
from bucoffea.helpers.dataset import (
                                      is_lo_w,
                                      is_lo_w_ewk,
                                      is_lo_z,
                                      is_lo_z_ewk,
                                      is_nlo_w,
                                      is_nlo_z
                                      )

def find_first_parent(in_mother, in_pdg, maxgen=10):
    """Finds the first parent with a PDG ID different from the daughter

    :param in_mother: Index of the mother particle for each gen. particle
    :type in_mother: JaggedArray
    :param pdg: PDG ID for each gen. particle
    :type pdg: JaggedArray
    :param maxgen: Number of maximal generations to go back, defaults to 10
    :type maxgen: int, optional
    :return: Index and PDG id of first parent with diff. PDG ID
    :rtype: tuple of JaggedArrays
    """
    out_mother = (in_mother>=0) * in_mother # Index of parent, output
    tmp_mother = out_mother # Index of parent, working copy
    found = np.zeros(in_pdg.size) # Mask

    # Loop over generations
    # If mother particle ID is same as daughter particle ID, do nothing
    # Otherwise, update parent_id, but only once!
    for i in range(maxgen):
        # Make sure we dont go negative
        tmp_mother = (tmp_mother>=0) * tmp_mother
        update = (in_pdg[tmp_mother]!=in_pdg) * (found==0)
        out_mother = update*tmp_mother + (~update) * out_mother
        found=found+update
        tmp_mother = in_mother[tmp_mother]

    return out_mother

def find_gen_dilepton(gen, pdgsum=0):
    """
    Finds and builds dilepton candidate from gen particles.

    Dilepton candidates are constructed from all charged and
    neutral leptons. A dilepton candidate is considered valid if
    the absolute value of the sum of the PDG IDs of the constituents
    is equal to the pdgsum parameter, and both constituents can be
    traced back to the same parent. It is *not* checked what exactly
    the parents are, as there will be no record for off-shell bosons.

    Choose pdgsum=0 for Z candidates, pdgsum=1 for W candidates.

    :param gen: Gen candidates
    :type gen: JaggedCandidateArray
    :param pdgsum: Absolute sum of PDG IDs to form a valid candidate
    :type pdgsum: int
    :return: Dilepton candidates
    :rtype: JaggedCandidateArray
    """
    leps = gen[(((gen.status==1) & islep(gen.pdg))) | ((gen.status==2) & (np.abs(gen.pdg)==15))]
    dileps = leps.distincts()

    dilepton_flavour = np.abs(dileps.i0.pdg + dileps.i1.pdg) == pdgsum
    good_dileps = dileps[dilepton_flavour]
    return good_dileps


def stat1_dilepton(df, gen):
    """Build a dilepton candidate from status 1 leptons

    :param df: Data frame
    :type df: dataframe
    :param gen: Gen. candidates
    :type gen: JaggedCandidateArray
    :return: pt and phi of dilepton
    :rtype: tuple of two 1D arrays
    """
    if is_lo_z(df['dataset']) or is_lo_z_ewk(df['dataset']) or is_nlo_z(df['dataset']):
        pdgsum = 0
    elif is_lo_w(df['dataset']) or  is_lo_w_ewk(df['dataset']) or is_nlo_w(df['dataset']):
        pdgsum = 1
    gen_dilep = find_gen_dilepton(gen, pdgsum)
    gen_dilep = gen_dilep[gen_dilep.mass.argmax()]
    return gen_dilep.pt.max(), gen_dilep.phi.max()


def merge_dileptons(dilepton1, dilepton2, dilepton3=None):
    """
    Choose highest mass dilepton from up to three option lists.

    :return: pt and phi of the chosen dilepton
    :rtype: tuple of two 1D arrays
    """

    mmax1 = dilepton1.mass.max()
    mmax2 = dilepton2.mass.max()
    if dilepton3 is not None:
        mmax3 = dilepton3.mass.max()
    else:
        mmax3 = -1 * np.ones(dilepton1.size)

    take2 = (mmax2 > mmax1) & (mmax2 > mmax3)
    take3 = (mmax3 > mmax1) & (mmax3 > mmax2)
    take1 = ~(take2 | take3)

    vpt1 = dilepton1.pt.max()
    vphi1 = dilepton1.phi.max()
    vpt1[~take1] = 0
    vphi1[~take1] = 0

    vpt2 = dilepton2.pt.max()
    vphi2 = dilepton2.phi.max()
    vpt2[~take2] = 0
    vphi2[~take2] = 0

    if dilepton3 is not None:
        vpt3 = dilepton3.pt.max()
        vphi3 = dilepton3.phi.max()
        vpt3[~take3] = 0
        vphi3[~take3] = 0
    else:
        vpt3 = np.zeros(dilepton1.size)
        vphi3 = np.zeros(dilepton1.size)
    vphi = vphi1 + vphi2 + vphi3
    vpt = vpt1 + vpt2 + vpt3

    return vpt, vphi


def dressed_dilep(df, gen, dressed):
    """
    Build a dilepton candidate from dressed leptons.

    :param df: Data frame
    :type df: dataframe
    :param gen: Gen. candidates
    :type gen: JaggedCandidateArray
    :param dressed: Dressed gen candidates
    :type dressed: JaggedCandidateArray
    :return: pt and phi of dilepton
    :rtype: tuple of two 1D arrays
    """
    # Dressed leptons
    neutrinos = gen[(gen.status==1) & isnu(gen.pdg)]
    if is_lo_z(df['dataset']) or is_nlo_z(df['dataset']) or is_lo_z_ewk(df['dataset']):
        # e, mu
        dilep_dress = find_gen_dilepton(dressed, 0)
        dilep_dress = dilep_dress[dilep_dress.mass.argmax()]

        #nu
        dilep_nu = find_gen_dilepton(neutrinos, 0)
        dilep_nu = dilep_nu[dilep_nu.mass.argmax()]

        # tau
        dilep_tau = find_gen_dilepton(gen[np.abs(gen.pdg)==15], 0)
        dilep_tau = dilep_tau[dilep_tau.mass.argmax()]

        # Merge by taking higher-mass
        return merge_dileptons(dilep_tau, dilep_dress, dilep_nu)

    elif is_lo_w(df['dataset']) or is_nlo_w(df['dataset']) or is_lo_w_ewk(df['dataset']):
        # e, mu
        dilep_dress = dressed.cross(neutrinos)
        dilep_dress = dilep_dress[np.abs(dilep_dress.i0.pdg + dilep_dress.i1.pdg)==1]
        dilep_dress = dilep_dress[dilep_dress.mass.argmax()]

        # tau
        dilep_tau = find_gen_dilepton(gen[(np.abs(gen.pdg)==15) | (np.abs(gen.pdg)==16)], 1)
        dilep_tau = dilep_tau[dilep_tau.mass.argmax()]
        return  merge_dileptons(dilep_tau, dilep_dress)

def fill_gen_v_info(df, gen, dressed):
    '''
    One-stop function to generate gen v pt info.

    For stat1, dressed and lhe V, the pt and phi
    information is written into the data frame.
    '''

    df['gen_v_pt_stat1'], df['gen_v_phi_stat1'] = stat1_dilepton(df, gen)
    df['gen_v_pt_dress'], df['gen_v_phi_dress'] = dressed_dilep(df, gen, dressed)

    # For events where we cannot find a dressed dilepton,
    # fill with stat 1 dileptons
    invalid = df['gen_v_pt_dress'] <= 0
    df['gen_v_pt_dress'][invalid] = df['gen_v_pt_stat1'][invalid]
    df['gen_v_phi_dress'][invalid] = df['gen_v_phi_stat1'][invalid]

    # LHE is just pass through
    df['gen_v_pt_lhe'] = df['LHE_Vpt']
    df['gen_v_phi_lhe'] = np.zeros(df.size)

def setup_gen_candidates(df):
    gen = JaggedCandidateArray.candidatesfromcounts(
        df['nGenPart'],
        pt=df['GenPart_pt'],
        eta=df['GenPart_eta'],
        phi=df['GenPart_phi'],
        mass=df['GenPart_mass'],
        charge=df['GenPart_pdgId'],
        pdg=df['GenPart_pdgId'],
        status=df['GenPart_status'],
        flag = df['GenPart_statusFlags'])
    return gen

def setup_gen_jets(df):
    genjets = JaggedCandidateArray.candidatesfromcounts(
        df['nGenJet'],
        pt=df['GenJet_pt'],
        eta=df['GenJet_eta'],
        phi=df['GenJet_phi'],
        mass=0*df['GenJet_pt']
        )
    return genjets

def setup_dressed_gen_candidates(df):
    dressed = JaggedCandidateArray.candidatesfromcounts(
        df['nGenDressedLepton'],
        pt=df['GenDressedLepton_pt'],
        eta=df['GenDressedLepton_eta'],
        phi=df['GenDressedLepton_phi'],
        mass=0*df['GenDressedLepton_pt'],
        status=np.ones(df['GenDressedLepton_pt'].size),
        pdg=df['GenDressedLepton_pdgId'])
    return dressed

def islep(pdg):
    """Returns True if the PDG ID represents a lepton."""
    abspdg = np.abs(pdg)
    return (11<=abspdg) & (abspdg<=16)

def setup_lhe_cleaned_genjets(df):
    genjets = JaggedCandidateArray.candidatesfromcounts(
            df['nGenJet'],
            pt=df['GenJet_pt'],
            eta=df['GenJet_eta'],
            abseta=np.abs(df['GenJet_eta']),
            phi=df['GenJet_phi'],
            mass=df['GenJet_mass']
        )
    lhe = JaggedCandidateArray.candidatesfromcounts(
                df['nLHEPart'],
                pt=df['LHEPart_pt'],
                eta=df['LHEPart_eta'],
                phi=df['LHEPart_phi'],
                mass=df['LHEPart_mass'],
                pdg=df['LHEPart_pdgId'],
            )

    lhe_leps_gams = lhe[(islep(lhe.pdg)) | (lhe.pdg==22)]

    return genjets[(~genjets.match(lhe_leps_gams,deltaRCut=0.4))]


def isnu(pdg):
    """Returns True if the PDG ID represents a neutrino."""
    abspdg = np.abs(pdg)
    return (12==abspdg) | (14==abspdg) | (16==abspdg)
