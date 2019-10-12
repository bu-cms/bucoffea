#!/usr/bin/env python

import numpy as np
from awkward import JaggedArray
from coffea.analysis_objects import JaggedCandidateArray


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


def setup_gen_candidates(df):
    # Find first ancestor with different PDG ID
    # before defining the gen candidates
    mothers = JaggedArray.fromcounts(df['nGenPart'],df['GenPart_genPartIdxMother'] )
    pdgids = JaggedArray.fromcounts(df['nGenPart'],df['GenPart_pdgId'] )
    # parent_index =  find_first_parent(mothers, pdgids)

    gen = JaggedCandidateArray.candidatesfromcounts(
        df['nGenPart'],
        pt=df['GenPart_pt'],
        eta=df['GenPart_eta'],
        phi=df['GenPart_phi'],
        mass=df['GenPart_mass'],
        charge=df['GenPart_pdgId'],
        pdg=df['GenPart_pdgId'],
        status=df['GenPart_status'],
        flag = df['GenPart_statusFlags'],
        mother = df['GenPart_genPartIdxMother'])
    return gen

def islep(pdg):
    """Returns True if the PDG ID represents a lepton."""
    abspdg = np.abs(pdg)
    return (11<=abspdg) & (abspdg<=16)
