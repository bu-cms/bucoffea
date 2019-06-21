import numpy as np

def dphi_jet_met(jets, met_phi, njet=4, ptmin=30):
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

