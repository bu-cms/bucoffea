from coffea import hist
Hist = hist.Hist
Bin = hist.Bin
Cat = hist.Cat

from coffea.analysis_objects import JaggedCandidateArray
import coffea.processor as processor
from awkward import JaggedArray
import numpy as np
from bucoffea.helpers import object_overlap
from bucoffea.helpers.paths import bucoffea_path
from bucoffea.helpers.gen import find_first_parent


from coffea.processor import LazyDataFrame

import uproot

fn = '/eos/cms/store/group/phys_exotica/monojet/aalbert/nanopost/16Jul19/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/DYJetsToLL_M-50_HT-400to600-MLM_2017/190717_212115/0000/tree_6.root'

file = uproot.open(fn)

tree = file['Events']
df = LazyDataFrame(tree, flatten=True)

ak4 = JaggedCandidateArray.candidatesfromcounts(
    df['nJet'],
    pt=df['Jet_pt'],
    eta=df['Jet_eta'],
    phi=df['Jet_phi'],
    mass=df['Jet_mass'],
)
muons = JaggedCandidateArray.candidatesfromcounts(
    df['nMuon'],
    pt=df['Muon_pt'],
    eta=df['Muon_eta'],
    phi=df['Muon_phi'],
    mass= 0 * df['Muon_pt'],
)

