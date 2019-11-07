import ROOT
from array import array

outfilename = "../../../data/sf/ak8/wtag_eff_SF_2017.root"
#outfilename = "wtag_eff_SF_2017.root"
outfile = ROOT.TFile.Open(outfilename, "recreate")

#ref: These numbers are read from https://indico.cern.ch/event/840827/contributions/3527925/attachments/1895214/3126510/DeepAK8_Top_W_SFs_2017_JMAR_PK.pdf page23
binning    = [200, 300, 400, 800]
SF_loose   = [1.01, 0.93, 1.02]
SF_loosemd = [0.81, 0.95, 0.85]
SF_tight   = [0.86, 0.86, 0.88]
SF_tightmd = [0.75, 0.61, 0.78]

h_loose = ROOT.TH1F("WTag_loose_ak8_pt", "WTag_loose_ak8_pt", len(binning)-1, array('d',binning))
for ibin, SF in enumerate(SF_loose):
    h_loose.SetBinContent(ibin+1, SF)
h_loosemd = ROOT.TH1F("WTag_loosemd_ak8_pt", "WTag_loosemd_ak8_pt", len(binning)-1, array('d',binning))
for ibin, SF in enumerate(SF_loosemd):
    h_loosemd.SetBinContent(ibin+1, SF)
h_tight = ROOT.TH1F("WTag_tight_ak8_pt", "WTag_tight_ak8_pt", len(binning)-1, array('d',binning))
for ibin, SF in enumerate(SF_tight):
    h_tight.SetBinContent(ibin+1, SF)
h_tightmd = ROOT.TH1F("WTag_tightmd_ak8_pt", "WTag_tightmd_ak8_pt", len(binning)-1, array('d',binning))
for ibin, SF in enumerate(SF_tightmd):
    h_tightmd.SetBinContent(ibin+1, SF)

outfile.Write()
outfile.Close()
