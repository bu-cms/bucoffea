import ROOT
from array import array

outfilename = "../../../data/sf/ak8/wtag_eff_SF.root"
#outfilename = "wtag_eff_SF_2017test.root"
outfile = ROOT.TFile.Open(outfilename, "recreate")

#ref: These numbers are read from https://indico.cern.ch/event/840827/contributions/3527925/attachments/1895214/3126510/DeepAK8_Top_W_SFs_2017_JMAR_PK.pdf page23
#loose: 5%, tight: 0.5%
binning    = [200, 300, 400, 800]
SFs = {}
SFs[(2017,'loose'  )] = [0.93, 1.00, 1.06]
SFs[(2017,'tight'  )] = [0.83, 0.79, 0.80]
SFs[(2018,'loose'  )] = [1.05, 0.82, 0.89]
SFs[(2018,'tight'  )] = [0.75, 0.76, 0.74]
SFs[(2017,'loosemd')] = [0.95, 0.99, 0.87]
SFs[(2017,'tightmd')] = [0.80, 0.66, 0.79]
SFs[(2018,'loosemd')] = [0.93, 0.95, 0.94]
SFs[(2018,'tightmd')] = [0.76, 0.71, 0.75]

for year in [2017,2018]:
    for wp in ['loose','loosemd','tight','tightmd']:
        htmp = ROOT.TH1F(f'WTag_{year}_{wp}_ak8_pt', f'WTag_{year}_{wp}_ak8_pt', len(binning)-1, array('d',binning))
        for ibin, SF in enumerate(SFs[(year, wp)]):
            htmp.SetBinContent(ibin+1, SF)
        htmp.Write()

outfile.Close()
