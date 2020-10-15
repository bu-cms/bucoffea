import ROOT
from array import array

outfilename = "../../../data/sf/ak8/wtag_eff_SF.root"
#outfilename = "wtag_eff_SF_2017test.root"
outfile = ROOT.TFile.Open(outfilename, "recreate")

#ref: These numbers are read from https://twiki.cern.ch/twiki/bin/viewauth/CMS/DeepAKXTagging#DeepAK8_V2_Working_Points_and_Sc
#loose: 5%, tight: 0.5%, medium: 2.5%
binning    = [200, 300, 400, 800]
SFs = {}
SFs[(2017,'tight'  )] = [ [1.00 , 0.08        , 0.07],
                          [0.78 , 0.10        , 0.11],
                          [0.82 , 0.08        , 0.08],]
SFs[(2017,'tightmd')] = [ [0.90 , 0.07        , 0.07],
                          [0.78 , 0.06        , 0.05],
                          [0.79 , 0.07        , 0.07],]
SFs[(2017,'loose'  )] = [ [1.05 , 0.04        , 0.04],
                          [1.04 , 0.11        , 0.11],
                          [1.09 , 0.08        , 0.08],]
SFs[(2017,'loosemd')] = [ [1.00 , 0.04        , 0.04],
                          [0.94 , 0.03        , 0.03],
                          [0.95 , 0.03        , 0.03],]
SFs[(2017,'medium'  )] = [[1.19 , 0.04        , 0.04],
                          [1.03 , 0.11        , 0.11],
                          [0.95 , 0.08        , 0.08],]
SFs[(2017,'mediummd')] = [[0.94 , 0.05        , 0.05],
                          [0.92 , 0.04        , 0.04],
                          [0.93 , 0.05        , 0.05],]
SFs[(2018,'tight'  )] = [ [0.87 , 0.06        , 0.06],
                          [0.86 , 0.05        , 0.05],
                          [0.85 , 0.06        , 0.05],]
SFs[(2018,'tightmd')] = [ [0.73 , 0.07        , 0.06],
                          [0.67 , 0.06        , 0.06],
                          [0.75 , 0.07        , 0.08],]
SFs[(2018,'loose'  )] = [ [1.04 , 0.02        , 0.02],
                          [1.03 , 0.09        , 0.09],
                          [1.07 , 0.04        , 0.04],]
SFs[(2018,'loosemd')] = [ [0.96 , 0.03        , 0.03],
                          [0.93 , 0.02        , 0.03],
                          [0.96 , 0.03        , 0.03],]
SFs[(2018,'medium'  )] = [[1.10 , 0.09        , 0.07],
                          [0.96 , 0.07        , 0.08],
                          [1.04 , 0.10        , 0.12],]
SFs[(2018,'mediummd')] = [[0.91 , 0.05        , 0.05],
                          [0.86 , 0.03        , 0.04],
                          [0.90 , 0.04        , 0.04],]

for year in [2017,2018]:
    for wp in ['loose','loosemd','tight','tightmd']:
        htmp = ROOT.TH1F(f'WTag_{year}_{wp}_ak8_pt', f'WTag_{year}_{wp}_ak8_pt', len(binning)-1, array('d',binning))
        htmp_up    = ROOT.TH1F(f'WTagUp_{year}_{wp}_ak8_pt', f'WTagUp_{year}_{wp}_ak8_pt', len(binning)-1, array('d',binning))
        htmp_down  = ROOT.TH1F(f'WTagDown_{year}_{wp}_ak8_pt', f'WTagDown_{year}_{wp}_ak8_pt', len(binning)-1, array('d',binning))
        for ibin, SF in enumerate(SFs[(year, wp)]):
            htmp.SetBinContent(ibin+1, SF[0])
            htmp_up.SetBinContent(ibin+1, SF[1])
            htmp_down.SetBinContent(ibin+1, SF[2])
        htmp.Write()
        htmp_up.Write()
        htmp_down.Write()

outfile.Close()
