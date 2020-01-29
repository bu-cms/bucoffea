import ROOT
from array import array

outfilename = "../../../data/sf/ak8/wtag_eff_SF.root"
#outfilename = "wtag_eff_SF_2017test.root"
outfile = ROOT.TFile.Open(outfilename, "recreate")

#ref: These numbers are read from https://www.dropbox.com/s/crezbefq0lbqfkp/DeepAK8V2_Top_W_SFs.csv?dl=0
#loose: 5%, tight: 0.5%
binning    = [200, 300, 400, 800]
SFs = {}
SFs[(2017,'tight'  )] = [ [0.94 , 0.05        , 0.04],
                          [0.69 , 0.15        , 0.09],
                          [0.96 , 0.07        , 0.04],]
SFs[(2017,'tightmd')] = [ [0.86 , 0.06        , 0.06],
                          [0.76 , 0.05        , 0.05],
                          [0.72 , 0.06        , 0.06],]
SFs[(2017,'loose'  )] = [ [0.93 , 0.07        , 0.07],
                          [0.91 , 0.06        , 0.06],
                          [1.07 , 0.08        , 0.04],]
SFs[(2017,'loosemd')] = [ [1.02 , 0.03        , 0.03],
                          [0.89 , 0.03        , 0.02],
                          [0.96 , 0.04        , 0.04],]
SFs[(2018,'tight'  )] = [ [0.94 , 0.07        , 0.05],
                          [0.84 , 0.06        , 0.06],
                          [0.97 , 0.05        , 0.03],]
SFs[(2018,'tightmd')] = [ [0.76 , 0.05        , 0.05],
                          [0.67 , 0.04        , 0.04],
                          [0.75 , 0.05        , 0.05],]
SFs[(2018,'loose'  )] = [ [1.04 , 0.05        , 0.10],
                          [0.85 , 0.05        , 0.05],
                          [0.95 , 0.06        , 0.06],]
SFs[(2018,'loosemd')] = [ [0.95 , 0.03        , 0.02],
                          [0.92 , 0.02        , 0.02],
                          [0.95 , 0.03        , 0.03],]

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
