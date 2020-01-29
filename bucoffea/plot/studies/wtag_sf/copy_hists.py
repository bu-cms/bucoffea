# Copy TH1Fs from here to data/sf/wtag/

import ROOT

input_filename = 'wtag_mistag_SF.root'
output_filename = '../../../data/sf/ak8/wtag_mistag_SF.root'

input_file = ROOT.TFile.Open(input_filename)
output_file = ROOT.TFile.Open(output_filename, "recreate")

for year in [2017,2018]:
    for wp in ['loose','loosemd','tightmd','tight']:
        hist_name = f'Wmistag_{year}_{wp}_ak8_pt'
        hist = input_file.Get(hist_name)
        output_file.cd()
        hist.Write()

input_file.Close()
output_file.Close()
