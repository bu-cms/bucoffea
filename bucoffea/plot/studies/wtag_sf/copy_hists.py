# Copy TH1Fs from here to data/sf/wtag/

import ROOT

input_filename = 'wtag_mistag_SF.root'
output_filename = '../../../../data/sf/ak8/wtag_mistag_SF.root'

input_file = ROOT.TFile.Open(input_filename)
output_file = ROOT.TFile.Open(output_filename, "recreate")

input_file.ls()

for year in [2017,2018]:
    for wp in ['loose','loosemd','tightmd','tight']:
        for lepton_tag in ['g','wz']:
            hist_name = f'mistag_SF_{lepton_tag}_{wp}_{year}'
            hist = input_file.Get(hist_name)
            output_file.cd()
            try:
                hist.Write()
            except ReferenceError:
                print(f"histogram {hist_name} not found")
                pass

input_file.Close()
output_file.Close()
