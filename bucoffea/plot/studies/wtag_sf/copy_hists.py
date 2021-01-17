# Copy TH1Fs from here to data/sf/wtag/

import ROOT

# use single bin for tight SF and multibin for loose SF
input_filename = {  'loose' : 'output_mistag/merged_2020_11_23_03Sep20v7_Wmass65To120/output_bin5_massden_massnum_nlogjet_splitsys_realVSF1p0/wtag_mistag_SF.root',
                    'tight' : 'output_mistag/merged_2020_11_23_03Sep20v7_Wmass65To120/output_bin2_massden_massnum_nlogjet_splitsys_realVSF1p0/wtag_mistag_SF.root',
        }
output_filename = '../../../data/sf/ak8/wtag_mistag_SF.root'

output_file = ROOT.TFile.Open(output_filename, "recreate")


for wp in ['loose','tight']:
    input_file = ROOT.TFile.Open(input_filename[wp])
    for year in [2017,2018]:
        for lepton_tag in ['g','w','z']:
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
print(f"copied histograms from file {input_filename} to {output_filename}")
