# Copy TH1Fs from here to data/sf/wtag/

import ROOT,os

def produce_pol1_fitted_hist(raw_hist, nbins=100, xlim=None):
    raw_name = raw_hist.GetName()
    fit_name = raw_name+"_pol1"
    print("########\n\n######### working on: "+fit_name)
    plot_dir_name = "pol1FitPlots"
    if not os.path.exists(plot_dir_name):
        os.makedirs(plot_dir_name)
    raw_nbins = raw_hist.GetNbinsX()
    xlow = raw_hist.GetBinLowEdge(1)
    xhigh = raw_hist.GetBinLowEdge(raw_nbins) + raw_hist.GetBinWidth(raw_nbins)
    if xlim:
        xlow,xhigh = xlim[0],xlim[1]
    canv = ROOT.TCanvas("canv","canv")
    raw_hist.Draw()
    if not raw_nbins==1:
        ftr = raw_hist.Fit("pol1","NS")
        p0 = ftr.Parameter(0)
        p1 = ftr.Parameter(1)
        cor01 = ftr.Correlation(0,1)
        print("p0-p1 correlation: \t  = "+str(cor01))
    else:
        print("Skip fit for single bin SF")
        p0 = raw_hist.GetBinContent(1)
        p1 = 0
    fit_hist = ROOT.TH1F(fit_name,"",nbins,xlow,xhigh)
    for ibin in range(nbins):
        bin_center = xlow + (0.5 + ibin) * (xhigh-xlow) / nbins
        y_center = p0 + p1 * bin_center
        fit_hist.SetBinContent(ibin+1, y_center)
    fit_hist.SetLineColor(ROOT.kBlue)
    fit_hist.Draw("same")
    canv.Print(f"{plot_dir_name}/{fit_name}.png")
    return fit_hist
        


# use single bin for tight SF and multibin for loose SF
input_filename = {  'loose' : 'output_mistag/merged_2020_11_23_03Sep20v7_Wmass65To120/output_bin5_massden_massnum_nlogjet_splitsys_realVSF1p0/wtag_mistag_SF.root',
                    'tight' : 'output_mistag/merged_2020_11_23_03Sep20v7_Wmass65To120/output_bin5_massden_massnum_nlogjet_splitsys_realVSF1p0/wtag_mistag_SF.root',
                    #'tight' : 'output_mistag/merged_2020_11_23_03Sep20v7_Wmass65To120/output_bin2_massden_massnum_nlogjet_splitsys_realVSF1p0/wtag_mistag_SF.root',
        }
output_filename = '../../../data/sf/ak8/wtag_mistag_SF.root'

output_file = ROOT.TFile.Open(output_filename, "recreate")


for wp in ['loose','tight']:
    input_file = ROOT.TFile.Open(input_filename[wp])
    for year in [2017,2018]:
        for lepton_tag in ['g','w','z']:
            hist_name = f'mistag_SF_{lepton_tag}_{wp}_{year}'
            hist = input_file.Get(hist_name)
            fit_hist = produce_pol1_fitted_hist(hist)
            output_file.cd()
            try:
                hist.Write()
                fit_hist.Write()
            except ReferenceError:
                print(f"histogram {hist_name} not found")
                pass
    input_file.Close()


output_file.Close()
print(f"copied histograms from file {input_filename} to {output_filename}")
