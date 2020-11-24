import uproot
import time
from mistagSF import outdir
import numpy as np
from matplotlib import pyplot as plt
import os

inputfilename = f"{outdir}/wtag_mistag_SF.root"
convertedfilename = f"{outdir}/wtag_mistag_SF_converted.root"
colors={
        '1m':6,
        '2m':7,
        '1e':8,
        '2e':9,
        'g':11,
        'combined':46,
        'wz':46,
        }

#def axis_limit(axis, wp=None, has_mass_cut = False):
#    if axis == "x":
#        return (0,800)
#    if axis == "y":
#        lim = (0,1)
#        if wp == "tight"

def convert_tefficiencies(inputFileName, outputFileName):
    # no need to convert if the converted file is newer than the origin
    if os.path.exists(outputFileName) and os.path.getmtime(outputFileName)>os.path.getmtime(inputFileName):
        print(f"Fount converted file {convertedfilename}")
        return
    import ROOT
    ROOT.gROOT.SetBatch()
    inputFile = ROOT.TFile.Open(inputFileName)
    outputFile = ROOT.TFile.Open(outputFileName,"recreate")
    listOfKeys = inputFile.GetListOfKeys()
    canv = ROOT.TCanvas("canv","canv",800,800)
    for key in listOfKeys:
        hist_name = key.GetName()
        hist = inputFile.Get(hist_name)
        if "TH1" in hist.ClassName():
            out_hist = hist
        elif "TEfficiency"==hist.ClassName():
            #print(f"Converting TEfficiency {hist_name}")
            hist.Draw('AP')
            gtmp = hist.GetPaintedGraph()
            out_hist = hist.GetTotalHistogram().Clone()
            nbins = out_hist.GetNbinsX()
            for ibin in range(1,nbins+1):
                out_hist.SetBinContent(ibin, hist.GetEfficiency(ibin))
                out_hist.SetBinError(ibin, hist.GetEfficiencyErrorUp(ibin)) # assuming rate<<1, treat as symmetric error
        outputFile.cd()
        out_hist.Write(hist_name)
    outputFile.Close()
    inputFile.Close()
    print(f"Converted file saved as {convertedfilename}")

def main():
    convert_tefficiencies(inputfilename, convertedfilename)
    inputfile = uproot.open(convertedfilename)

    # first plot the rate plots
    fig, ax = plt.subplots(1,1, sharex=True)
    for year in [2017,2018]:
        for wp in ['loose', 'tight','medium']:
            # plot the mistag rates 
            for prefix in ['mistag_rate_data','mistag_rate_mc','mistag_SF']:
                def plot_from_tags(lepton_tags):
                    ax.clear()
                    for lepton_tag in lepton_tags:
                        hist_name = f'{prefix}_{lepton_tag}_{wp}_{year}'
                        hist = inputfile[hist_name]
                        edges = hist.edges
                        centers = 0.5*(edges[:-1]+edges[1:])
                        widths = (edges[1:]-edges[:-1])
                        p = ax.errorbar(centers, hist.values, yerr=np.sqrt(hist.variances), capsize=5, markeredgewidth=1, label=f"{lepton_tag} with stat. uncertainties")
                        if ("rate_data" in prefix) or ("SF" in prefix):
                            hist_sysup = inputfile[hist_name+"_sysUp"]
                            hist_sysdn = inputfile[hist_name+"_sysDn"]
                            ax.bar(centers ,hist_sysup.values-hist_sysdn.values, bottom=hist_sysdn.values, width=widths, label=f"{lepton_tag} with syst. variation on top/diboson", color=p[0].get_color(), alpha=0.5)
                    ax.set_xlim(0,900)
                    ax.set_xlabel(r"AK8 Jet $p_T$")
                    ax.set_title(prefix)
                    if "rate" in prefix:
                        if wp=="loose":  ylim = np.array([0, 0.10])
                        if wp=="medium": ylim = np.array([0, 0.02])
                        if wp=="tight":  ylim = np.array([0, 0.01])
                        if "massden" in convertedfilename:
                            ylim = ylim*10
                        ax.set_ylabel("efficiency")
                    else:
                        ylim = np.array([0.0,3.0])
                        ax.set_ylabel("SF")
                    ax.set_ylim(ylim[0],ylim[1])
                    ax.minorticks_on()
                    ax.grid(b=True, which='major', color='#666666', linestyle='-')
                    ax.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
                    ax.legend()
                plot_from_tags(['g','1e','2e','1m','2m'])
                fig.savefig(f'{outdir}/{prefix}_{wp}_{year}.png')
                plot_from_tags(['g','wz'])
                fig.savefig(f'{outdir}/{prefix}_{wp}_gvswz_{year}.png')
                plot_from_tags(['all'])
                fig.savefig(f'{outdir}/{prefix}_{wp}_all_{year}.png')
    print(f'Saved plots in {outdir}')



if __name__ == "__main__":
    main()
exit()



inputfile = uproot.open(inputfilename)
fig, ax = plt.subplots(2,1,sharex=True)
for year in [2017,2018]:
    for wp in ['loose', 'tight','medium']:
        # plot the mistag rates 
        for prefix in ['mistag_rate_data','mistag_rate_mc']:
            canv.Clear()
            canvEmpty=True
            for lepton_flag in ['g','wz']:
                htmp = inputfile[f'{prefix}_{lepton_flag}_{wp}_{year}']
                #htmp.SetLineStyle(2)
                try:
                    htmp.SetTitle(lepton_flag) # for the legend builder
                except:
                    print(inputfilename)
                    inputfile.ls()
                    print(f'{prefix}_{lepton_flag}_{wp}_{year}')
                    continue
                htmp.SetLineColor(colors[lepton_flag])
                htmp.SetLineWidth(3)
                if canvEmpty:
                    htmp.Draw('AP')
                    canv.Update()
                    gtmp=htmp.GetPaintedGraph()
                    print(htmp)
                    while not gtmp:
                        time.sleep(1)
                        gtmp=htmp.GetPaintedGraph()
                    if 'tight' in wp:
                        if 'nomass' in outdir:
                            gtmp.GetYaxis().SetRangeUser(0,0.020)
                        else:
                            gtmp.GetYaxis().SetRangeUser(0,0.20)
                    elif 'medium' in wp:
                        if 'nomass' in outdir:
                            gtmp.GetYaxis().SetRangeUser(0,0.040)
                        else:
                            gtmp.GetYaxis().SetRangeUser(0,0.40)
                    else:
                        if 'nomass' in outdir:
                            gtmp.GetYaxis().SetRangeUser(0,0.15)
                        else:
                            gtmp.GetYaxis().SetRangeUser(0,1)
                    gtmp.GetYaxis().SetTitle('mistagging rate')
                    gtmp.GetXaxis().SetTitle('AK8 Jet p_{T}')
                    canv.Draw()
                    canvEmpty=False
                else:
                    htmp.Draw('P same')
            canv.BuildLegend()
            canv.Update()
            canv.Print(f'{outdir}/{prefix}_{wp}_{year}.png')
        # plot the mistag SF 
        for prefix in ['mistag_SF']:
            canv.Clear()
            canvEmpty=True
            for lepton_flag in ['g','wz']:
                htmp = inputfile.Get(f'{prefix}_{lepton_flag}_{wp}_{year}')
                if not htmp:
                    continue
                #htmp.SetLineStyle(2)
                htmp.SetTitle(lepton_flag) # for the legend builder
                htmp.SetLineColor(colors[lepton_flag])
                htmp.SetLineWidth(3)
                if canvEmpty:
                    if 'tight' in wp:
                        htmp.GetYaxis().SetRangeUser(0,2.5)
                    else:
                        htmp.GetYaxis().SetRangeUser(0,2.5)
                    htmp.GetYaxis().SetTitle('mistagging SF')
                    htmp.GetXaxis().SetTitle('AK8 Jet p_{T}')
                    htmp.Draw('e1')
                    canvEmpty=False
                else:
                    htmp.Draw('e1 same')
            canv.BuildLegend()
            canv.Update()
            canv.Print(f'{outdir}/{prefix}_{wp}_{year}.png')
