import uproot
import time
from mistagSF import outdir
import numpy as np
from matplotlib import pyplot as plt
import os
import pandas
import ROOT

inputfilename = f"{outdir}/wtag_mistag_SF.root"
convertedfilename = f"{outdir}/wtag_mistag_SF_converted.root"
colors={
        '1m'       : '#1f77b4' , 
        '2m'       : '#ff7f0e' , 
        '1e'       : '#2ca02c' , 
        '2e'       : '#d62728' , 
        'g'        : '#9467bd' , 
        'ave'      : '#17becf' , 
        'w'        : '#8c564b' , 
        'z'        : '#e377c2' , 
        'all' : '#7f7f7f' , 
        'wz'       : '#bcbd22' , 
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
        print(f"Found converted file {convertedfilename}")
        return False
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
    return True

def compute_weighted_average(convertedfilename):
    # will only compute the averaged SF for now
    # uproot.update is not yet implemented in version 2, will need to create a temporary file
    tmpfilename = "tmpfile.root"
    inputFile = uproot.open(convertedfilename)
    outputFile = uproot.recreate(tmpfilename)
    # first copy all existing histograms
    for key in inputFile.keys():
        outputFile[key] = inputFile[key]
    # close it for now but then open with pyROOT (having trouble to save histograms with uproot)
    outputFile.close()
    outputFile = ROOT.TFile.Open(tmpfilename,"update")
    # compute the average
    lepton_tags = ['w','z','g'] # to calculate average from these
    for year in [2017,2018]:
        for wp in ['loose', 'tight']:
            for prefix in ['mistag_SF']:
                h_centers = {x : inputFile[f'{prefix}_{x}_{wp}_{year}'].values             for x in lepton_tags}
                h_stat    = {x : np.sqrt(inputFile[f'{prefix}_{x}_{wp}_{year}'].variances) for x in lepton_tags}
                h_systUp  = {x : inputFile[f'{prefix}_{x}_{wp}_{year}_sysUp'].values       for x in lepton_tags}
                h_systDn  = {x : inputFile[f'{prefix}_{x}_{wp}_{year}_sysDn'].values       for x in lepton_tags}
                h_syst    = {x : 0.5*np.abs(h_systUp[x]-h_systDn[x])                      for x in lepton_tags}
                h_weights = {x : 1/(h_stat[x]*h_stat[x]+h_syst[x]*h_syst[x])              for x in lepton_tags}
                ave_values = []
                ave_variances = []
                edges = inputFile[f'{prefix}_{lepton_tags[0]}_{wp}_{year}'].edges
                for ibin in range(len(edges)-1):
                    ave_value = sum(h_centers[x][ibin] * h_weights[x][ibin] for x in lepton_tags) / sum(h_weights[x][ibin] for x in lepton_tags)
                    ave_variance = 1. / sum(h_weights[x][ibin] for x in lepton_tags)
                    ave_values.append(ave_value)
                    ave_variances.append(ave_variance)
                #df = inputFile[f'{prefix}_{lepton_tags[0]}_{wp}_{year}'].pandas()
                #df = df[1:-1] # remove overflow bins
                #df.index = pandas.IntervalIndex.from_breaks(edges)
                #df["count"] = ave_values
                #df["variance"] = ave_variances
                #outputFile[f'{prefix}_ave_{wp}_{year}'] = df
                htmp = ROOT.TH1F(f'{prefix}_ave_{wp}_{year}','',len(edges)-1, np.array(edges))
                for ibin in range(len(edges)-1):
                    htmp.SetBinContent(ibin+1, ave_values[ibin])
                    htmp.SetBinError(ibin+1, np.sqrt(ave_variances[ibin]))
                htmp.Write("", ROOT.TFile.kWriteDelete)
    # now replace the original file with the tmp file
    inputFile.close()
    outputFile.Close()
    os.rename(tmpfilename, convertedfilename)

def main():
    if convert_tefficiencies(inputfilename, convertedfilename):
        compute_weighted_average(convertedfilename)
    inputfile = uproot.open(convertedfilename)

    # first plot the rate plots
    fig, ax = plt.subplots(1,1, sharex=True)
    for year in [2017,2018]:
        for wp in ['loose', 'tight']:
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
                        if lepton_tag == "ave":
                            p = ax.errorbar(centers, hist.values, xerr=0.5*widths, yerr=np.sqrt(hist.variances), capsize=5, markeredgewidth=1, label=f"{lepton_tag} with all uncertainties", color=colors[lepton_tag])
                        else:
                            p = ax.errorbar(centers, hist.values, xerr=0.5*widths, yerr=np.sqrt(hist.variances), capsize=5, markeredgewidth=1, label=f"{lepton_tag} with stat. uncertainties", color=colors[lepton_tag])
                        if ("rate_data" in prefix) or ("SF" in prefix):
                            try:
                                hist_sysup = inputfile[hist_name+"_sysUp"]
                                hist_sysdn = inputfile[hist_name+"_sysDn"]
                                ax.bar(centers ,hist_sysup.values-hist_sysdn.values, bottom=hist_sysdn.values, width=widths, label=f"{lepton_tag} with syst. variation on top/diboson", color=p[0].get_color(), alpha=0.5)
                                # for print
                                if "SF" in prefix:
                                    msg = f"{wp} {prefix} {lepton_tag} {year}: , "
                                    for i in range(len(newbin.edges())-1):
                                        msg += f"{hist.values[i]} , {np.sqrt(hist.variances[i])} , {hist_sysup.values[i]-hist.values[i]} , "
                                    #print(msg)
                            except KeyError:
                                pass
                    ax.set_xlim(100,1100)
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
                plot_from_tags(['wz'])
                fig.savefig(f'{outdir}/{prefix}_{wp}_wz_{year}.png')
                plot_from_tags(['g'])
                fig.savefig(f'{outdir}/{prefix}_{wp}_g_{year}.png')
                if prefix == 'mistag_SF':
                    plot_from_tags(['g','w','z','ave'])
                    fig.savefig(f'{outdir}/{prefix}_{wp}_avegwz_{year}.png')
                    plot_from_tags(['ave'])
                    fig.savefig(f'{outdir}/{prefix}_{wp}_ave_{year}.png')
    print(f'Saved plots in {outdir}')



if __name__ == "__main__":
    main()
