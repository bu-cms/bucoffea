import ROOT
import time
from mistagSF import outdir

inputfilename = f"{outdir}/wtag_mistag_SF.root"
inputfile = ROOT.TFile.Open(inputfilename,'read')

colors={
        '1m':6,
        '2m':7,
        '1e':8,
        '2e':9,
        'g':11,
        'combined':46,
        }
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(ROOT.kFALSE)

canv = ROOT.TCanvas('canv','canv',800,800)
for year in [2017,2018]:
    for wp in ['loose', 'tight']:
        # plot the mistag rates 
        for prefix in ['mistag_rate_data','mistag_rate_mc']:
            canv.Clear()
            canvEmpty=True
            for lepton_flag in ['g','combined']:
                if lepton_flag=='combined':
                    htmp = inputfile.Get(f'{prefix}_{wp}_{year}')
                    htmp.SetTitle('combined except g') # for the legend builder
                else:
                    htmp = inputfile.Get(f'{prefix}_{lepton_flag}_{wp}_{year}')
                    htmp.SetLineStyle(2)
                    htmp.SetTitle(lepton_flag) # for the legend builder
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
                        gtmp.GetYaxis().SetRangeUser(0,0.14)
                    else:
                        gtmp.GetYaxis().SetRangeUser(0,0.8)
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
            for lepton_flag in ['g','combined']:
                if lepton_flag=='combined':
                    htmp = inputfile.Get(f'Wmistag_{year}_{wp}_ak8_pt')
                    htmp.SetTitle("combined except g")
                else:
                    htmp = inputfile.Get(f'{prefix}_{lepton_flag}_{wp}_{year}')
                    htmp.SetLineStyle(2)
                    htmp.SetTitle(lepton_flag) # for the legend builder
                htmp.SetLineColor(colors[lepton_flag])
                htmp.SetLineWidth(3)
                if canvEmpty:
                    if 'tight' in wp:
                        htmp.GetYaxis().SetRangeUser(0.5,1.7)
                    else:
                        htmp.GetYaxis().SetRangeUser(1,1.3)
                    htmp.GetYaxis().SetTitle('mistagging SF')
                    htmp.GetXaxis().SetTitle('AK8 Jet p_{T}')
                    htmp.Draw('e1')
                    canvEmpty=False
                else:
                    htmp.Draw('e1 same')
            canv.BuildLegend()
            canv.Update()
            canv.Print(f'{outdir}/{prefix}_{wp}_{year}.png')
