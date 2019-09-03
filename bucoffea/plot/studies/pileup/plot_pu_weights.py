#!/usr/bin/env python

import ROOT as r
from bucoffea.helpers.paths import bucoffea_path
from bucoffea.plot.root_util import create_tdr_style, apply_style_to_axis, setup_canvas

def main():
    create_tdr_style()
    r.gROOT.SetBatch(r.kTRUE)
    fname=bucoffea_path('data/sf/pileup/pileup.root')

    f = r.TFile(fname)
    for year in [2017, 2018]:
        c,t1,t2=setup_canvas(want_ratio=True)
        hdata = f.Get(f'data{year}_nominal')
        hmc = f.Get(f'mc{year}')

        hmc.SetLineWidth(4)
        hmc.SetFillStyle(4001)
        hmc.SetFillColorAlpha(r.kAzure,0.25)
        
        hdata.SetMarkerStyle(20)
        hdata.SetMarkerSize(0.7)
        ratio = hdata.Clone("ratio")
        ratio.Divide(hmc)

        s = r.THStack()
        s.Add(hmc)
        s.Add(hdata,'PE')

        
        t1.cd()
        s.Draw('nostack')
        s.SetMaximum(.05)
        s.GetXaxis().SetRangeUser(0,80)

        leg = r.TLegend(0.6,0.7,0.9,0.9)
        leg.SetHeader(f'{year}')
        leg.AddEntry(hmc,f"MC, #mu = {hmc.GetMean():.1f}")
        leg.AddEntry(hdata,f"Data #mu = {hdata.GetMean():.1f}")
        leg.Draw()

        t2.cd()
        ratio.Draw('PE')

        apply_style_to_axis(s,is_ratio=False)
        apply_style_to_axis(ratio,is_ratio=True,ymin=0.,ymax=2.)
        ratio.SetMarkerStyle(20)

        s.GetXaxis().SetRangeUser(0,80)
        ratio.GetXaxis().SetRangeUser(0,80)
        ratio.GetXaxis().SetTitle("True number of PU events")
        c.SaveAs(f"pu_weights_{year}.pdf")

if __name__ == "__main__":
    main()

# for i in range(h11.)