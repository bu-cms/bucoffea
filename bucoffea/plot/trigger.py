#!/usr/bin/env python

import os
from bucoffea.plot.util import acc_from_dir, merge_datasets
from bucoffea.plot.style import markers
from coffea import hist
from coffea.hist.plot import clopper_pearson_interval
from matplotlib import pyplot as plt
import numpy as np
from tabulate import tabulate
pjoin = os.path.join

def trgname(year, tag):
    if year==2018 :
        if tag=='120pfht':
            return 'HLT_PFMETNoMu120_PFMHTNoMu120_IDTight(_PFHT60)'
        elif tag=='120only':
            return 'HLT_PFMETNoMu120_PFMHTNoMu120_IDTight'
    if year==2017:
        if tag=='120pfht':
            return 'HLT_PFMETNoMu120_PFMHTNoMu120_IDTight(_PFHT60)'
        elif tag=='120only':
            return 'HLT_PFMETNoMu120_PFMHTNoMu120_IDTight'


xmax = 1e3

def content_table(hnum, hden):
    table = []

    for x,ynum, yden in zip(hnum.axis('recoil').identifiers(),hnum.values()[()],hden.values()[()]):
        eff =  ynum/ yden if yden != 0 else 0
        unc = clopper_pearson_interval(ynum, yden, 0.68)
        line = [x, ynum, yden,eff,unc]
        table.append(line)
    return tabulate(table, headers=['Recoil', 'Numerator', 'Denominator',"Efficiency", "Uncertainty"])

def plot_recoil(acc, region_tag="1m", dataset='SingleMuon', year=2018, lumi=59.7, tag="test", distribution="recoil"):
    h = merge_datasets(acc[distribution])
    newbin = hist.Bin(distribution,"{distribution} (GeV)",np.array(list(range(0,400,20)) + list(range(400,1100,100))))
    h = h.rebin(h.axis(distribution), newbin)
    ds = f'{dataset}_{year}'
    h = h.project(h.axis('dataset'), ds)
    
    hnum = h.project(h.axis('region'),f'tr_{region_tag}_num')
    hden = h.project(h.axis('region'),f'tr_{region_tag}_den')
    
    
    # Recoil plot
    fig, ax,_ = hist.plot1d(hnum)
    hist.plot1d(hden, ax=ax, clear=False)

    outdir = f"./output/{tag}"
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    fig.savefig(pjoin(outdir, f'{distribution}_{region_tag}_{dataset}_{year}.pdf'))
    with open(pjoin(outdir,f'table_{region_tag}_{dataset}_{year}.txt'),"w") as f:
        f.write(content_table(hnum, hden))
    # Efficiency
    fig, ax,_ = hist.plotratio(hnum, hden,
                guide_opts={},
                unc='clopper-pearson',
                error_opts=markers('data')
                )
    ax.set_ylim(0,1.1)
    ax.set_xlim(0,xmax)
    ax.set_ylabel("Efficiency")

    plt.text(1., 1., r"%.1f fb$^{-1}$ (13 TeV)" % lumi,
                fontsize=16, 
                horizontalalignment='right', 
                verticalalignment='bottom', 
                transform=ax.transAxes
               )
    plt.text(0., 1., f'{region_tag}, {year}',
                fontsize=16, 
                horizontalalignment='left', 
                verticalalignment='bottom', 
                transform=ax.transAxes
               )
    plt.text(1., 0., f'{trgname(year, tag)}',
                fontsize=12, 
                horizontalalignment='right', 
                verticalalignment='bottom', 
                transform=ax.transAxes
               )

    plt.plot([0,xmax],[0.95,0.95],'r-')
    fig.savefig(pjoin(outdir, f'eff_{region_tag}_{dataset}_{year}.pdf'))



def main():
    # indir = "/home/albert/repos/bucoffea/bucoffea/plot/input/eff/test"

    tag = '120pfht'
    indir = f"/home/albert/repos/bucoffea/bucoffea/plot/input/eff/{tag}"
    acc = acc_from_dir(indir)

    # dataset="SingleMuon"
    # for region in ['1m', '2m']:
    #     plot_recoil(acc,region,dataset=dataset,year=2018, tag=tag)
    #     plot_recoil(acc,region,dataset=dataset,year=2017, lumi=41, tag=tag)
    
    
    for region in ['1e', '2e']:
        plot_recoil(acc,region,distribution="met",dataset="EGamma",year=2018, tag=tag)
        plot_recoil(acc,region,distribution="met",dataset="SingleElectron",year=2017, lumi=41, tag=tag)



if __name__ == "__main__":
    main()