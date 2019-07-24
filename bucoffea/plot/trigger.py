#!/usr/bin/env python

import os
from bucoffea.plot.util import acc_from_dir, merge_datasets, merge_extensions, scale_xs_lumi, fig_ratio, lumi
from bucoffea.plot.style import markers
from coffea import hist
from coffea.hist.plot import clopper_pearson_interval
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np
from tabulate import tabulate
from pprint import pprint
import copy
pjoin = os.path.join

def trgname(year, tag):
    if year==2018 :
        if '120pfht' in tag:
            return 'HLT_PFMETNoMu120_PFMHTNoMu120_IDTight(_PFHT60)'
        elif tag=='120only':
            return 'HLT_PFMETNoMu120_PFMHTNoMu120_IDTight'
    if year==2017:
        if '120pfht' in tag:
            return 'HLT_PFMETNoMu120_PFMHTNoMu120_IDTight(_PFHT60)'
        elif tag=='120only':
            return 'HLT_PFMETNoMu120_PFMHTNoMu120_IDTight'


xmax = 1e3

def content_table(hnum, hden, axis_name):
    table = []
    for x,ynum, yden in zip(hnum.axis(axis_name).identifiers(),hnum.values()[()],hden.values()[()]):
        eff =  ynum/ yden if yden != 0 else 0
        unc = clopper_pearson_interval(ynum, yden, 0.68)
        line = [(x.lo + x.hi)/2, ynum, yden, eff,unc[0], unc[1]]
        table.append(line)
    return tabulate(table, headers=['Recoil', 'Numerator', 'Denominator',"Efficiency", "Eff-sigma","Eff+sigma"])

def plot_recoil(acc, region_tag="1m", dataset='SingleMuon', year=2018, tag="test", distribution="recoil",axis_name=None):
    # Select and prepare histogram
    h = copy.deepcopy(acc[distribution])
    h = merge_extensions(h, acc)
    scale_xs_lumi(h)
    h = merge_datasets(h)

    # Rebinning
    axis_name = distribution if not axis_name else axis_name
    newbin = hist.Bin(axis_name,f"{axis_name} (GeV)",np.array(list(range(0,400,20)) + list(range(400,1100,100))))
    h = h.rebin(h.axis(axis_name), newbin)
    ds = f'{dataset}_{year}'

    # Pick dataset and regions
    h = h.project(h.axis('dataset'), ds)
    hnum = h.project(h.axis('region'),f'tr_{region_tag}_num')
    hden = h.project(h.axis('region'),f'tr_{region_tag}_den')

    # Recoil plot
    try:
        fig, ax,_ = hist.plot1d(hnum, binwnorm=True)
    except KeyError:
        print(f'ERROR: {region_tag}, {dataset}, {year}')
        return
    hist.plot1d(hden, ax=ax, clear=False, binwnorm=True)
    plt.yscale('log')
    plt.gca().set_ylim(0.1,1e6)
    outdir = f"./output/{tag}"
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    fig.savefig(pjoin(outdir, f'{distribution}_{region_tag}_{dataset}_{year}.pdf'))
    with open(pjoin(outdir,f'table_{region_tag}_{dataset}_{year}.txt'),"w") as f:
        f.write(content_table(hnum, hden, axis_name) + "\n")
    plt.close(fig)

    # Efficiency plot
    fig, ax,_ = hist.plotratio(hnum, hden,
                guide_opts={},
                unc='clopper-pearson',
                error_opts=markers('data')
                )
    ax.set_ylim(0,1.1)
    ax.set_xlim(0,xmax)
    ax.set_ylabel("Efficiency")

    plt.text(1., 1., r"$\approx$ %.1f fb$^{-1}$ (13 TeV)" % lumi(year),
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
                fontsize=10,
                horizontalalignment='right',
                verticalalignment='bottom',
                transform=ax.transAxes
               )

    plt.plot([0,xmax],[0.95,0.95],'r-')
    fig.savefig(pjoin(outdir, f'eff_{region_tag}_{dataset}_{year}.pdf'))
    plt.close(fig)

def get_xy(file):
    data=np.loadtxt(file,skiprows=2)
    x = np.array(data[:,0])
    y = np.array(data[:,3])
    yerr = np.array(data[:,4:6])
    return x.T, y.T, np.abs(yerr.T-y.T)


colors = {
            '1m' : 'darkslategrey',
            '1m_hlt' : 'blue',
            '2m' : 'darkorange',
            '2m_hlt' : 'red',
            '1e' : 'skyblue',
            '2e' : 'magenta'
        }
region_marker = {
            '1m' : '.',
            '1m_hlt' : 'o',
            '2m' : 's',
            '2m_hlt' : '^',
            '1e' : 'v',
            '2e' : '>'
        }


def region_comparison_plot(tag):
    for year in [2017,2018]:
        regions = ['1m', '2m', '1e','1m_hlt','2m_hlt']
        opts = markers('data')
        opts['markersize'] = 5.
        opts['fillstyle'] = 'none'
        emarker = opts.pop('emarker', '')

        fig, ax, rax = fig_ratio()

        x, y, yerr = {}, {}, {}
        for region in regions:
            if region.endswith('e'):
                file = f'output/{tag}/table_{region}_EGamma_{year}.txt'
            else:
                file = f'output/{tag}/table_{region}_SingleMuon_{year}.txt'
            x[region], y[region], yerr[region] = get_xy(file)
            opts['color'] = colors[region]
            ax.errorbar(x[region], y[region], yerr=yerr[region],label=f'{region} region', **opts)

        # opts.pop('elinewidth')
            if region=='1m':
                continue

            rax.errorbar(x['1m'], y[region]/y['1m'], yerr[region]/y['1m'], **opts)

        # for f in files: plot(f)
        outdir = f"./output/{tag}"
        # ax.set_ylim(0.9,1)
        ax.legend()
        ax.set_ylabel("Efficiency")
        ax.xaxis.set_major_locator(MultipleLocator(200))
        ax.xaxis.set_minor_locator(MultipleLocator(50))
        ax.yaxis.set_major_locator(MultipleLocator(0.05))
        ax.yaxis.set_minor_locator(MultipleLocator(0.01))
        ax.set_ylim(0.9,1.02)
        ax.grid(1)
        rax.set_ylim(0.9,1.1)
        rax.grid(1)
        rax.set_xlabel("Recoil or $p_{T}^{miss}$ (GeV)")
        rax.set_ylabel(r"Ratio to single-$\mu$")
        plt.text(1., 1., r"$\approx$ %.1f fb$^{-1}$ (13 TeV)" % lumi(year),
                fontsize=16,
                horizontalalignment='right',
                verticalalignment='bottom',
                transform=ax.transAxes
               )
        plt.text(0., 1., f'{year}',
                fontsize=16,
                horizontalalignment='left',
                verticalalignment='bottom',
                transform=ax.transAxes
               )
        fig.savefig(pjoin(outdir, f'region_comparison_data_{year}.pdf'))
        fig.clear()
        plt.close(fig)

def ratio_unc(num, den, numunc, denunc):
    return np.hypot(numunc / den, num*denunc/den**2)


def sf_comparison_plot(tag):
    for year in [2017,2018]:
        regions = ['1m', '2m', '1m_hlt', '2m_hlt']
        opts = markers('data')
        opts['markersize'] = 5
        opts['fillstyle'] = 'none'
        emarker = opts.pop('emarker', '')

        fig, ax, rax = fig_ratio()

        x, y, yerr = {}, {}, {}
        for region in regions:
            if '1e' in region:
                fnum = f'output/{tag}/table_{region}_EGamma_{year}.txt'
                fden = f'output/{tag}/table_{region}_WJetsToLNu-MLM_{year}.txt'
            elif '1m' in region:
                fnum = f'output/{tag}/table_{region}_SingleMuon_{year}.txt'
                fden = f'output/{tag}/table_{region}_WJetsToLNu-MLM_{year}.txt'
            elif '2m' in region:
                fnum = f'output/{tag}/table_{region}_SingleMuon_{year}.txt'
                fden = f'output/{tag}/table_{region}_DYNJetsToLL_M-50-MLM_{year}.txt'


            xnum, ynum, yerrnum = get_xy(fnum)
            xden, yden, yerrden = get_xy(fden)
            x[region] = xnum
            y[region] = ynum / yden
            yerr[region] = ratio_unc(ynum, yden, yerrnum, yerrden)
            opts['color'] = colors[region]
            opts['marker'] = region_marker[region]
            ax.errorbar(x[region], y[region], yerr=yerr[region],label=f'{region} region', **opts)

        # opts.pop('elinewidth')
            if region=='1m':
                continue

            rax.errorbar(x['1m'], y[region]/y['1m'], ratio_unc(y[region],y['1m'],yerr[region],yerr['1m']), **opts)

        # for f in files: plot(f)
        outdir = f"./output/{tag}"
        # ax.set_ylim(0.9,1)
        ax.legend()
        ax.set_ylabel("Data / MC SF")
        ax.xaxis.set_major_locator(MultipleLocator(200))
        ax.xaxis.set_minor_locator(MultipleLocator(50))
        ax.yaxis.set_major_locator(MultipleLocator(0.05))
        ax.yaxis.set_minor_locator(MultipleLocator(0.01))
        ax.set_ylim(0.9,1.1)
        ax.grid(1)
        rax.set_ylim(0.95,1.05)
        rax.yaxis.set_major_locator(MultipleLocator(0.05))
        rax.yaxis.set_minor_locator(MultipleLocator(0.01))
        rax.grid(1)
        rax.set_xlabel("Recoil or $p_{T}^{miss}$ (GeV)")
        rax.set_ylabel(r"Ratio to single-$\mu$")
        plt.text(1., 1., r"$\approx$ %.1f fb$^{-1}$ (13 TeV)" % lumi(year),
                fontsize=16,
                horizontalalignment='right',
                verticalalignment='bottom',
                transform=ax.transAxes
               )
        plt.text(0., 1., f'{year}',
                fontsize=16,
                horizontalalignment='left',
                verticalalignment='bottom',
                transform=ax.transAxes
               )
        fig.savefig(pjoin(outdir, f'sf_comparison_{year}.pdf'))
        fig.clear()
        plt.close(fig)

def met_triggers():
        tag = '120pfht_hltmu'
        indir = f"/home/albert/repos/bucoffea/bucoffea/plot/input/eff/{tag}"
        acc = acc_from_dir(indir)

        for year in [2017, 2018]:
            region = '1m'
            for dataset in ["WJetsToLNu-MLM", "SingleMuon"]:
                plot_recoil(acc,region,dataset=dataset,year=year, tag=tag)
            region = '2m'
            for dataset in ["DYNJetsToLL_M-50-MLM", "SingleMuon"]:
                plot_recoil(acc,region,dataset=dataset,year=year, tag=tag)
            region = '1m_hlt'
            for dataset in ["WJetsToLNu-MLM", "SingleMuon"]:
                plot_recoil(acc,region,dataset=dataset,year=year, tag=tag)
            region = '2m_hlt'
            for dataset in ["DYNJetsToLL_M-50-MLM", "SingleMuon"]:
                plot_recoil(acc,region,dataset=dataset,year=year, tag=tag)
            region = '1e'
            for dataset in ["WJetsToLNu-MLM", "EGamma"]:
                plot_recoil(acc,region,dataset=dataset,year=year, tag=tag, distribution='met')
            region = '2e'
            for dataset in ["DYNJetsToLL_M-50-MLM", "EGamma"]:
                plot_recoil(acc,region,dataset=dataset,year=year, tag=tag, distribution='met')

        region_comparison_plot(tag)
        sf_comparison_plot(tag)

def photon_triggers():
    tag = 'gamma'
    indir = f"/home/albert/repos/bucoffea/bucoffea/plot/input/eff/{tag}/sel"
    acc = acc_from_dir(indir)

    for year in [2017, 2018]:
            region = 'g'
            for dataset in ["GJets_HT_MLM", "JetHT"]:
                plot_recoil(acc,region,dataset=dataset,year=year, tag=tag, distribution='photonpt0',axis_name='pt')
def main():
    # indir = "/home/albert/repos/bucoffea/bucoffea/plot/input/eff/test"
    photon_triggers()


if __name__ == "__main__":
    main()