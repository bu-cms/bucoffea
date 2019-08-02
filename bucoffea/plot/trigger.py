#!/usr/bin/env python

import copy
import os
import re
from collections import defaultdict
from pprint import pprint

import numpy as np
from coffea import hist
from coffea.hist.plot import clopper_pearson_interval
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
from tabulate import tabulate

from bucoffea.plot.style import markers
from bucoffea.plot.util import (acc_from_dir, fig_ratio, lumi, merge_datasets,
                                merge_extensions, scale_xs_lumi)

pjoin = os.path.join

def trgname(year, tag):
    if year==2018 :
        if '120pfht' in tag:
            return 'HLT_PFMETNoMu120_PFMHTNoMu120_IDTight(_PFHT60)'
        elif tag=='120only':
            return 'HLT_PFMETNoMu120_PFMHTNoMu120_IDTight'
        elif tag=='gamma':
            return 'HLT_Photon200'
    if year==2017:
        if '120pfht' in tag:
            return 'HLT_PFMETNoMu120_PFMHTNoMu120_IDTight(_PFHT60)'
        elif tag=='120only':
            return 'HLT_PFMETNoMu120_PFMHTNoMu120_IDTight'
        elif tag=='gamma':
            return 'HLT_Photon200'

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
        pprint(h.axis('region').identifiers())
        print(f'ERROR: {region_tag}, {dataset}, {year}')
        return
    hist.plot1d(hden, ax=ax, clear=False, binwnorm=True)
    plt.yscale('log')
    plt.gca().set_ylim(0.1,1e6)
    outdir = f"./output/{tag}"
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    fig.savefig(pjoin(outdir, f'{region_tag}_{distribution}_{dataset}_{year}.pdf'))
    with open(pjoin(outdir,f'table_{region_tag}_{distribution}_{dataset}_{year}.txt'),"w") as f:
        f.write(content_table(hnum, hden, axis_name) + "\n")
    plt.close(fig)

    # Efficiency plot
    fig, ax,_ = hist.plotratio(hnum, hden,
                guide_opts={},
                unc='clopper-pearson',
                error_opts=markers('data')
                )
    ax.set_ylim(0.8,1.1)
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
    fig.savefig(pjoin(outdir, f'eff_{region_tag}_{distribution}_{dataset}_{year}.pdf'))
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
            '1m' : 'o',
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
        opts['markersize'] = 1.
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
        regions = ['1m', '2m', '2m_hlt']
        opts = markers('data')
        opts['markersize'] = 5
        # opts['fillstyle'] = 'none'
        emarker = opts.pop('emarker', '')

        fig, ax, rax = fig_ratio()

        x, y, yerr = {}, {}, {}
        for region in regions:
            if '1e' in region:
                fnum = f'output/{tag}/table_{region}_EGamma_{year}.txt'
                fden = f'output/{tag}/table_{region}_WJetsToLNu_HT_MLM_{year}.txt'
            elif '1m' in region:
                fnum = f'output/{tag}/table_{region}_SingleMuon_{year}.txt'
                fden = f'output/{tag}/table_{region}_WJetsToLNu_HT_MLM_{year}.txt'
            elif '2m' in region:
                fnum = f'output/{tag}/table_{region}_SingleMuon_{year}.txt'
                fden = f'output/{tag}/table_{region}_DYJetsToLL_M-50_HT_MLM_{year}.txt'


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
        rax.plot([0,1000],[0.99,0.99],color='blue')
        rax.plot([0,1000],[1.01,1.01],color='blue')
        rax.plot([250,250],[0.95,0.95],color='blue')
        ax.plot([250,250],[0.9,1.1],color='blue')
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

def data_mc_comparison_plot(tag):
    regions = ['1m', '2m', '1e', '2m_hlt']
    opts = markers('data')
    # opts['markersize'] = 5
    # opts['fillstyle'] = 'none'
    emarker = opts.pop('emarker', '')

    for year in [2017,2018]:
        for region in regions:
            fig, ax, rax = fig_ratio()
            if '1e' in region:
                fnum = f'output/{tag}/table_{region}_EGamma_{year}.txt'
                fden = f'output/{tag}/table_{region}_WJetsToLNu_HT_MLM_{year}.txt'
            elif '1m' in region:
                fnum = f'output/{tag}/table_{region}_SingleMuon_{year}.txt'
                fden = f'output/{tag}/table_{region}_WJetsToLNu_HT_MLM_{year}.txt'
            elif '2m' in region:
                fnum = f'output/{tag}/table_{region}_SingleMuon_{year}.txt'
                fden = f'output/{tag}/table_{region}_DYJetsToLL_M-50_HT_MLM_{year}.txt'

            if not os.path.exists(fnum):
                print(f"File not found {fnum}")
                continue
            if not os.path.exists(fden):
                print(f"File not found {fden}")
                continue

            xnum, ynum, yerrnum = get_xy(fnum)
            xden, yden, yerrden = get_xy(fden)

            xsf = xnum
            ysf = ynum / yden
            ysferr = ratio_unc(ynum, yden, yerrnum, yerrden)

            opts['color'] = 'k'
            ax.errorbar(xnum, ynum, yerr=yerrnum,label=f'Data, {region} region', **opts)
            opts['color'] = 'r'
            ax.errorbar(xden, yden, yerr=yerrden,label=f'MC, {region} region', **opts)
            rax.plot([0,1000],[0.98,0.98],color='blue')
            rax.plot([0,1000],[0.99,0.99],color='blue',linestyle='--')
            ax.plot([250,250],[0.9,1.1],color='blue')
            rax.plot([250,250],[0.95,1.05],color='blue')
            opts['color'] = 'k'
            rax.errorbar(xsf, ysf, ysferr, **opts)



            outdir = f"./output/{tag}"
        # ax.set_ylim(0.9,1)
            ax.legend()
            ax.set_ylabel("Efficiency")
            rax.set_ylabel("Data / MC SF")
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
            fig.savefig(pjoin(outdir, f'data_mc_comparison_{region}_{year}.pdf'))
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

def met_triggers_ht():
        tag = 'gamma'
        indir = f"/home/albert/repos/bucoffea/bucoffea/plot/input/eff/{tag}/"
        acc = acc_from_dir(indir)

        for year in [2017, 2018]:
            region = '1m'
            for dataset in ["WJetsToLNu_HT_MLM", "SingleMuon"]:
                plot_recoil(acc,region,dataset=dataset,year=year, tag=tag)
            region = '2m'
            for dataset in ["DYJetsToLL_M-50_HT_MLM", "SingleMuon"]:
                plot_recoil(acc,region,dataset=dataset,year=year, tag=tag)
        #     region = '1m_hlt'
        #     for dataset in ["WJetsToLNu_HT_MLM", "SingleMuon"]:
        #         plot_recoil(acc,region,dataset=dataset,year=year, tag=tag)
        #     region = '2m_hlt'
        #     for dataset in ["DYJetsToLL_M-50_HT_MLM", "SingleMuon"]:
        #         plot_recoil(acc,region,dataset=dataset,year=year, tag=tag)
            region = '1e'
            for dataset in ["WJetsToLNu_HT_MLM", "EGamma"]:
                plot_recoil(acc,region,dataset=dataset,year=year, tag=tag, distribution='met')
        #     region = '2e'
        #     for dataset in ["DYJetsToLL_M-50_HT_MLM", "EGamma"]:
        #         plot_recoil(acc,region,dataset=dataset,year=year, tag=tag, distribution='met')

        region_comparison_plot(tag)
        sf_comparison_plot(tag)

def photon_triggers_merged():
    tag = 'gamma'
    indir = f"/home/albert/repos/bucoffea/bucoffea/plot/input/eff/{tag}/sel"
    acc = acc_from_dir(indir)

    # All regions
    regions = []
    for k in acc['photon_pt0'].axis('region').identifiers():
        k = str(k)
        if not k.startswith('tr_g'):
            continue
        regions.append(k)

    lumi_by_trig = {
        2017: {
        "HLT_PFHT1050" : 41.527192272,
        "HLT_PFHT590" : 0.444896852,
        "HLT_PFHT680" : 0.785538800,
        "HLT_PFHT780" : 1.461679046,
        "HLT_PFHT890" : 2.802018098,
        },
        2018: {
        "HLT_PFHT1050" : 59.735969368,
        "HLT_PFHT590" :  0.467058719,
        "HLT_PFHT680" :  0.924441113,
        "HLT_PFHT780" :  1.837625206,
        "HLT_PFHT890" :  3.661338636,
        }
    }
    overall_lumi = {
        2017 : 41.527192272,
        2018 : 59.735969368
    }
    for year in [2017, 2018]:
        acc_cp = copy.deepcopy(acc)

        # Region combination mapping
        # 1. Scale each region according to effective lumi
        # 2. Combine regions with different triggers
        mapping = defaultdict(set)
        lumi_mapping = {}
        for rname in regions:
            base = re.sub(r'_HLT_PFHT(\d+)','', rname)
            mapping[base].add(rname)

            for hlt, lumi in lumi_by_trig[year].items():
                if hlt in rname:
                    lumi_mapping[rname] = lumi / overall_lumi[year]

        pprint(lumi_mapping)


        pprint(mapping)
        region_ax = hist.Cat("region", "Selection region")
        year_regex = re.compile(f'.*{year}')




        for distribution in ['photon_pt0', 'recoil']:
            h = copy.deepcopy(acc[distribution])
            h = h[year_regex]
            h.scale(lumi_mapping, 'region')
            acc_cp[distribution] = h.group(region_ax, "region", {k:list(v) for k, v in mapping.items()})


            for region in mapping.keys():
                print(f"Region {region}")
                region = region.replace("tr_","").replace("_num", "").replace("_den","")
                for dataset in ["GJets_HT_MLM", "JetHT"]:
                    if 'photon_pt_' in region:
                        distribution = 'recoil'
                        axis_name = 'recoil'
                    else:
                        distribution = 'photon_pt0'
                        axis_name = 'pt'
                    plot_recoil(
                                acc_cp,
                                region,
                                dataset=dataset,
                                year=year,
                                tag=tag,
                                distribution=distribution,
                                axis_name=axis_name
                                )
                    # plot_recoil(acc,region,dataset=dataset,year=year, tag=tag, distribution='recoil',axis_name='recoil')

def photon_triggers():
    tag = 'gamma'
    indir = f"/home/albert/repos/bucoffea/bucoffea/plot/input/eff/{tag}/sel"
    acc = acc_from_dir(indir)

    # All regions
    regions = []
    for k in acc['photon_pt0'].axis('region').identifiers():
        k = str(k)
        if not k.startswith('tr_g'):
            continue
        regions.append(k.replace("tr_","").replace("_num", "").replace("_den",""))

    for year in [2017, 2018]:
        # for distribution in ['photon_pt0', 'recoil']:
        for region in set(regions):
            for dataset in ["GJets_HT_MLM", "JetHT"]:
                if 'photon_pt_' in region:
                    distribution = 'recoil'
                    axis_name = 'recoil'
                else:
                    distribution = 'photon_pt0'
                    axis_name = 'pt'
                acc_cp = copy.deepcopy(acc)
                plot_recoil(
                            acc_cp,
                            region,
                            dataset=dataset,
                            year=year,
                            tag=tag,
                            distribution=distribution,
                            axis_name=axis_name
                            )
                    # plot_recoil(acc,region,dataset=dataset,year=year, tag=tag, distribution='recoil',axis_name='recoil')
def main():
    # indir = "/home/albert/repos/bucoffea/bucoffea/plot/input/eff/test"
    met_triggers_ht()
    # data_mc_comparison_plot('gamma')
    # photon_triggers()


if __name__ == "__main__":
    main()
