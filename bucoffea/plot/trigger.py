#!/usr/bin/env python

import copy
import os
import re
from collections import defaultdict
from pprint import pprint

import numpy as np
import uproot
from coffea import hist
from coffea.hist.plot import clopper_pearson_interval
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
from tabulate import tabulate

from bucoffea.plot.style import markers, matplotlib_rc
from bucoffea.plot.util import (acc_from_dir, fig_ratio, lumi, merge_datasets,
                                merge_extensions, scale_xs_lumi)

from klepto.archives import dir_archive

matplotlib_rc()

pjoin = os.path.join

def lumi_by_region(region, year):
    if 'HLT' not in region:
        return lumi(year)
    else:
        if year == 2017:
            if 'HLT_PFHT590' in region: return 0.445
            if 'HLT_PFHT680' in region: return 0.786
            if 'HLT_PFHT780' in region: return 1.462
            if 'HLT_PFHT890' in region: return 2.802
            if 'HLT_PFHT1050' in region: return 41.527
        elif year==2018:
            if 'HLT_PFHT590' in region: return 0.467
            if 'HLT_PFHT680' in region: return 0.924
            if 'HLT_PFHT780' in region: return 1.838
            if 'HLT_PFHT890' in region: return 3.661
            if 'HLT_PFHT1050' in region: return 59.736

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
        line = [(x.lo + x.hi)/2, x.lo, x.hi, ynum, yden, eff,unc[0], unc[1]]
        table.append(line)
    return tabulate(table, headers=['Recoil', 'Numerator', 'Denominator',"Efficiency", "Eff-sigma","Eff+sigma"])

def plot_recoil(acc, region_tag="1m", dataset='SingleMuon', year=2018, tag="test", distribution="recoil",axis_name=None, noscale=False, jeteta_config=None, output_format='pdf'):
    # Select and prepare histogram
    h = copy.deepcopy(acc[distribution])
    h = merge_extensions(h, acc,reweight_pu=('nopu' in distribution), noscale=noscale)
    if not noscale:
        scale_xs_lumi(h)
    h = merge_datasets(h)

    # Rebinning
    axis_name = distribution if not axis_name else axis_name
    if 'photon' in distribution:
        newbin = hist.Bin(axis_name,f"{axis_name} (GeV)",np.array(list(range(0,250,10)) + list(range(250,400,50))+ list(range(400,1100,100))))
    elif distribution == 'mjj':
        newbin = hist.Bin(axis_name,r'$M_{jj}$ (GeV)',np.array(list(range(200,600,200)) + list(range(600,1500,300)) + [1500,2000,2750,3500]))
    else:
        newbin = hist.Bin(axis_name,f"{axis_name} (GeV)",np.array(list(range(0,400,20)) + list(range(400,1100,100))))
    h = h.rebin(h.axis(axis_name), newbin)
    ds = f'{dataset}_{year}'

    # Pick dataset and regions
    h = h.integrate(h.axis('dataset'), ds)
    if jeteta_config:
        hnum = h.integrate(h.axis('region'),f'tr_{region_tag}_num_{jeteta_config}')
        hden = h.integrate(h.axis('region'),f'tr_{region_tag}_den_{jeteta_config}')
    else:
        hnum = h.integrate(h.axis('region'),f'tr_{region_tag}_num')
        hden = h.integrate(h.axis('region'),f'tr_{region_tag}_den')

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

    outname = f'{region_tag}{"_noscale_" if noscale else "_"}{distribution}_{dataset}_{year}{"_"+jeteta_config if jeteta_config else ""}'

    fig.savefig(pjoin(outdir, f'{outname}.{output_format}'))
    with open(pjoin(outdir,f'table_{outname}.txt'),"w") as f:
        f.write(content_table(hnum, hden, axis_name) + "\n")
    plt.close(fig)

    # Efficiency plot
    fig, ax,_ = hist.plotratio(hnum, hden,
                guide_opts={},
                unc='clopper-pearson',
                error_opts=markers('data')
                )
    if distribution == 'recoil':
        ax.set_ylim(0.0,1.1)
    elif distribution == 'mjj':
        ax.set_ylim(0.8,1.1)
    ax.set_xlim(0,xmax)
    ax.set_ylabel("Efficiency")

    plt.text(1., 1., r"%.1f fb$^{-1}$ (13 TeV)" % lumi_by_region(region_tag,year),
                fontsize=16,
                horizontalalignment='right',
                verticalalignment='bottom',
                transform=ax.transAxes
                )
    plt.text(1., 0.95, f'{jeteta_config if jeteta_config else ""}',
                fontsize=12,
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

    if 'g_' in region_tag:
        plt.plot([215,215],[0.8,1.1],'r-')

    plt.plot([0,xmax],[0.95,0.95],'r-')
    fig.savefig(pjoin(outdir, f'eff_{outname}.pdf'))
    plt.close(fig)

def get_xy(file):
    data=np.loadtxt(file,skiprows=2)
    x = np.array(data[:,0])
    xedges = np.array(data[:,1:3])
    y = np.array(data[:,5])
    yerr = np.array(data[:,6:8])
    return x.T, xedges.T, y.T, np.abs(yerr.T-y.T)


colors = {
            '1m' : 'darkslategrey',
            '1m_hlt' : 'blue',
            '2m' : 'skyblue',
            '2m_hlt' : 'darkorange',
            '1e' : 'darkorange',
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
        regions = ['1m', '2m','2m_hlt']
        opts = markers('data')
        # opts['markersize'] = 1.
        # opts['fillstyle'] = 'none'
        emarker = opts.pop('emarker', '')

        fig, ax, rax = fig_ratio()

        x, y, yerr = {}, {}, {}
        for region in regions:
            if region.endswith('e'):
                file = f'output/{tag}/table_{region}_met_EGamma_{year}.txt'
            else:
                file = f'output/{tag}/table_{region}_recoil_SingleMuon_{year}.txt'
            x[region], _, y[region], yerr[region] = get_xy(file)
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
        plt.text(1., 1., r"%.1f fb$^{-1}$ (13 TeV)" % lumi_by_region(region, year),
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

    lo = np.hypot(numunc[0] / den, num*denunc[1]/den**2)
    high = np.hypot(numunc[1] / den, num*denunc[0]/den**2)
    return np.vstack((lo, high))

def sf_comparison_plot(tag):
    for year in [2017,2018]:
        regions = ['1m', '2m','1e']
        opts = markers('data')
        opts['markersize'] = 5
        # opts['fillstyle'] = 'none'
        emarker = opts.pop('emarker', '')

        fig, ax, rax = fig_ratio()

        x, y, yerr = {}, {}, {}
        for region in regions:
            if '1e' in region:
                fnum = f'output/{tag}/table_{region}_met_EGamma_{year}.txt'
                fden = f'output/{tag}/table_{region}_met_WJetsToLNu_HT_MLM_{year}.txt'
            elif '1m' in region:
                fnum = f'output/{tag}/table_{region}_recoil_SingleMuon_{year}.txt'
                fden = f'output/{tag}/table_{region}_recoil_WJetsToLNu_HT_MLM_{year}.txt'
            elif '2m' in region:
                fnum = f'output/{tag}/table_{region}_recoil_SingleMuon_{year}.txt'
                fden = f'output/{tag}/table_{region}_recoil_DYJetsToLL_M-50_HT_MLM_{year}.txt'


            xnum, _, ynum, yerrnum = get_xy(fnum)
            xden, _, yden, yerrden = get_xy(fden)
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
        rax.plot([250,250],[0.95,1.05],color='blue')
        ax.plot([250,250],[0.9,1.1],color='blue')
        plt.text(1., 1., r"%.1f fb$^{-1}$ (13 TeV)" % lumi_by_region(region, year),
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
    if 'gamma' in tag:
        regions = ['g_HLT_PFHT1050','g_HLT_PFHT590','g_HLT_PFHT680','g_HLT_PFHT780','g_HLT_PFHT890']
    else:
        regions = ['1m', '2m', '1e', '2m_hlt']
    opts = markers('data')
    # opts['markersize'] = 5
    # opts['fillstyle'] = 'none'
    emarker = opts.pop('emarker', '')
    outdir = f"./output/{tag}"
    outpath = pjoin(outdir,f'trig_sf.root')
    try:
        outfile = uproot.recreate(outpath)
    except OSError:
        outfile = uproot.update(outpath)

    for year in [2017,2018]:
        for region in regions:
            fig, ax, rax = fig_ratio()
            if '1e' in region:
                fnum = f'output/{tag}/table_{region}_met_EGamma_{year}.txt'
                fden = f'output/{tag}/table_{region}_met_WJetsToLNu_HT_MLM_{year}.txt'
                xlabel = "$p_{T}^{miss}$ (GeV)"
            elif '1m' in region:
                fnum = f'output/{tag}/table_{region}_recoil_SingleMuon_{year}.txt'
                fden = f'output/{tag}/table_{region}_recoil_WJetsToLNu_HT_MLM_{year}.txt'
                xlabel = "Recoil (GeV)"
            elif '2m' in region:
                fnum = f'output/{tag}/table_{region}_recoil_SingleMuon_{year}.txt'
                fden = f'output/{tag}/table_{region}_recoil_DYJetsToLL_M-50_HT_MLM_{year}.txt'
                xlabel = "Recoil (GeV)"
            elif 'g_' in region:
                fnum = f'output/{tag}/table_{region}_photon_pt0_JetHT_{year}.txt'
                fden = f'output/{tag}/table_{region}_photon_pt0_GJets_HT_MLM_{year}.txt'
                xlabel = "Photon $p_{T}$ (GeV)"

            if not os.path.exists(fnum):
                print(f"File not found {fnum}")
                continue
            if not os.path.exists(fden):
                print(f"File not found {fden}")
                continue

            xnum, xedgnum, ynum, yerrnum = get_xy(fnum)
            xden, xedgden, yden, yerrden = get_xy(fden)

            xsf = xnum
            ysf = ynum / yden
            ysferr = ratio_unc(ynum, yden, yerrnum, yerrden)

            opts['color'] = 'k'
            ax.errorbar(xnum, ynum, yerr=yerrnum,label=f'Data, {region} region', **opts)
            opts['color'] = 'r'
            ax.errorbar(xden, yden, yerr=yerrden,label=f'MC, {region} region', **opts)
            rax.plot([0,1000],[0.98,0.98],color='blue')
            rax.plot([0,1000],[0.99,0.99],color='blue',linestyle='--')

            if 'g_' in region:
                ax.plot([215,215],[0.9,1.1],color='blue')
                rax.plot([215,215],[0.95,1.05],color='blue')
            else:
                ax.plot([250,250],[0.9,1.1],color='blue')
                rax.plot([250,250],[0.95,1.05],color='blue')
            opts['color'] = 'k'
            rax.errorbar(xsf, ysf, ysferr, **opts)




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

            rax.set_xlabel(xlabel)
            plt.text(1., 1., r"$\approx$ %.1f fb$^{-1}$ (13 TeV)" % lumi_by_region(region, year),
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


            vals = np.array(sorted(list(set(list(xedgnum.flatten())))))
            ysf[np.isnan(ysf) | np.isinf(np.abs(ysf))] = 1
            outfile[f'{tag}_{region}_{year}'] = (ysf, vals)


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

def met_trigger_eff(distribution):
        if distribution == 'mjj':
            tag = '120pfht_mu_mjj'
        elif distribution == 'recoil':
            tag = '120pfht_mu_recoil'
            indir = '/afs/cern.ch/user/a/aakpinar/bucoffea/bucoffea/submission/2019-11-13_vbf_trigger_recoil'

        acc = dir_archive(
                          indir,
                          serialized=True,
                          compression=0,
                          memsize=1e3
                          )

        # Pre-load neccessary information 
        acc.load('recoil')      
        acc.load('sumw')      
        acc.load('sumw2')      

        for year in [2017, 2018]:
            for jeteta_config in ['two_central_jets', 'two_forward_jets', 'one_jet_forward_one_jet_central']:
                # Single muon CR
                region_tag='1m'
                for dataset in ['WJetsToLNu_HT_MLM', 'SingleMuon']:
                    plot_recoil(acc, region_tag=region_tag,
                                distribution=distribution,
                                axis_name=distribution,
                                dataset=dataset,
                                year=year,
                                tag=tag,
                                jeteta_config=jeteta_config,
                                output_format='pdf')        
                # Double muon CR
                region_tag='2m'
                for dataset in ['VDYJetsToLL_M-50_HT_MLM', 'SingleMuon']:
                    plot_recoil(acc, region_tag=region_tag,
                                distribution=distribution,
                                axis_name=distribution,
                                dataset=dataset,
                                year=year,
                                tag=tag,
                                jeteta_config=jeteta_config,
                                output_format='pdf')        

def met_triggers_ht():
        tag = '120pfht_hltmu'
        indir = f"/home/albert/repos/bucoffea/bucoffea/plot/input/16Jul19_incomplete_v7"
        acc = acc_from_dir(indir)

        # for noscale in [False]:
        #     distribution = 'recoil_noweight'  if noscale else 'recoil'
        #     for year in [2017, 2018]:
        #         region = '1m'
        #         for dataset in ["WJetsToLNu_HT_MLM", "SingleMuon"]:
        #             if dataset=='SingleMuon' and noscale:
        #                 continue
        #             plot_recoil(acc,region,distribution=distribution,axis_name='recoil',dataset=dataset,year=year, tag=tag, noscale=noscale)
        #         region = '2m'
        #         for dataset in ["DYJetsToLL_M-50_HT_MLM", "SingleMuon"]:
        #             if dataset=='SingleMuon' and noscale:
        #                 continue
        #             plot_recoil(acc,region,distribution=distribution,axis_name='recoil',dataset=dataset,year=year, tag=tag, noscale=noscale)
        #         # region = '1m_hlt'
        #         # for dataset in ["WJetsToLNu_HT_MLM", "SingleMuon"]:
        #         #     plot_recoil(acc,region,dataset=dataset,year=year, tag=tag)
        #         region = '2m_hlt'
        #         for dataset in ["DYJetsToLL_M-50_HT_MLM", "SingleMuon"]:
        #             if dataset=='SingleMuon' and noscale:
        #                 continue
        #             plot_recoil(acc,region,distribution=distribution,axis_name='recoil',dataset=dataset,year=year, tag=tag, noscale=noscale)
        #         region = '1e'
        #         for dataset in ["WJetsToLNu_HT_MLM", "EGamma"]:
        #             plot_recoil(acc,region,dataset=dataset,year=year, tag=tag, distribution='met')
        #         # region = '2e'
        #         # for dataset in ["DYJetsToLL_M-50_HT_MLM", "EGamma"]:
        #         #     plot_recoil(acc,region,dataset=dataset,year=year, tag=tag, distribution='met')

        # region_comparison_plot(tag)
        # sf_comparison_plot(tag)
        data_mc_comparison_plot(tag)

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
    indir = f"/home/albert/repos/bucoffea/bucoffea/plot/input/gamma_alltrig"
    acc = acc_from_dir(indir)

    # All regions
    regions = []
    for k in acc['photon_pt0'].axis('region').identifiers():
        k = str(k)
        if not k.startswith('tr_g'):
            continue
        if 'pt_trig_cut' in k:
            continue
        regions.append(k.replace("tr_","").replace("_num", "").replace("_den",""))

    for year in [2017, 2018]:
        # for distribution in ['photon_pt0', 'recoil']:
        for region in set(regions):
            for dataset in ["GJets_HT_MLM", "JetHT"]:
                if 'photon_pt_' in region:
                    distribution = 'recoil_noweight'
                    axis_name = 'recoil'
                else:
                    distribution = 'photon_pt0'
                    axis_name = 'pt'

                for noscale in [False]:
                    if ('JetHT' in dataset) and noscale:
                        continue
                    plot_recoil(
                                copy.deepcopy(acc),
                                region,
                                dataset=dataset,
                                year=year,
                                tag=tag,
                                distribution=distribution,
                                axis_name=axis_name,
                                noscale=noscale
                                )

                    # plot_recoil(acc,region,dataset=dataset,year=year, tag=tag, distribution='recoil',axis_name='recoil')
    data_mc_comparison_plot(tag)

def photon_sf_plot(tag):

    for year in [2017,2018]:
        opts = markers('data')
        # opts['markersize'] = 5
        # opts['fillstyle'] = 'none'
        opts.pop('emarker', '')

        fig = plt.gcf()
        ax = plt.gca()
        x, y, yerr = {}, {}, {}
        fnum = f'output/{tag}/table_g_HLT_PFHT1050_photon_pt0_JetHT_{year}.txt'
        fden = f'output/{tag}/table_g_HLT_PFHT590_photon_pt0_GJets_HT_MLM_{year}.txt'

        xnum, _, ynum, yerrnum = get_xy(fnum)
        xden, _, yden, yerrden = get_xy(fden)
        x = xnum
        y = ynum / yden
        dy = ratio_unc(ynum, yden, yerrnum, yerrden)

        ax.errorbar(x, y, yerr=dy,label=f'{year}', **opts)

        # for f in files: plot(f)
        outdir = f"./output/{tag}"
        # ax.set_ylim(0.9,1)
        ax.legend()
        ax.set_ylabel("Data / MC SF")
        ax.xaxis.set_major_locator(MultipleLocator(200))
        ax.xaxis.set_minor_locator(MultipleLocator(50))
        ax.yaxis.set_major_locator(MultipleLocator(0.05))
        ax.yaxis.set_minor_locator(MultipleLocator(0.01))
        ax.set_ylim(0.95,1.02)
        ax.grid(1)
        # rax.set_ylim(0.95,1.05)
        # rax.yaxis.set_major_locator(MultipleLocator(0.05))
        # rax.yaxis.set_minor_locator(MultipleLocator(0.01))
        # rax.grid(1)
        # rax.set_xlabel("Recoil or $p_{T}^{miss}$ (GeV)")
        # rax.set_ylabel(r"Ratio to single-$\mu$")
        # rax.plot([0,1000],[0.99,0.99],color='blue')
        # rax.plot([0,1000],[1.01,1.01],color='blue')
        # rax.plot([250,250],[0.95,0.95],color='blue')
        ax.plot([215,215],[0.9,1.1],color='red',zorder=-1)
        plt.text(1., 1., r"%.1f fb$^{-1}$ (13 TeV)" % lumi(year),
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
        fig.savefig(pjoin(outdir, f'photon_sf_{year}.pdf'))
        fig.clear()
        plt.close(fig)

def main():
    # indir = "/home/albert/repos/bucoffea/bucoffea/plot/input/eff/test"
    #met_triggers_ht()
    met_trigger_eff('recoil')
    # photon_triggers()
    # photon_sf_plot('gamma')
    # photon_triggers()


if __name__ == "__main__":
    main()
