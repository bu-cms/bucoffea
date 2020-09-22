#!/usr/bin/env python

import copy
import os
import re
from collections import defaultdict
from pprint import pprint

import matplotlib.ticker
import numpy as np
from coffea import hist
from coffea.util import load
from matplotlib import pyplot as plt

from bucoffea.execute.dataset_definitions import short_name
from bucoffea.helpers.dataset import is_data
from bucoffea.helpers.paths import bucoffea_path
from bucoffea.plot.util import (acc_from_dir, lumi, merge_datasets,
                                merge_extensions, scale_xs_lumi)
from bucoffea.plot import style

pjoin = os.path.join
#suppress true_divide warnings
np.seterr(divide='ignore', invalid='ignore')
import matplotlib.pylab as pylab
params = {'legend.fontsize': 'small',
        #   'figure.figsize': (15, 5),
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
pylab.rcParams.update(params)
colors = {
    'WN*J.*' : '#feb24c',
    '.*DY.*' : '#ffffcc',
    'EWKW.*' : '#c6dbef',
    'EWKZ.*ZToLL.*' : '#d5bae2',
    'EWKZ.*ZToNuNu.*' : '#c4cae2',
    'EWK_V.*' : '#deebf7',
    '.*Diboson.*' : '#4292c6',
    '.*WZ.*' : '#a6bddb',
    '.*ZZ.*' : '#3690c0',
    '.*WW.*' : '#034e7b',
    'Top.*' : '#6a51a3',
    '.*TT.*' : '#6a51a3',
    '.*ST.*' : '#9e9ac8',
    '.*QCD.*' : '#08306b',
    'GJets_HT.*' : '#fc4e2a',
    'GJets_DR-0p4.*' : '#fc4e2a',
    'GJets_SM.*' : '#a76b51',


    'VQQGamma.*' : '#51b84f',
    'ZQQGamma.*' : '#3a8739',
    'WQQGamma.*' : '#51b84f',
    'ZJetsToNuNu.*' : '#31a354',
    'ZNuNuGJets_.*' : '#0050ec'
}
legend_labels = {
    'VBF': {
        'GJets_DR-0p4.*' : "QCD $\\gamma$+jets",
        'GJets_SM.*' : "EWK $\\gamma$+jets",
        'DY.*' : "QCD Z$\\rightarrow\\ell\\ell$",
        'EWKZ.*ZToLL.*' : "EWK Z$\\rightarrow\\ell\\ell$",
        'WN*J.*LNu.*' : "QCD W$\\rightarrow\\ell\\nu$",
        'EWKW.*LNu.*' : "EWK W$\\rightarrow\\ell\\nu$",
        'ZJetsToNuNu.*.*' : "QCD Z$\\rightarrow\\nu\\nu$",
        'EWKZ.*ZToNuNu.*' : "EWK Z$\\rightarrow\\nu\\nu$"
    },
    'Monojet/Mono-V': {
        'GJets_DR-0p4.*' : "$\\gamma$+jets",
        'DY.*' : "Z$\\rightarrow\\ell\\ell$",
        'WN*J.*LNu.*' : "W$\\rightarrow\\ell\\nu$",
        'ZJetsToNuNu.*.*' : "Z$\\rightarrow\\nu\\nu$"
    },
    'Common': {
        'QCD.*' : "QCD",
        'Top.*' : "Top quark",
        'Diboson.*' : "WW/WZ/ZZ",
        'ZQQGamma.*' : 'Z(qq)$\gamma$',
        'WQQGamma.*' : 'W(qq)$\gamma$',
        'VQQGamma.*' : 'V(qq)$\gamma$',
        'WW.*' : 'WW',
        'ZZ.*' : 'ZZ',
        'WZ.*' : 'WZ',
        'MET|Single(Electron|Photon|Muon)|EGamma.*' : "Data"
    }
}
class Style():
    def __init__(self):
        self.region_names = {
            'cr_1m_j' : 'Single-$\mu$ CR, monojet',
            'cr_2m_j' : 'Di-$\mu$ CR, monojet',
            'cr_1e_j' : 'Single-e CR, monojet',
            'cr_2e_j' : 'Di-e CR, monojet',
            'cr_g_j' : 'Single-Photon CR, monojet',
            'sr_j' : 'Signal region, monojet',
            'cr_1m_v' : 'Single-$\mu$ CR, mono-v',
            'cr_2m_v' : 'Di-$\mu$ CR, mono-v',
            'cr_1e_v' : 'Single-e CR, mono-v',
            'cr_2e_v' : 'Di-e CR, mono-v',
            'cr_g_v' : 'Single-Photon CR, mono-v',
            'sr_vbf' : 'Signal region, vbfhinv',
            'cr_1m_vbf' : 'Single-$\mu$ CR, vbfhinv',
            'cr_2m_vbf' : 'Di-$\mu$ CR, vbfhinv',
            'cr_1e_vbf' : 'Single-e CR, vbfhinv',
            'cr_2e_vbf' : 'Di-e CR, vbfhinv',
            'cr_g_vbf' : 'Single-Photon CR, vbfhinv'
        }
        recoil_bins_2016 = [ 250,  280,  310,  340,  370,  400,  430,  470,  510, 550,  590,  640,  690,  740,  790,  840,  900,  960, 1020, 1090, 1160, 1250, 1400]
        recoil_monov_bins_2016 = [250,300,350,400,500,600,750,1000]
        self.binnings = {
            'default': {
                    'dimuon_mass' : hist.Bin('dilepton_mass',r'M($\mu^{+}\mu^{-}$)',30,60,120),
                    'dielectron_mass' : hist.Bin('dilepton_mass',r'M($e^{+}e^{-}$)',30,60,120),
                    # 'recoil' : hist.Bin('recoil','recoil',list(range(250,300,50)) + list(range(300,500,50)) + list(range(500,1000,100)) + list(range(1000,2000,200))),
                    'recoil' : hist.Bin('recoil','Recoil (GeV)', recoil_bins_2016),
                    'recoil_hardbveto' : hist.Bin('recoil','Recoil (GeV)', recoil_bins_2016),
                    'recoil_nopog' : hist.Bin('recoil','Recoil (GeV)', recoil_bins_2016),
                    'recoil_nopu' : hist.Bin('recoil','Recoil (GeV)', recoil_bins_2016),
                    'recoil_nopref' : hist.Bin('recoil','Recoil (GeV)', recoil_bins_2016),
                    'recoil_notheory' : hist.Bin('recoil','Recoil (GeV)', recoil_bins_2016),
                    'recoil_recosfrecoil' : hist.Bin('recoil','Recoil (GeV)', recoil_bins_2016),
                    'recoil_nodibosonnnlo' : hist.Bin('recoil','Recoil (GeV)', recoil_bins_2016),
                    'met' : hist.Bin('met',r'$p_{T}^{miss}$ (GeV)',list(range(0,500,50)) + list(range(500,1000,100)) + list(range(1000,2000,250))),
                    'calomet' : hist.Bin('met',r'$p_{T}^{miss}$ (GeV)',list(range(0,500,50)) + list(range(500,1000,100)) + list(range(1000,2000,250))),
                    'puppimet' : hist.Bin('met',r'$p_{T}^{miss}$ (GeV)',list(range(0,500,50)) + list(range(500,1000,100)) + list(range(1000,2000,250))),
                    'tkmet' : hist.Bin('met',r'$p_{T}^{miss}$ (GeV)',list(range(0,500,50)) + list(range(500,1000,100)) + list(range(1000,2000,250))),
                    'ak4_pt0' : hist.Bin('jetpt',r'Leading AK4 jet $p_{T}$ (GeV)',list(range(80,600,20)) + list(range(600,1000,20)) ),
                    'ak4_pt1' : hist.Bin('jetpt',r'Trailing AK4 jet $p_{T}$ (GeV)',list(range(40,600,20)) + list(range(600,1000,20)) ),
                    'ak4_pt' : hist.Bin('jetpt',r'All AK4 jet $p_{T}$ (GeV)',list(range(100,600,20)) + list(range(600,1000,20)) ),
                    'ak4_ptraw0' : hist.Bin('jetpt',r'Leading AK4 raw jet $p_{T}$ (GeV)',list(range(100,600,20)) + list(range(600,1000,20)) ),
                    'ak4_pt0_eta0' : hist.Bin('jetpt',r'AK4 jet $p_{T}$ (GeV)',list(range(100,600,20)) + list(range(600,1000,20)) ),
                    'ak8_pt0' : hist.Bin('jetpt',r'Leading AK8 jet $p_{T}$ (GeV)',list(range(80,600,20)) + list(range(600,1000,40)) ),
                    'ak8_pt1' : hist.Bin('jetpt',r'Trailing AK8 jet $p_{T}$ (GeV)',list(range(40,600,20)) + list(range(600,1000,40)) ),
                    'ak8_pt' : hist.Bin('jetpt',r'All AK8 jet $p_{T}$ (GeV)',list(range(100,600,20)) + list(range(600,1000,40)) ),
                    'ak8_ptraw0' : hist.Bin('jetpt',r'Leading AK8 raw jet $p_{T}$ (GeV)',list(range(100,600,20)) + list(range(600,1000,40)) ),
                    'ak8_pt0_eta0' : hist.Bin('jetpt','jetpt',list(range(100,600,20)) + list(range(600,1000,40)) ),
                    'photon_pt0' : hist.Bin('pt',r'Photon $p_{T}$ (GeV)',list(range(200,600,20)) + list(range(600,1000,20)) ),
                    'electron_pt0' : hist.Bin('pt',r'Leading electron $p_{T}$ (GeV)',list(range(0,600,20))),
                    'electron_pt1' : hist.Bin('pt',r'Trailing electron $p_{T}$ (GeV)',list(range(0,600,20))),
                    'electron_pt' : hist.Bin('pt',r'All electron $p_{T}$ (GeV)',list(range(0,600,20))),
                    'muon_pt0' : hist.Bin('pt',r'Leading muon $p_{T}$ (GeV)',list(range(0,600,20))),
                    'muon_pt1' : hist.Bin('pt',r'Trailing muon $p_{T}$ (GeV)',list(range(0,600,20))),
                    'muon_pt' : hist.Bin('pt',r'All muon $p_{T}$ (GeV)',list(range(0,600,20))),
                    'dielectron_pt' : hist.Bin('pt',r'Dielectron $p_{T}$ (GeV)',list(range(0,400,20)) + list(range(400,800,40)) + list(range(800,1100,100))),
                    'dimuon_pt' : hist.Bin('pt',r'Dimuon $p_{T}$ (GeV)',list(range(0,400,20)) + list(range(400,800,40)) + list(range(800,1100,100))),
                    'mjj' : hist.Bin('mjj', r'$M_{jj}$ (GeV)', list(range(200,800,300)) + list(range(800,2000,400)) + [2000, 2750, 3500]),
                    'softjet_pt' : hist.Bin("pt","Sum $p_{T}$ (GeV)", list(range(80,600,20))),
                    'softjet_hem_pt_sum' : hist.Bin("pt","Sum $p_{T}$ (GeV)", 20,0,400),
                    'softjet_anti_hem_pt_sum' : hist.Bin("pt","Sum $p_{T}$ (GeV)", 20,0,400),
                    'isotrack_hem_pt_sum' : hist.Bin("pt","Sum $p_{T}$ (GeV)", [0,10,20,30,40,50]),
                    'isotrack_anti_hem_pt_sum' : hist.Bin("pt","Sum $p_{T}$ (GeV)", [0,10,20,30,40,50]),
                    }
        }
        # binning for all monov regions:
        for region in ['cr_1m_v','cr_2m_v','cr_1e_v','cr_2e_v','cr_g_v','sr_v','cr_nobveto_v']:
            for wp in ['','_inclusive','_loose','_tight','_loosemd','_tightmd']:
                new_region_name = region.replace('_v', wp+'_v')
                self.binnings[new_region_name] = {
                        'recoil' : hist.Bin('recoil','Recoil (GeV)', recoil_monov_bins_2016),
                        'recoil_nodibosonnlo' : hist.Bin('recoil','Recoil (GeV)', recoil_monov_bins_2016),
                        }

    def get_binning(self, distribution, region='default'):
        if region in self.binnings and distribution in self.binnings[region]:
            return self.binnings[region][distribution]
        else:
            return self.binnings['default'][distribution]
def channel_name(region):
    if '_vbf' in region:
        return 'VBF'
    if '_j' in region:
        return 'Monojet'
    if '_v' in region:
        return 'Mono-V'

def make_plot(acc, region, distribution, year,  data, mc, signal=None, outdir='./output/stack/', integrate=None, ylim=None, xlim=None, rylim=None, tag=None, output_format='pdf', ratio=True):
    """Creates a data vs MC comparison plot

    :param acc: Accumulator (processor output)
    :type acc: coffea.processor.accumulator
    """
    # Rebin
    s = Style()
    h = copy.deepcopy(acc[distribution])
    assert(h)
    try:
        newax = s.get_binning(distribution, region)
        h = h.rebin(h.axis(newax.name), newax)
    except KeyError:
        pass

    # Integrate over an extra axis
    inte_tag=""
    if integrate:
        (inte_axis, inte_low, inte_high) = integrate
        h = h.integrate(inte_axis,slice(inte_low,inte_high)) #can add an overflow option here
        inte_tag+="_"+inte_axis+"_"+str(inte_low)+"_"+str(inte_high)

    # Pick the region we want to look at
    # E.g. cr_2m_j = Di-Muon control region with monojet selection
    h = h.integrate(h.axis('region'),region)

    # Plotting
    # Add ratio plot at the bottom if specified (default)
    # Otherwise just plot the histogram
    if ratio:
        fig, (ax, rax) = plt.subplots(2, 1, figsize=(7,7), gridspec_kw={"height_ratios": (3, 1)}, sharex=True)

    else:
        fig, ax = plt.subplots(1, 1, figsize=(7,5))

    data_err_opts = {
        'linestyle':'none',
        'marker': '.',
        'markersize': 10.,
        'color':'k',
        'elinewidth': 1,
    }
    signal_err_opts = {
        'linestyle':'-',
        'color':'crimson',
        'elinewidth': 1,
    }

    # Plot single muon data
    # Note the syntax we use to pick the data set
    if data:
        hist.plot1d(
            h[data],
            overlay='dataset',
            error_opts=data_err_opts,
            ax=ax,
            overflow='all',
            binwnorm=1)

    if signal:
        hist.plot1d(
            h[signal],
            overlay='dataset',
            error_opts=signal_err_opts,
            ax=ax,
            overflow='all',
            binwnorm=1,
            clear=False)

    # Plot MC background samples
    # Here we use a regular expression to match
    # data sets we want
    hist.plot1d(
        h[mc],
        overlay='dataset',
        stack=True,
        clear=False,
        overflow='all',
        ax=ax,
        binwnorm=1)

    # Apply correct colors to BG histograms
    handles, labels = ax.get_legend_handles_labels()
    new_labels = []
    for handle, label in zip(handles, labels):
        col = None
        for k, v in colors.items():
            if re.match(k, label):
                col = v
                break
        if col:
            handle.set_color(col)
            handle.set_linestyle('-')
            handle.set_edgecolor('k')

        l = None

        channel = channel_name(region)
        # Pick the proper legend labels for the channel
        if channel == 'VBF':
            legend_labels_to_use = legend_labels['VBF']
        elif channel in ['Monojet', 'Mono-V']:
            legend_labels_to_use = legend_labels['Monojet/Mono-V']

        # Add in the common labels
        legend_labels_to_use.update(legend_labels['Common'])

        for k, v in legend_labels_to_use.items():
            if re.match(k, label):
                l = v
        new_labels.append(l if l else label)

    # Update legend
    try:
        region_name = s.region_names[region]
    except KeyError:
        region_name = region
    ax.legend(title=region_name,ncol=2,handles=handles, labels=new_labels)

    # Ratio plot
    if data:
        hist.plotratio(h[data].integrate('dataset'), h[mc].integrate('dataset'),
                    ax=rax,
                    denom_fill_opts={},
                    guide_opts={},
                    unc='num',
                    overflow='all',
                    error_opts=data_err_opts
                    )

    ax.text(1., 0., distribution,
                fontsize=10,
                horizontalalignment='right',
                verticalalignment='bottom',
                transform=ax.transAxes
               )
    fig.text(0., 1., '$\\bf{CMS}$ internal',
                fontsize=14,
                horizontalalignment='left',
                verticalalignment='bottom',
                transform=ax.transAxes
               )

    fig.text(1., 1., f'{channel_name(region)}, {lumi(year):.1f} fb$^{{-1}}$ ({year})',
                fontsize=14,
                horizontalalignment='right',
                verticalalignment='bottom',
                transform=ax.transAxes
               )
    # Aesthetics
    ax.set_yscale("log")
    ax.set_ylabel('Events / GeV')
    plot_settings=style.plot_settings()
    if region in plot_settings.keys():
        plot_settings=plot_settings[region]
    if distribution in plot_settings.keys():
        plot_settings=plot_settings[distribution]
    if ylim:
        if ylim=="auto":
            width = np.diff([x for x in h.axes() if "dataset" not in x.name][0].edges())
            vmc = h[mc].integrate("dataset").values()[()] / width
            try:
                vdata = h[data].integrate("dataset").values()[()] / width
            except:
                vdata = vmc
            if signal:
                vsig = h[signal].integrate("dataset").values()[()] / width
            else:
                vsig = vmc


            ax.set_ylim(
                0.5*min([np.min(vmc[vmc>0]), np.min(vdata[vdata>0]), np.min(vsig[vsig>0])]),
                1e2*max([np.max(vmc), np.max(vdata), np.min(vsig)]),
            )

        else:
            ax.set_ylim(ylim[0],ylim[1])
    elif 'ylim' in plot_settings.keys():
        ax.set_ylim(plot_settings['ylim'])
    else:
        ax.set_ylim(1e-1,1e6)

    if xlim:
        ax.set_xlim(xlim[0],xlim[1])
    elif 'xlim' in plot_settings.keys():
        ax.set_xlim(plot_settings['xlim'])

    if ratio:
        if rylim:
            rax.set_ylim(*rylim)
        else:
            rax.set_ylim(0.5,1.5)
        loc1 = matplotlib.ticker.MultipleLocator(base=0.2)
        loc2 = matplotlib.ticker.MultipleLocator(base=0.1)
        rax.yaxis.set_major_locator(loc1)
        rax.yaxis.set_minor_locator(loc2)
        rax.grid(axis='y',which='minor',linestyle='--')
        rax.grid(axis='y',which='major',linestyle='--')
        rax.set_ylabel('Data / MC')
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    for form in output_format.split(','):
        outpath = pjoin(outdir, f"{region}_{distribution}{inte_tag}_{tag + '_' if tag else ''}{year}.{form}")
        fig.savefig(outpath)
        print(f"Saved plot file in {outpath}")
    plt.close('all')

def main():
    # The input is saved in individual *.coffea files
    # in the directory given here.
    indir = "./input/eff/gamma"

    # 'acc' is short for 'accumulator', which is the output
    # produced by a coffea processor. It behaves like a python dict,
    # so you can import it also in an interactive python shell to play with it
    # Note that 'acc_from_indir' uses caching, so subsequent readings will be much faster!
    acc = acc_from_dir(indir)

    # The make_plot function currently just makes a dimuon
    # mass plot for the dimuon control region.
    # TODO:
    #   * Make more flexible: More regions, more plots, etc
    #   * Selection of input processes is currently just hardcoded -> handle better!
    for year in [2017]:
        data = re.compile(f'SingleMuon_{year}')
        mc = re.compile(f'DY.*HT.*{year}')
        region='cr_2m_j'
        for distribution in ['recoil', 'dimuon_mass']:
            make_plot(copy.deepcopy(acc), region=region,distribution=distribution, year=year, data=data, mc=mc)

    for year in [2017, 2018]:
        data = re.compile(f'SingleMuon_{year}')
        mc = re.compile(f'W.*HT.*{year}')
        region='cr_1m_j'
        for distribution in ['recoil']:
            make_plot(copy.deepcopy(acc), region=region,distribution=distribution, year=year, data=data, mc=mc)



if __name__ == "__main__":
    main()
