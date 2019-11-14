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

colors = {
    'WN*J.*' : '#feb24c',
    '.*DY.*' : '#ffffcc',
    '.*EWK.*V.*' : '#c6dbef',
    '.*Diboson.*' : '#4292c6',
    '.*TT.*' : '#6a51a3',
    '.*ST.*' : '#9e9ac8',
    '.*QCD.*' : '#08306b',
    '.*GJets_HT.*' : '#fc4e2a',
    '.*GJets_SM.*' : '#a76b51',
    'ZJetsToNuNu.*' : '#0050ec',
    'ZNuNuGJets_.*' : '#0050ec'
}
class Style():
    def __init__(self):
        self.region_names = {
            'cr_1m_j' : 'Single-$\mu$ CR, monojet',
            'cr_2m_j' : 'Di-$\mu$ CR, monojet',
            'cr_1e_j' : 'Single-e CR, monojet',
            'cr_2e_j' : 'Di-e CR, monojet',
            'cr_g_j' : 'Single-Photon CR, monojet',
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
        recoil_bins_2016 = [ 250.,  280.,  310.,  340.,  370.,  400.,  430.,  470.,  510., 550.,  590.,  640.,  690.,  740.,  790.,  840.,  900.,  960., 1020., 1090., 1160., 1250., 1400., 1600., 1800., 2000.]
        recoil_monov_bins_2016 = [250,300,350,400,500,600,750,1000]
        self.binnings = {
            'default': {
                    'dimuon_mass' : hist.Bin('dilepton_mass',r'M($\mu^{+}\mu^{-}$)',30,60,120),
                    'dielectron_mass' : hist.Bin('dilepton_mass',r'M($e^{+}e^{-}$)',30,60,120),
                    # 'recoil' : hist.Bin('recoil','recoil',list(range(250,300,50)) + list(range(300,500,50)) + list(range(500,1000,100)) + list(range(1000,2000,200))),
                    'recoil' : hist.Bin('recoil','Recoil (GeV)', recoil_bins_2016),
                    'recoil_nopog' : hist.Bin('recoil','Recoil (GeV)', recoil_bins_2016),
                    'recoil_nopu' : hist.Bin('recoil','Recoil (GeV)', recoil_bins_2016),
                    'recoil_nopref' : hist.Bin('recoil','Recoil (GeV)', recoil_bins_2016),
                    'met' : hist.Bin('met',r'$p_{T}^{miss}$ (GeV)',list(range(0,500,50)) + list(range(500,1000,100)) + list(range(1000,2000,250))),
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
                    'dielectron_pt' : hist.Bin('pt',r'Dielectron $p_{T}$ (GeV)',list(range(0,400,25)) + list(range(400,800,50)) + list(range(800,1100,100))),
                    'dimuon_pt' : hist.Bin('pt',r'Dimuon $p_{T}$ (GeV)',list(range(0,400,25)) + list(range(400,800,50)) + list(range(800,1100,100))),
                    'mjj' : hist.Bin('mjj', r'$M_{jj}$ (GeV)', list(range(200,800,300)) + list(range(800,2000,400)) + [2000, 2750, 3500])
                    }
        }
        # binning for all monov regions:
        for region in ['cr_1m_v','cr_2m_v','cr_1e_v','cr_2e_v','cr_g_v','sr_v','cr_nobveto_v']:
            for wp in ['','_inclusive','_loose','_tight','_loosemd','_tightmd']:
                new_region_name = region.replace('_v', wp+'_v')
                self.binnings[new_region_name] = {
                        'recoil' : hist.Bin('recoil','Recoil (GeV)', recoil_monov_bins_2016),
                        }

    def get_binning(self, distribution, region='default'):
        if region in self.binnings and distribution in self.binnings[region]:
            return self.binnings[region][distribution]
        else:
            return self.binnings['default'][distribution]


def make_plot(acc, region, distribution, year,  data, mc, signal=None, outdir='./output/stack/', integrate=None, ylim=None, xlim=None, rylim=None, tag=None, output_format='pdf'):
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
    if not region.startswith('sr'):
        fig, (ax, rax) = plt.subplots(2, 1, figsize=(7,7), gridspec_kw={"height_ratios": (3, 1)}, sharex=True)

    else:
        fig, ax = plt.subplots(1, 1, figsize=(7,5))

    data_err_opts = {
        'linestyle':'none',
        'marker': '.',
        'markersize': 10.,
        'color':'k',
        'elinewidth': 1,
        'emarker': '_'
    }
    signal_err_opts = {
        'linestyle':'none',
        'marker': '.',
        'markersize': 10.,
        'color':'r',
        'elinewidth': 1,
        'emarker': '_'
    }

    # Plot single muon data
    # Note the syntax we use to pick the data set
    if data:
        fig, ax, _ = hist.plot1d(
            h[data],
            overlay='dataset',
            error_opts=data_err_opts,
            ax=ax,
            overflow='all',
            binwnorm=True)

    if signal:
        fig, ax, _ = hist.plot1d(
            h[signal],
            overlay='dataset',
            error_opts=signal_err_opts,
            ax=ax,
            overflow='all',
            binwnorm=True)

    # Plot MC background samples
    # Here we use a regular expression to match
    # data sets we want
    _, _, primitives = hist.plot1d(
        h[mc],
        overlay='dataset',
        stack=True,
        clear=False,
        overflow='all',
        ax=ax,
        binwnorm=True)

    for name, ps in primitives.items():
        name = str(name)
        col = None
        for k, v in colors.items():
            if re.match(k, name):
                col = v
                break
        for item in ps:
            if col:
                item.set_facecolor(col)
            item.set_linestyle('-')
            item.set_edgecolor('k')
    # Legend
    try:
        region_name = s.region_names[region]
    except KeyError:
        region_name = region
    ax.legend(title=region_name,ncol=1)

    # Ratio plot
    if not region.startswith('sr'):
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
    fig.text(1., 1., f'{lumi(year)} fb$^{{-1}}$ ({year})',
                fontsize=14,
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
    # Aesthetics
    ax.set_yscale("log")
    ax.set_ylabel('Events / Bin width')
    plot_settings=style.plot_settings()
    if region in plot_settings.keys(): 
        plot_settings=plot_settings[region]
    if distribution in plot_settings.keys(): 
        plot_settings=plot_settings[distribution]
    if ylim:
        ax.set_ylim(ylim[0],ylim[1])
    elif 'ylim' in plot_settings.keys():
        ax.set_ylim(plot_settings['ylim'])
    else:
        ax.set_ylim(1e-1,1e6)

    if xlim:
        ax.set_xlim(xlim[0],xlim[1])
    elif 'xlim' in plot_settings.keys():
        ax.set_xlim(plot_settings['xlim'])
    
    if not region.startswith('sr'):
        if rylim:
            rax.set_ylim(*rylim)
        else:
            rax.set_ylim(0.75,1.25)
        loc1 = matplotlib.ticker.MultipleLocator(base=0.2)
        loc2 = matplotlib.ticker.MultipleLocator(base=0.1)
        rax.yaxis.set_major_locator(loc1)
        rax.yaxis.set_minor_locator(loc2)
        rax.grid(axis='y',which='minor',linestyle='--')
        rax.grid(axis='y',which='major',linestyle='--')
        rax.set_ylabel('Data / MC')
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    outpath = pjoin(outdir, f"{region}_{distribution}{inte_tag}_{tag + '_' if tag else ''}{year}.{output_format}")
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
