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
    '.*GJet.*' : '#fc4e2a',
}
class Style():
    def __init__(self):
        self.region_names = {
            'cr_1m_j' : 'Single-$\mu$ CR, monojet',
            'cr_2m_j' : 'Di-$\mu$ CR, monojet',
            'cr_1e_j' : 'Single-e CR, monojet',
            'cr_2e_j' : 'Di-e CR, monojet',
            'cr_g_j' : 'Single-Photon CR, monojet',
        }
        recoil_bins_2016 = [ 250.,  280.,  310.,  340.,  370.,  400.,  430.,  470.,  510., 550.,  590.,  640.,  690.,  740.,  790.,  840.,  900.,  960., 1020., 1090., 1160., 1250., 1400., 1600., 1800., 2000.]
        self.rebin_axes = {
            'dimuon_mass' : hist.Bin('dilepton_mass','dilepton_mass',30,60,120),
            'dielectron_mass' : hist.Bin('dilepton_mass','dilepton_mass',30,60,120),
            # 'recoil' : hist.Bin('recoil','recoil',list(range(250,300,50)) + list(range(300,500,50)) + list(range(500,1000,100)) + list(range(1000,2000,200))),
            'recoil' : hist.Bin('recoil','Recoil (GeV)', recoil_bins_2016),
            'recoil_nopog' : hist.Bin('recoil','Recoil (GeV)', recoil_bins_2016),
            'recoil_nopu' : hist.Bin('recoil','Recoil (GeV)', recoil_bins_2016),
            'recoil_nopref' : hist.Bin('recoil','Recoil (GeV)', recoil_bins_2016),
            'met' : hist.Bin('met','met',list(range(0,500,50)) + list(range(500,1000,100)) + list(range(1000,2000,250))),
            'ak4_pt0' : hist.Bin('jetpt','jetpt',list(range(100,600,20)) + list(range(600,1000,20)) ),
            'ak4_pt' : hist.Bin('jetpt','jetpt',list(range(100,600,20)) + list(range(600,1000,20)) ),
            'ak4_ptraw0' : hist.Bin('jetpt','jetpt',list(range(100,600,20)) + list(range(600,1000,20)) ),
            'ak4_pt0_eta0' : hist.Bin('jetpt','jetpt',list(range(100,600,20)) + list(range(600,1000,20)) ),
            'photon_pt0' : hist.Bin('pt','pt',list(range(200,600,20)) + list(range(600,1000,20)) ),
            'electron_pt0' : hist.Bin('pt','pt',list(range(0,600,20))),
            'electron_pt1' : hist.Bin('pt','pt',list(range(0,600,20))),
            'electron_pt' : hist.Bin('pt','pt',list(range(0,600,20))),
            'muon_pt0' : hist.Bin('pt','pt',list(range(0,600,20))),
            'muon_pt1' : hist.Bin('pt','pt',list(range(0,600,20))),
            'muon_pt' : hist.Bin('pt','pt',list(range(0,600,20))),
            'dielectron_pt' : hist.Bin('dilepton_pt','dilepton_pt',list(range(0,400,25)) + list(range(400,800,50)) + list(range(800,1100,100))),
            'dimuon_pt' : hist.Bin('dilepton_pt','dilepton_pt',list(range(0,400,25)) + list(range(400,800,50)) + list(range(800,1100,100))),
        }


def make_plot(acc, region, distribution, year,  data, mc, outdir='./output/stack/', integrate=None, ylim=None, xlim=None, rylim=None, tag=None, output_format='pdf'):
    """Creates a data vs MC comparison plot

    :param acc: Accumulator (processor output)
    :type acc: coffea.processor.accumulator
    """
    # Rebin
    s = Style()
    h=copy.deepcopy(acc[distribution])
    try:
        newax = s.rebin_axes[distribution]
        h = h.rebin(h.axis(newax.name), newax)
    except KeyError:
        pass

    # Integrate over an extra axis
    inte_tag=""
    if integrate:
        (inte_axis, inte_low, inte_high) = integrate
        h = h.integrate(inte_axis,slice(inte_low,inte_high)) #can add an overflow option here
        inte_tag+="_"+inte_axis+"_"+str(inte_low)+"_"+str(inte_high)


    # Step 1: Merge extension samples together
    # Extension samples are separate MC samples for identical physics processes
    # (E.g. when people realize that they need more events for an existing sample,
    # they produce an 'extension')
    h = merge_extensions(h, acc, reweight_pu=('nopu' in distribution))

    # Step 2: Scale each dataset according to its cross section
    # and the luminosity for the corresponding year
    # TODO: We may need this to be more flexible for control regions
    # If a CR has a prescaled trigger, lumi will not be total year lumi
    scale_xs_lumi(h)

    # Step 3: Merge datasets together
    # E.g. merge all ZToNuNu HT binned samples into just one sample
    h = merge_datasets(h)

    # Step 4: Pick the region we want to look at
    # E.g. cr_2m_j = Di-Muon control region with monojet selection
    h = h.integrate(h.axis('region'),region)


    # Plotting

    fig, (ax, rax) = plt.subplots(2, 1, figsize=(7,7), gridspec_kw={"height_ratios": (3, 1)}, sharex=True)

    data_err_opts = {
        'linestyle':'none',
        'marker': '.',
        'markersize': 10.,
        'color':'k',
        'elinewidth': 1,
        'emarker': '_'
    }

    # Plot single muon data
    # Note the syntax we use to pick the data set
    fig, ax, _ = hist.plot1d(
        h[data],
        overlay='dataset',
        error_opts=data_err_opts,
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
    if ylim:
        ax.set_ylim(ylim[0],ylim[1])
    else:
        ax.set_ylim(1e-1,1e6)
    if xlim:
        ax.set_xlim(xlim[0],xlim[1])

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
    ax.set_ylabel('Events / Bin width')
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
