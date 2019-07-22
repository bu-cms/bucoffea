#!/usr/bin/env python

import os
import re
from collections import defaultdict

import numpy as np
from coffea import hist
from coffea.util import load
from matplotlib import pyplot as plt

from bucoffea.execute.dataset_definitions import short_name
from bucoffea.helpers.dataset import is_data
from bucoffea.helpers.paths import bucoffea_path
from bucoffea.plot.util import (acc_from_dir, lumi, merge_datasets,
                                merge_extensions, normalize_mc)

pjoin = os.path.join

def make_plot(acc):
    """Creates a data vs MC comparison plot

    :param acc: Accumulator (processor output)
    :type acc: coffea.processor.accumulator
    """
    year = 2018

    # Rebin
    h=acc['dimuon_mass']
    # h=h.rebin(h.axis('recoil'),hist.Bin('recoil','recoil',25,0,1000))

    # The histograms are binned in the dataset name
    # Step 1: Only keep datasets that correspond to the right year
    mapping ={str(k):str(k) for k in h.axis('dataset').identifiers() if str(year) in str(k)}
    h=h.group(h.axis('dataset'), 'dataset', mapping)

    # Step 2: Merge extension samples together
    # Extension samples are separate MC samples for identical physics processes
    # (E.g. when people realize that they need more events for an existing sample,
    # they produce an 'extension')
    h = merge_extensions(h)

    # Step 3: Scale each dataset according to its cross section
    # and the luminosity for the corresponding year
    normalize_mc(h, acc)

    # Step 4: Merge datasets together
    # E.g. merge all ZToNuNu HT binned samples into just one sample
    h = merge_datasets(h)

    # Step 5: Pick the region we want to look at
    # E.g. cr_2m_j = Di-Muon control region with monojet selection
    h = h.project(h.axis('region'),'cr_2m_j')


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
    fig, ax, _ = hist.plot1d(
        h['SingleMuon_2018'],
        overlay='dataset',
        error_opts=data_err_opts,
        ax=ax,
        overflow='over')

    # Plot MC background samples
    mc=re.compile("DYNJets.*")
    hist.plot1d(
        h[mc],
        overlay='dataset',
        stack=True,
        clear=False,
        overflow='over',
        ax=ax)

    # Legend
    leg = ax.legend(title='Di-$\mu$ CR',ncol=1)

    # Ratio plot
    data_err_opts = {
        'linestyle':'none',
        'marker': '.',
        'markersize': 10.,
        'color':'k',
        'elinewidth': 1,
        'emarker': '_'
    }

    hist.plotratio(h['SingleMuon_2018'].project('dataset'), h[mc].project('dataset'),
                ax=rax,
                denom_fill_opts={},
                guide_opts={},
                unc='num',
                error_opts=data_err_opts
                )

    # Aesthetics
    # ax.set_xlim(60,120)
    ax.set_yscale("log")
    ax.set_ylim(1e-1,1e6)
    rax.set_ylim(0.5,1.5)
    fig.savefig("dimu_mass.pdf")

def main():
    # The input is saved in individual *.coffea files
    # in the directory given here.
    indir = "./input/eff/120pfht_hltmu"

    # 'acc' is short for 'accumulator', which is the output
    # produced by a coffea processor. It behaves like a python dict,
    # so you can import it also in an interactive python shell to play with it
    acc = acc_from_dir(indir)

    # The make_plot function currently just makes a dimuon
    # mass plot for the dimuon control region.
    # TODO:
    #   * Make more flexible: More regions, more plots, etc
    #   * Selection of input processes is currently just hardcoded -> handle better!
    make_plot(acc)



if __name__ == "__main__":
    main()
