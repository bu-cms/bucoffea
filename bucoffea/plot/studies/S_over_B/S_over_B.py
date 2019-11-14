#!/usr/bin/env python

import copy
import os
import re
from collections import defaultdict
from pprint import pprint
from klepto.archives import dir_archive

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

def S_over_B(acc,distribution,region,mc,signal,unc=0.05,outname="S_over_B.png", newax=None, cutlim=None):

    # Assuming files are extensions are merged and xs corrected
    h = copy.deepcopy(acc[distribution])
    assert(h)
    if newax:
        h = h.rebin(h.axis(newax.name), newax)
    binax = h.axis(distribution)

    h_mc = h[mc].integrate(h.axis('dataset'))
    h_signal = h[signal].integrate(h.axis('dataset'))

    values_mc = h_mc.values()[(region,)]
    values_signal = h_signal.values()[(region,)]

    cut_points = binax.edges()[:-1]
    nbins = len(cut_points)
    bins = list(range(nbins))

    if cutlim:
        new_points = []
        new_bins = []
        for ibin,ipoint in enumerate(cut_points):
            if ipoint >= cutlim[0] and ipoint <= cutlim[1]:
                new_points.append(ipoint)
                new_bins.append(ibin)
        cut_points = new_points
        bins = new_bins

    values_SB = np.zeros(len(bins))
    for ix,ibin in enumerate(bins):
        Nsig = sum(values_signal[ibin:])
        Nbkg = sum(values_mc[ibin:])
        values_SB[ix] = Nsig / max(1,np.sqrt(Nbkg + pow(unc*Nbkg,2))) 
        print(cut_points[ix],round(Nsig,2), round(Nbkg,2), round(values_SB[ix],2))

    fig = plt.figure()
    plt.plot(cut_points, values_SB)
    fig.savefig('output/'+outname)


def main():
    inpath = "../../input/merged"
    year=2017
    mc = re.compile(f'(VDY.*HT.*|QCD.*|W.*HT.*|ST_|TTJets-FXFX_|Diboson_|GJets.*HT.*|ZJetsToNuNu.*){year}')
    signal=re.compile(f'WH.*{year}')
    distribution="recoil"
    acc = dir_archive(
        inpath,
        serialized=True,
        compression=0,
        memsize=1e3,
        )
    acc.load(distribution)
    acc.load('sumw')
    acc.load('sumw_pileup')
    acc.load('nevents')
    try:
        acc[distribution] = merge_extensions(acc[distribution], acc, reweight_pu=not ('nopu' in distribution))
        scale_xs_lumi(acc[distribution])
        acc[distribution] = merge_datasets(acc[distribution])
        S_over_B(acc, distribution, 'sr_tight_v',mc=mc,signal=signal,unc=0.05,outname="SB_unc005.png",cutlim=(250,750))
        S_over_B(acc, distribution, 'sr_tight_v',mc=mc,signal=signal,unc=0.10,outname="SB_unc010.png",cutlim=(250,750))
    except KeyError:
        print("key error ")
        return -2

if __name__ == "__main__":
    main()
