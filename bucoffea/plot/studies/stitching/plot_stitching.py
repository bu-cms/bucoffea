#!/usr/bin/env python
import copy
import os
import re

from matplotlib import pyplot as plt

from coffea import hist

from bucoffea.plot.util import (acc_from_dir, merge_datasets, merge_extensions,
                                scale_xs_lumi,klepto_load)

pjoin = os.path.join
def plot_ht_stitching(acc, tag, regex):
    outdir = './output/ht/'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    for dist in ['lhe_ht']:
        h=copy.deepcopy(acc[dist])
        h = merge_extensions(h, acc)
        scale_xs_lumi(h)
        h = merge_datasets(h)

        ax = hist.plot1d(
            h[re.compile(regex)],
            overlay='dataset',
            overflow='all',
            binwnorm=True)
        plt.yscale('log')
        plt.ylim(1e-3,1e6)

        for ext in ['pdf','png']:
            ax.figure.savefig(pjoin(outdir,f'{tag}_{dist}.{ext}'))
        ax.figure.clf()
def main():
    # acc = acc_from_dir("./input/")

    acc = klepto_load('../lo_vs_nlo/input/2020-01-21_correction_facs_updated')
    acc.load('nevents')
    acc.load("sumw")
    acc.load("sumw_pileup")
    acc.load('lhe_ht')
    plot_ht_stitching(acc, tag='wjet', regex='WJet.*LNu.*HT.*')
    plot_ht_stitching(acc, tag='dy', regex='DYJetsToLL.*HT.*')
    # # acc = acc_from_dir("./input/znunu")
    # plot_ht_stitching(acc, tag='znunu', regex='ZJets.*NuNu.*HT.*')
    # # acc = acc_from_dir("./input/gjets")
    # plot_ht_stitching(acc, tag='gjets', regex='GJets_HT.*')
    # plot_ht_stitching(acc, tag='qcd', regex='QCD.*HT.*')





if __name__ == "__main__":
    main()
