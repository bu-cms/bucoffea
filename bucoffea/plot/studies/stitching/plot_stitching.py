#!/usr/bin/env python
import copy
import os
import re

from matplotlib import pyplot as plt

from coffea import hist

from bucoffea.plot.util import (acc_from_dir, merge_datasets, merge_extensions,
                                scale_xs_lumi)

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

        fig, ax, _ = hist.plot1d(
            h[re.compile(regex)],
            overlay='dataset',
            overflow='all',
            binwnorm=True)
        plt.yscale('log')
        plt.ylim(1e-3,1e6)

        fig.savefig(pjoin(outdir,f'{tag}_{dist}.pdf'))
def main():
    acc = acc_from_dir("./input/")
    # plot_ht_stitching(acc, tag='dy', regex='DY.*HT(?!(.*new.*))')
    # acc = acc_from_dir("./input/wjet")
    plot_ht_stitching(acc, tag='wjet', regex='W.*HT(?!(.*new.*))')
    # acc = acc_from_dir("./input/znunu")
    plot_ht_stitching(acc, tag='znunu', regex='ZJets.*NuNu.*HT(?!(.*new.*))')
    # acc = acc_from_dir("./input/gjets")
    plot_ht_stitching(acc, tag='gjets', regex='GJets.*HT(?!(.*new.*))')
    plot_ht_stitching(acc, tag='qcd', regex='QCD.*HT(?!(.*new.*))')





if __name__ == "__main__":
    main()
