#!/usr/bin/env python

import os
import sys
import re
import numpy as np

from coffea import hist
from bucoffea.plot.util import merge_datasets, merge_extensions, scale_xs_lumi
from matplotlib import pyplot as plt
from klepto.archives import dir_archive

pjoin = os.path.join

def plot_dphi_hfjet_met(acc, outtag, region='sr_vbf_fail_hf_cuts', distribution='dphi_hfjet_met'):
    acc.load(distribution)
    h = acc[distribution]

    h = merge_extensions(h,acc,reweight_pu=False)
    scale_xs_lumi(h)
    h = merge_datasets(h)

    h = h.integrate('dataset', 'MET_2017').integrate('region', region)

    fig, ax = plt.subplots()
    hist.plot1d(h, ax=ax)

    ax.set_yscale('log')
    ax.set_ylim(1e-1,1e4)
    ax.set_xlabel(r'$\Delta\phi(HF\ jet, MET)$')

    ax.get_legend().remove()

    ax.text(0.,1.,'ULv8 MET 2017',
        fontsize=14,
        ha='left',
        va='bottom',
        transform=ax.transAxes
    )

    ax.text(1.,1.,'QCD CR (HF cuts inv.)',
        fontsize=14,
        ha='right',
        va='bottom',
        transform=ax.transAxes
    )

    # Calculate the % events where dphi(j,MET) > 2.5
    mask = h.axis('dphi').centers() > 2.5
    total_events = np.sum(h.values()[()])
    high_dphi_events = np.sum(h.values()[()][mask])
    r = high_dphi_events / total_events * 100

    ax.text(1., 0.9, f'$\\Delta\\phi(j_{{HF}}, MET) > 2.5$: {r:.2f}%',
        fontsize=12,
        ha='right',
        va='bottom',
        transform=ax.transAxes
    )

    outdir = f'./output/{outtag}'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outpath = pjoin(outdir, f'{distribution}.pdf')
    fig.savefig(outpath)
    plt.close(fig)

    print(f'File saved: {outpath}')

def main():
    inpath = sys.argv[1]
    acc = dir_archive(inpath)
    acc.load('sumw')
    acc.load('sumw2')

    outtag = re.findall('merged_.*', inpath)[0].replace('/','')

    plot_dphi_hfjet_met(acc, outtag)

if __name__ == '__main__':
    main()