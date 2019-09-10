#!/usr/bin/env python

from matplotlib import pyplot as plt
import os
from coffea import hist
import numpy as np

def debug_plot_output(output, region='inclusive', outdir='out'):
    """Dump all histograms as PDF."""
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    for name in output.keys():
        if name.startswith("_"):
            continue
        # if any([x in name for x in ['sumw','cutflow','selected_events','kinematics','weights']]):
        #     continue
        try:
            if np.sum(output[name].values().values()) == 0:
                continue
        except:
            continue
        try:
            h = output[name].integrate("region",region)
        except:
            continue
        print(name)
        try:
            fig, ax, _ = hist.plot1d(
                h,
                overlay='dataset',
                overflow='all',
                )
        except:
            continue
        fig.suptitle(f'{region}, {name}')
        # ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_ylim(0.1, 1e8)
        fig.savefig(os.path.join(outdir, f"{region}_{name}.pdf"))
        plt.close(fig)
