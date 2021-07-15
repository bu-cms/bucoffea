#!/usr/bin/env python
import os
import sys
import re
import uproot
import numpy as np
import mplhep as hep

from matplotlib import pyplot as plt
from coffea.hist import poisson_interval
from bucoffea.plot.util import lumi, fig_ratio
from pprint import pprint
from tqdm import tqdm

pjoin = os.path.join

def get_region_label(region):
    mapping = {
        'signal' : 'VBF SR',
        'dielec' : r'VBF $Z(ee)$',
        'dimuon' : r'VBF $Z(\mu\mu)$',
        'singleel' : r'VBF $W(e\nu)$',
        'singlemu' : r'VBF $W(\mu\nu)$',
        'photon' : r'VBF $\gamma+jets$',
    }
    return mapping[region]

def plot_year_ratio(tree, outtag):
    '''
    Plot 2017 to 2018 ratio in data and total background.
    Takes the "shapes_prefit" tree from fit diagnostics file as an input.
    '''
    regions = [
        'signal',
        'dielec',
        'dimuon',
        'singleel',
        'singlemu',
        'photon',
    ]

    lumi_2017 = lumi(2017)
    lumi_2018 = lumi(2018)

    # Luminosity normalization factor
    lumi_factor = lumi_2018 / lumi_2017

    data_err_opts = {
        'color' : 'k',
        'linestyle': '',
    }

    outdir = f'./output/{outtag}/2017_over_2018'
    try:
        os.makedirs(outdir)
    except FileExistsError:
        pass

    for region in tqdm(regions):
        t_2017 = tree[f'vbf_2017_{region}']
        t_2018 = tree[f'vbf_2018_{region}']

        # Extract data and total background
        h_data_2017 = t_2017['data']
        h_bkg_2017 = t_2017['total_background']
        
        h_data_2018 = t_2018['data']
        h_bkg_2018 = t_2018['total_background']

        r_data = h_data_2017.yvalues / h_data_2018.yvalues * lumi_factor
        rerr_data = np.vstack([
            h_data_2017.yerrorslow / h_data_2018.yvalues * lumi_factor,
            h_data_2017.yerrorshigh / h_data_2018.yvalues * lumi_factor
        ])

        r_bkg = h_bkg_2017.values / h_bkg_2018.values  * lumi_factor
        rerr_bkg = np.sqrt(h_bkg_2017.variances) / h_bkg_2018.values * lumi_factor

        fig, ax, rax = fig_ratio()
        bins = h_bkg_2017.edges
        hep.histplot(r_data, 
            bins=bins, 
            ax=ax, 
            histtype='errorbar', 
            yerr=rerr_data, 
            label='Ratio in Data', 
            **data_err_opts)
        
        hep.histplot(r_bkg, 
            bins=bins, 
            ax=ax, 
            histtype='step', 
            yerr=rerr_bkg, 
            label='Ratio in Total Bkg')
        
        ax.set_xlim(0,5000)
        ax.set_xlabel(r'$M_{jj} \ (GeV)$')
        ax.set_ylabel(r'Ratio')
        ax.legend(title='2017 / 2018')

        ax.yaxis.set_ticks_position('both')

        ax.text(0,1,get_region_label(region),
            fontsize=14,
            ha='left',
            va='bottom',
            transform=ax.transAxes
        )

        ax.text(1,1,'Lumi Normalized',
            fontsize=14,
            ha='right',
            va='bottom',
            transform=ax.transAxes
        )
        
        # Plot ratio of ratios
        rr = r_data / r_bkg
        rr_err = rerr_data / r_bkg
        hep.histplot(rr,
            ax=rax,
            bins=bins,
            yerr=rr_err,
            histtype='errorbar',
            **data_err_opts
        )

        rax.set_xlabel(r'$M_{jj} \ (GeV)$')
        rax.set_xlim(0,5000)
        rax.set_ylabel('Data / Prediction')
        rax.set_ylim(0.5,1.5)
        rax.grid(True)

        # Stat + sys uncertainty on predictions
        # unity = np.ones_like(rr)
        rr_err_bkg = (r_data / r_bkg ** 2) * rerr_bkg

        opts = {"step": "post", "facecolor": (0, 0, 0, 0.3), "linewidth": 0}
        rax.fill_between(
            bins,
            1-np.r_[rr_err_bkg, rr_err_bkg[-1]],
            1+np.r_[rr_err_bkg, rr_err_bkg[-1]],
            **opts
        )
        rax.axhline(1, xmin=0, xmax=1, color='red')

        outpath = pjoin(outdir, f'{region}.pdf')
        outpath = pjoin(outdir, f'{region}.png')
        fig.savefig(outpath)
        plt.close(fig)

def main():
    # Input fit diagnostic file
    inpath = sys.argv[1]
    f = uproot.open(inpath)

    outtag = os.path.basename(os.path.dirname(inpath))
    
    # Tree with pre-fit shapes
    t = f['shapes_prefit']

    plot_year_ratio(t, outtag)

if __name__ == '__main__':
    main()