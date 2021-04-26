#!/usr/bin/env python

import os
import sys
import re
import numpy as np
import uproot

from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
from bucoffea.helpers.paths import bucoffea_path
from bucoffea.plot.trigger import get_xy, ratio_unc
from bucoffea.plot.util import fig_ratio
from bucoffea.plot.style import matplotlib_rc
from matplotlib.ticker import MultipleLocator
from numpy import linalg

pjoin = os.path.join

matplotlib_rc()

def get_pretty_legend_label(tag):
    pretty_legend_label = {
        'two_central_jets' : 'Two Central Jets',
        'one_jet_forward_one_jet_central' : 'Mixed',
        'inclusive_nohfhf' : 'Inclusive (No HF-HF)'
    }
    return pretty_legend_label[tag]

def check_files(fnum, fden):
    if not os.path.exists(fnum):
        raise RuntimeError(f"File not found {fnum}")
    if not os.path.exists(fden):
        raise RuntimeError(f"File not found {fden}")

def sigmoid(x,a,b,c):
    return c / (1 + np.exp(-a * (x-b)) )

def do_fit(xsf, ysf, ysferr):
    '''Fit a sigmoid function to scale factor data.'''

    # Initial guess
    p0 = (0.05, 200, 1)

    popt, pcov = curve_fit(
        sigmoid,
        xsf, ysf,
        p0=p0
    )

    # Diagonalize the covariance matrix + calculate fit variations
    fit_params = {
        'fit' : popt
    }
    pcov_w, pcov_v = linalg.eig(pcov)
    for idx in range(len(pcov_w)):
        variation = np.sign(pcov_w[idx]) * np.sqrt(np.abs(pcov_w[idx])) * pcov_v[:,idx]
        fit_params[f'fit_{idx}_up'] = popt + variation
        fit_params[f'fit_{idx}_dn'] = popt - variation

    return fit_params

def fit_efficiencies(fdata, fmc, jeteta_config, year, outputrootfile, outdir):
    '''Fit the efficiency with a sigmoid function.'''
    x_data, x_edg_data, y_data, yerr_data = get_xy(fdata)
    x_mc, _, y_mc, yerr_mc = get_xy(fmc)

    fig, ax, rax = fig_ratio()
    # Plot original data and MC
    ax.errorbar(x_data, y_data, yerr=yerr_data, marker='o', ls='', label='Data')
    ax.errorbar(x_mc, y_mc, yerr=yerr_mc, marker='o', ls='', label='MC')

    # Get the sigmoid fit for data and MC
    fit_params_data = do_fit(x_data, y_data, yerr_data)
    f_data = sigmoid(x_data, *fit_params_data['fit'])
    ax.plot(x_data, f_data, lw=2, label='Data fit')

    fit_params_mc = do_fit(x_mc, y_mc, yerr_mc)
    f_mc = sigmoid(x_mc, *fit_params_mc['fit'])
    ax.plot(x_mc, f_mc, lw=2, label='MC fit', ls='--')

    ax.set_ylabel('Trigger Efficiency')
    ax.legend()

    # Show the fit parameter values on the plots
    num_params = len(fit_params_data['fit'])
    ax.text(0.25, 0.35, 'Data fit:', transform=ax.transAxes)
    labels = ['a', 'b', 'c']
    for idx in range(num_params):
        x = 0.25 
        y = 0.35 - (idx+1) * 0.08
        ax.text(x,y, f'{labels[idx]}: {fit_params_data["fit"][idx]:.3f}', transform=ax.transAxes)

    ax.text(0.48, 0.35, 'MC fit:', transform=ax.transAxes)
    labels = ['a', 'b', 'c']
    for idx in range(num_params):
        x = 0.48
        y = 0.35 - (idx+1) * 0.08
        ax.text(x,y, f'{labels[idx]}: {fit_params_mc["fit"][idx]:.3f}', transform=ax.transAxes)

    # Calculate scale factors
    sf_orig = y_data / y_mc
    sf_fit = f_data / f_mc

    sf_orig_err = yerr_data / y_mc

    opts = {
        'markersize' : 12.,
        'linestyle' : 'none',
        'marker' : '.',
        'color' : 'C1'
    }

    # Plot the original and fitted SF
    rax.errorbar(x_data, sf_orig, yerr=sf_orig_err, label='Data/MC', **opts)
    rax.plot(x_data, sf_fit, label='Data/MC Fit Ratio', color='k', lw=2)

    rax.set_xlabel('Recoil (GeV)')
    rax.set_ylabel('Data / MC SF')
    rax.grid(True)
    rax.set_ylim(0.95,1.05)
    rax.legend(prop={'size': 10.})

    loc1 = MultipleLocator(0.05)
    loc2 = MultipleLocator(0.01)
    rax.yaxis.set_major_locator(loc1)
    rax.yaxis.set_minor_locator(loc2)

    ylim = rax.get_ylim()
    rax.plot([250., 250.], ylim, color='red', lw=2)
    rax.set_ylim(ylim)

    plt.text(0., 1., f'{get_pretty_legend_label(jeteta_config)} {year}',
        fontsize=16,
        horizontalalignment='left',
        verticalalignment='bottom',
        transform=ax.transAxes
        )

    outpath = pjoin(outdir, f'eff_fit_{jeteta_config}_{year}.pdf')
    fig.savefig(outpath)
    plt.close(fig)
    print(f'File saved: {outpath}')

    x_edges = np.array(list(x_edg_data[0]) + [x_edg_data[1][-1]])

    outputrootfile[f'sf_{jeteta_config}_{year}'] = (sf_fit, x_edges)

def compare_with_current_sf(new_rootfile, outdir):
    '''Plot the newly calculated SF and compare with the original one.'''
    # Get the current scale factors for comparison
    orig_sf_file = bucoffea_path('data/sf/trigger/met_trigger_sf.root')
    f_orig = uproot.open(orig_sf_file)

    # Get the newly calculated SFs
    f_new = uproot.open(new_rootfile)

    for year in [2017, 2018]:
        fig, ax, rax = fig_ratio()

        old_sf = f_orig[f'120pfht_hltmu_1m_{year}'].values
        old_sf_edges = f_orig[f'120pfht_hltmu_1m_{year}'].edges

        old_sf_centers = ((old_sf_edges + np.roll(old_sf_edges, -1)) / 2)[:-1]

        # Take values where recoil > 250 GeV
        recoilmask =  old_sf_centers > 240
        old_sf = old_sf[recoilmask]
        old_sf_centers = old_sf_centers[recoilmask]

        ax.plot(old_sf_centers, old_sf, marker='o', ls='', label='Current SF', color='k')

        ratios = {}

        opts = {
            'markersize' : 12.,
            'linestyle' : 'none',
            'marker' : '.'
        }

        jeteta_configs = ['two_central_jets', 'one_jet_forward_one_jet_central', 'inclusive_nohfhf']
        for idx, jeteta_config in enumerate(jeteta_configs):
            new_sf = f_new[f'sf_{jeteta_config}_{year}'].values
            new_sf_edges = f_new[f'sf_{jeteta_config}_{year}'].edges

            new_sf_centers = ((new_sf_edges + np.roll(new_sf_edges, -1)) / 2)[:-1]

            new_sf = new_sf[recoilmask]
            new_sf_centers = new_sf_centers[recoilmask]

            ax.plot(new_sf_centers, new_sf, label=get_pretty_legend_label(jeteta_config), color=f'C{idx}', lw=2, marker='o')

            # Ratio with the current SF
            r = new_sf / old_sf

            rax.plot(new_sf_centers, r, color=f'C{idx}', **opts)

        ax.legend()
        ax.set_ylabel('Data / MC SF')

        loc1 = MultipleLocator(0.01)
        loc2 = MultipleLocator(0.005)
        ax.yaxis.set_major_locator(loc1)
        ax.yaxis.set_minor_locator(loc2)

        rax.grid(True)
        rax.set_ylim(0.95,1.05)
        rax.set_xlabel('Recoil (GeV)')
        rax.set_ylabel('Ratio to current SF')
            
        loc3 = MultipleLocator(0.01)
        rax.yaxis.set_minor_locator(loc3)

        plt.text(0., 1., year,
            fontsize=16,
            horizontalalignment='left',
            verticalalignment='bottom',
            transform=ax.transAxes
            )


        # Save figure
        outpath = pjoin(outdir, f'sf_comparison_{year}.pdf')
        fig.savefig(outpath)
        plt.close(fig)
        print(f'File saved: {outpath}')

def main():
    # Input directory to read txt files from
    input_dir = './output/120pfht_mu_recoil/merged_2020-10-20_vbfhinv_03Sep20v7_trigger_study'
    outtag = input_dir.split('/')[-1]

    jeteta_configs = [
        'two_central_jets',
        'one_jet_forward_one_jet_central',
        'inclusive_nohfhf'
    ]

    # Save the SFs to an output ROOT file
    outdir = f'./output/120pfht_mu_recoil/{outtag}/sigmoid_fit/with_normalization'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    outputrootpath = pjoin(outdir, 'fitted_sf.root')
    outputrootfile = uproot.recreate(outputrootpath)

    for year in [2017,2018]:
        for jeteta_config in jeteta_configs:
            # Figure out the input txt files
            f_data = pjoin(input_dir, f'table_1m_recoil_SingleMuon_{year}_{jeteta_config}.txt')
            f_mc = pjoin(input_dir, f'table_1m_recoil_WJetsToLNu_HT_MLM_{year}_{jeteta_config}.txt')
        
            # Check the file paths
            check_files(f_data, f_mc)

            fit_efficiencies(f_data, f_mc, 
                jeteta_config=jeteta_config, 
                year=year, 
                outdir=outdir,
                outputrootfile=outputrootfile)

    compare_with_current_sf(outputrootpath, outdir)

if __name__ == '__main__':
    main()