#!/usr/bin/env python
import argparse
import os
import re
import sys
import uproot
import numpy as np
import mplhep as hep

from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
from coffea import hist
from coffea.hist import poisson_interval
from bucoffea.plot.util import merge_datasets, merge_extensions, scale_xs_lumi, fig_ratio, lumi
from bucoffea.helpers.paths import bucoffea_path
from klepto.archives import dir_archive
from pprint import pprint

pjoin = os.path.join

# Suppress true_divide warnings
np.seterr(divide='ignore', invalid='ignore')
import matplotlib.pylab as pylab
params = {'legend.fontsize': 'small',
        #   'figure.figsize': (15, 5),
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
pylab.rcParams.update(params)

Bin = hist.Bin

distributions = [
    'mjj',
    'ak4_eta0',
    'ak4_eta1',
    'ak4_phi0',
    'ak4_phi1',
    'ak4_pt0',
    'ak4_pt1',
    'ak4_nef0',
    'ak4_nef1',
    'ak4_nhf0',
    'ak4_nhf1',
    'ak4_chf0',
    'ak4_chf1',
    'ak4_mt0',
    'ak4_mt1',
    'vecb',
    'vecdphi',
    'dphitkpf',
    'met',
    'met_phi',
    'ak4_mult',
    'ak4_mt0',
    'ak4_mt1',
]

binnings = {
    'mjj': Bin('mjj', r'$M_{jj} \ (GeV)$', [200., 400., 600., 900., 1200., 1500., 2000., 2750., 3500.]),
    'ak4_pt0': Bin('jetpt',r'Leading AK4 jet $p_{T}$ (GeV)',list(range(80,600,20)) + list(range(600,1000,20)) ),
    'ak4_pt1': Bin('jetpt',r'Trailing AK4 jet $p_{T}$ (GeV)',list(range(40,600,20)) + list(range(600,1000,20)) ),
    'ak4_phi0' : Bin("jetphi", r"Leading AK4 jet $\phi$", 50,-np.pi, np.pi),
    'ak4_phi1' : Bin("jetphi", r"Trailing AK4 jet $\phi$", 50,-np.pi, np.pi),
    'ak4_nef0' : Bin('frac', 'Leading Jet Neutral EM Frac', 50, 0, 1),
    'ak4_nef1' : Bin('frac', 'Trailing Jet Neutral EM Frac', 50, 0, 1),
    'ak4_nhf0' : Bin('frac', 'Leading Jet Neutral Hadronic Frac', 50, 0, 1),
    'ak4_nhf1' : Bin('frac', 'Trailing Jet Neutral Hadronic Frac', 50, 0, 1),
    'ak4_chf0' : Bin('frac', 'Leading Jet Charged Hadronic Frac', 50, 0, 1),
    'ak4_chf1' : Bin('frac', 'Trailing Jet Charged Hadronic Frac', 50, 0, 1),
    'dphitkpf' : Bin('dphi', r'$\Delta\phi_{TK,PF}$', 50, 0, 3.5),
    'met' : Bin('met',r'$p_{T}^{miss}$ (GeV)',list(range(0,500,50)) + list(range(500,1000,100)) + list(range(1000,2000,250))),
    'met_phi' : Bin("phi", r"$\phi_{MET}$", 50, -np.pi, np.pi),
    'ak4_mult' : Bin("multiplicity", r"AK4 multiplicity", 10, -0.5, 9.5),
    'ak4_mt0' : Bin("mt", r"Leading AK4 $M_{T}$ (GeV)", 50, 0, 1000),
    'ak4_mt1' : Bin("mt", r"Trailing AK4 $M_{T}$ (GeV)", 50, 0, 1000),
}

ylims = {
    'ak4_eta0' : (1e-3,1e8),
    'ak4_eta1' : (1e-3,1e8),
    'ak4_nef0' : (1e0,1e8),
    'ak4_nef1' : (1e0,1e8),
    'ak4_nhf0' : (1e0,1e8),
    'ak4_nhf1' : (1e0,1e8),
    'ak4_chf0' : (1e0,1e8),
    'ak4_chf1' : (1e0,1e8),
    'vecb' : (1e-1,1e9),
    'vecdphi' : (1e0,1e9),
    'dphitkpf' : (1e0,1e9),
    'met' : (1e-3,1e5),
    'ak4_mult' : (1e-1,1e8),
}

legend_labels = {
    'GJets_DR-0p4.*' : "QCD $\\gamma$+jets",
    'GJets_SM.*' : "EWK $\\gamma$+jets",
    'DY.*' : "QCD Z$\\rightarrow\\ell\\ell$",
    'EWKZ.*ZToLL.*' : "EWK Z$\\rightarrow\\ell\\ell$",
    'WN*J.*LNu.*' : "QCD W$\\rightarrow\\ell\\nu$",
    'EWKW.*LNu.*' : "EWK W$\\rightarrow\\ell\\nu$",
    'ZJetsToNuNu.*.*' : "QCD Z$\\rightarrow\\nu\\nu$",
    'EWKZ.*ZToNuNu.*' : "EWK Z$\\rightarrow\\nu\\nu$",
    'QCD.*' : "QCD Estimation",
    'Top.*' : "Top quark",
    'Diboson.*' : "WW/WZ/ZZ",
    'MET|Single(Electron|Photon|Muon)|EGamma.*' : "Data",
    'VBF_HToInv.*' : "VBF H(inv)",
}

legend_titles = {
    'sr_vbf' : 'VBF Signal Region',
    'cr_1m_vbf' : r'VBF $1\mu$ Region',
    'cr_2m_vbf' : r'VBF $2\mu$ Region',
    'cr_1e_vbf' : r'VBF $1e$ Region',
    'cr_2e_vbf' : r'VBF $2e$ Region',
    'cr_g_vbf' : r'VBF $\gamma$ Region',
}

colors = {
    'DY.*' : '#ffffcc',
    'EWKW.*' : '#c6dbef',
    'EWKZ.*ZToLL.*' : '#d5bae2',
    'EWKZ.*ZToNuNu.*' : '#c4cae2',
    '.*Diboson.*' : '#4292c6',
    'Top.*' : '#6a51a3',
    '.*HF Noise.*' : '#08306b',
    '.*TT.*' : '#6a51a3',
    '.*ST.*' : '#9e9ac8',
    'ZJetsToNuNu.*' : '#31a354',
    'WJets.*' : '#feb24c',
}

def plot_data_mc(acc, outtag, year, data, mc, data_region, mc_region, distribution='mjj', plot_signal=True, mcscale=1, fformat='pdf', qcd_file=None):
    '''Plot data/MC comparison with the QCD template included.'''
    acc.load(distribution)
    h = acc[distribution]

    if distribution == 'mjj':
        overflow = 'over'
    else:
        overflow = 'none'

    h = merge_extensions(h, acc, reweight_pu=False)
    scale_xs_lumi(h)
    h = merge_datasets(h)

    if distribution in binnings.keys():
        new_ax = binnings[distribution]
        h = h.rebin(new_ax.name, new_ax)

    h.axis('dataset').sorting = 'integral'

    h.scale({
        ds : (mcscale  if mc.match(ds) else 1) for ds in map(str,h.axis("dataset").identifiers())
    }, axis='dataset')

    h_data = h.integrate('region', data_region)
    h_mc = h.integrate('region', mc_region)

    # Get the QCD template (estimation from HF), only to be used in SR for now
    if data_region == 'sr_vbf':
        # If a path to QCD file has been given, use it!
        # Otherwise, take the one from the relevant output directory
        if qcd_file:
            qcdfilepath = qcd_file
        else:    
            qcdfilepath = f'output/{outtag}/qcd_estimate/hf_qcd_estimate.root'    
        h_qcd = uproot.open(qcdfilepath)[f'qcd_estimate_{distribution}_{year}']

    data_err_opts = {
        'linestyle':'none',
        'marker': '.',
        'markersize': 10.,
        'color':'k',
        'elinewidth': 1,
    }

    datasets = list(map(str, h[mc].identifiers('dataset')))

    plot_info = {
        'label' : datasets,
        'sumw' : [],
    }

    # Temporary fix for sorting
    if data_region == 'sr_vbf':
        if 'WJetsToLNu' in datasets[-1]:
            tmp = datasets[-2]
            datasets[-2] = datasets[-1]
            datasets[-1] = tmp

    for dataset in datasets:
        sumw = h_mc.integrate('dataset', dataset).values(overflow=overflow)[()]

        plot_info['sumw'].append(sumw)

    # Add the QCD contribution (for SR only)
    if data_region == 'sr_vbf':
        plot_info['label'].insert(6, 'HF Noise Estimate')
        plot_info['sumw'].insert(6, h_qcd.values * mcscale)

    fig, ax, rax = fig_ratio()
    hist.plot1d(h_data[data], ax=ax, overflow=overflow, overlay='dataset', binwnorm=1, error_opts=data_err_opts)

    xedges = h_data.integrate('dataset').axes()[0].edges(overflow=overflow)

    hep.histplot(plot_info['sumw'], xedges, 
        ax=ax,
        label=plot_info['label'], 
        histtype='fill',
        binwnorm=1,
        stack=True
        )

    if plot_signal:
        signal = re.compile(f'VBF_HToInvisible.*withDipoleRecoil.*{year}')

        signal_line_opts = {
            'linestyle': '-',
            'color': 'crimson',
        }

        h_signal = h.integrate('region', mc_region)[signal]
        h_signal.scale(mcscale)

        hist.plot1d(
            h_signal,
            ax=ax,
            overlay='dataset',
            overflow=overflow,
            line_opts=signal_line_opts,
            binwnorm=1,
            clear=False
        )

    ax.set_yscale('log')
    if distribution == 'mjj':
        ax.set_ylim(1e-3,1e5)
        ax.set_ylabel('Events / GeV')
    else:
        if distribution in ylims.keys():
            ax.set_ylim(ylims[distribution])
        else:
            ax.set_ylim(1e-2,1e6)
        ax.set_ylabel('Events / Bin Width')
    
    if distribution == 'mjj':
        ax.set_xlim(left=0.)

    ax.yaxis.set_ticks_position('both')

    handles, labels = ax.get_legend_handles_labels()
    for handle, label in zip(handles, labels):
        for datasetregex, new_label in legend_labels.items():
            col = None
            if re.match(datasetregex, label):
                handle.set_label(new_label)
            for k, v in colors.items():
                if re.match(k, label):
                    col = v
                    break

            if col:
                handle.set_color(col)
                handle.set_linestyle('-')
                handle.set_edgecolor('k')


    ax.legend(title=legend_titles[data_region], handles=handles, ncol=2)

    # Plot ratio
    h_data = h_data.integrate('dataset', data)
    h_mc = h_mc.integrate('dataset', mc)

    sumw_data, sumw2_data = h_data.values(overflow=overflow, sumw2=True)[()]
    sumw_mc = h_mc.values(overflow=overflow)[()]
    # Add the QCD contribution to the MC
    if data_region == 'sr_vbf':
        sumw_mc = sumw_mc + h_qcd.values * mcscale

    r = sumw_data / sumw_mc
    rerr = np.abs(poisson_interval(r, sumw2_data / sumw_mc**2) - r)

    r[np.isnan(r) | np.isinf(r)] = 0.
    rerr[np.isnan(rerr) | np.isinf(rerr)] = 0.

    hep.histplot(
        r,
        xedges,
        yerr=rerr,
        ax=rax,
        histtype='errorbar',
        **data_err_opts
    )

    xlabels = {
        'mjj': r'$M_{jj} \ (GeV)$',
        'ak4_eta0': r'Leading Jet $\eta$',
        'ak4_eta1': r'Trailing Jet $\eta$',
    }

    if distribution in xlabels.keys():
        ax.set_xlabel(xlabels[distribution])
        rax.set_xlabel(xlabels[distribution])
    
    rax.set_ylabel('Data / MC')
    rax.set_ylim(0.5,1.5)
    loc1 = MultipleLocator(0.2)
    loc2 = MultipleLocator(0.1)
    rax.yaxis.set_major_locator(loc1)
    rax.yaxis.set_minor_locator(loc2)

    rax.yaxis.set_ticks_position('both')

    sumw_denom, sumw2_denom = h_mc.values(overflow=overflow, sumw2=True)[()]

    unity = np.ones_like(sumw_denom)
    denom_unc = poisson_interval(unity, sumw2_denom / sumw_denom ** 2)
    opts = {"step": "post", "facecolor": (0, 0, 0, 0.3), "linewidth": 0}
    
    rax.fill_between(
        xedges,
        np.r_[denom_unc[0], denom_unc[0, -1]],
        np.r_[denom_unc[1], denom_unc[1, -1]],
        **opts
    )

    rax.grid(axis='y',which='both',linestyle='--')

    rax.axhline(1., xmin=0, xmax=1, color=(0,0,0,0.4), ls='--')

    fig.text(0., 1., 'UL Nanov8',
                fontsize=14,
                horizontalalignment='left',
                verticalalignment='bottom',
                transform=ax.transAxes
               )

    fig.text(1., 1., f'VBF, {lumi(year, mcscale):.1f} fb$^{{-1}}$ ({year})',
                fontsize=14,
                horizontalalignment='right',
                verticalalignment='bottom',
                transform=ax.transAxes
               )

    outdir = f'./output/{outtag}/{data_region}'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outpath = pjoin(outdir, f'{data_region}_data_mc_{distribution}_{year}.{fformat}')
    fig.savefig(outpath)
    plt.close(fig)

    print(f'File saved: {outpath}')

