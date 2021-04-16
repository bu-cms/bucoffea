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
    'MET|Single(Electron|Photon|Muon)|EGamma.*' : "Data"
}

colors = {
    'DY.*' : '#ffffcc',
    'EWKW.*' : '#c6dbef',
    'EWKZ.*ZToLL.*' : '#d5bae2',
    'EWKZ.*ZToNuNu.*' : '#c4cae2',
    '.*Diboson.*' : '#4292c6',
    'Top.*' : '#6a51a3',
    '.*QCD.*' : '#08306b',
    '.*TT.*' : '#6a51a3',
    '.*ST.*' : '#9e9ac8',
    'ZJetsToNuNu.*' : '#31a354',
    'WJets.*' : '#feb24c',
}

def plot_datamc_with_qcd(acc, outtag, year, region='sr_vbf', distribution='mjj', mcscale=1):
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

    if distribution == 'mjj':
        mjj_ax = hist.Bin('mjj', r'$M_{jj} \ (GeV)$', [200., 400., 600., 900., 1200., 1500., 2000., 2750., 3500.])
        h = h.rebin('mjj', mjj_ax)

    h.axis('dataset').sorting = 'integral'

    mc_region_suffix = ''
    if region == 'sr_vbf':
        mc_region_suffix = '_no_veto_all'

    regions = {
        'data': f'{region}',
        'mc': f'{region}{mc_region_suffix}',
    }

    data = f'MET_{year}'
    mc = re.compile(f'(ZJetsToNuNu.*|EW.*|Top_FXFX.*|Diboson.*|DYJetsToLL_M-50_HT_MLM.*|WJetsToLNu.*HT.*).*{year}')
    
    h.scale({
        ds : (mcscale  if mc.match(ds) else 1) for ds in map(str,h.axis("dataset").identifiers())
    }, axis='dataset')

    h_data = h.integrate('region', regions['data'])
    h_mc = h.integrate('region', regions['mc'])

    # Get the QCD template (estimation from HF)
    # qcdfilepath = bucoffea_path('data/templates/qcd_estimate_sr.root')    
    qcdfilepath = 'output/merged_2021-04-16_vbfhinv_ULv8_05Feb21/qcd_estimate/hf_qcd_estimate.root'    
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

    for dataset in datasets:
        sumw = h_mc.integrate('dataset', dataset).values(overflow=overflow)[()]

        plot_info['sumw'].append(sumw)

    # Add the QCD contribution
    plot_info['label'].insert(6, 'QCD Estimation')
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

    ax.set_yscale('log')
    if distribution == 'mjj':
        ax.set_ylim(1e-3,1e5)
        ax.set_ylabel('Events / GeV')
    elif 'ak4_eta' in distribution:
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

    ax.legend(title='VBF Signal Region', handles=handles, ncol=2)

    # Plot ratio
    h_data = h_data.integrate('dataset', data)
    h_mc = h_mc.integrate('dataset', mc)

    sumw_data, sumw2_data = h_data.values(overflow=overflow, sumw2=True)[()]
    # Add the QCD contribution to the MC
    sumw_mc = h_mc.values(overflow=overflow)[()] + h_qcd.values * mcscale

    r = sumw_data / sumw_mc
    rerr = np.abs(poisson_interval(r, sumw2_data / sumw_mc**2) - r)

    r[np.isnan(r)] = 1.
    rerr[np.isnan(rerr)] = 0.

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

    outdir = f'./output/{outtag}'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outpath = pjoin(outdir, f'{region}_data_mc_{distribution}_{year}.pdf')
    fig.savefig(outpath)
    plt.close(fig)

    print(f'File saved: {outpath}')

def commandline():
    parser = argparse.ArgumentParser(prog='Plotter.')
    parser.add_argument('inpath', type=str, help='Input folder to use.')
    parser.add_argument('--region', type=str, default='sr_vbf', help='Region to plot.')
    parser.add_argument('--distribution', type=str, default='.*', help='Regex specifying the distributions to plot.')
    parser.add_argument('--years', nargs='*', default=[2017,2018], help='Years to run on.')
    parser.add_argument('--one_fifth_unblind', action='store_true', help='1/5th unblinded data.')
    args = parser.parse_args()
    return args

def main():
    args = commandline()
    
    acc = dir_archive(args.inpath)
    acc.load('sumw')
    acc.load('sumw2')

    outtag = re.findall('merged_.*', args.inpath)[0].replace('/','')

    if args.one_fifth_unblind:
        mcscale = 0.2
    else:
        mcscale = 1

    distributions = ['mjj','ak4_eta0','ak4_eta1']

    for year in args.years:
        for distribution in distributions:
            if not re.match(args.distribution, distribution):
                continue
            plot_datamc_with_qcd(acc, outtag, year=year, mcscale=mcscale, distribution=distribution)

if __name__ == "__main__":
    main()