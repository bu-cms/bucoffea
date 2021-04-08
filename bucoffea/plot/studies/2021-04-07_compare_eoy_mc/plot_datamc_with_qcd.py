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
    'QCD.*' : "QCD",
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
    '.*TT.*' : '#6a51a3',
    '.*ST.*' : '#9e9ac8',
    'ZJetsToNuNu.*' : '#31a354',
    'WJets.*' : '#feb24c',
}

def plot_datamc_with_qcd(acc, outtag, year, region='sr_vbf', distribution='mjj', mcscale=1):
    '''Plot data/MC comparison with the QCD template included.'''
    acc.load(distribution)
    h = acc[distribution]

    h = merge_extensions(h, acc, reweight_pu=False)
    scale_xs_lumi(h)
    h = merge_datasets(h)

    if distribution == 'mjj':
        mjj_ax = hist.Bin('mjj', r'$M_{jj} \ (GeV)$', [200., 400., 600., 900., 1200., 1500., 2000., 2750., 3500., 5000.])
        h = h.rebin('mjj', mjj_ax)

    h.axis('dataset').sorting = 'integral'

    h = h.integrate('region', region)
    data = f'MET_{year}'
    mc = re.compile(f'(ZJetsToNuNu.*|EW.*|Top_FXFX.*|Diboson.*|.*DYJetsToLL_M-50_HT_MLM.*|.*WJetsToLNu.*HT.*).*{year}')
    
    h.scale({
        ds : (mcscale  if mc.match(ds) else 1) for ds in map(str,h.axis("dataset").identifiers())
    }, axis='dataset')

    # Get the QCD template (estimation from HF)
    qcdfilepath = bucoffea_path('data/templates/qcd_estimate_sr.root')    
    h_qcd = uproot.open(qcdfilepath)[f'sr_template_{distribution}']

    data_err_opts = {
        'linestyle':'none',
        'marker': '.',
        'markersize': 10.,
        'color':'k',
        'elinewidth': 1,
    }

    qcd_err_opts = {
        'color':'crimson',
    }

    fig, ax, rax = fig_ratio()
    hist.plot1d(h[data], ax=ax, overlay='dataset', binwnorm=1, error_opts=data_err_opts)
    hist.plot1d(h[mc], ax=ax, overlay='dataset', stack=True, clear=False, binwnorm=1)

    hep.histplot(h_qcd.values * mcscale, h_qcd.edges, ax=ax, label='HF QCD Estimate', binwnorm=1, **qcd_err_opts)

    ax.set_yscale('log')
    ax.set_ylim(1e-4,1e4)
    ax.set_ylabel('Events / GeV')

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
    h_data = h.integrate('dataset', data)
    h_mc = h.integrate('dataset', mc)

    sumw_data, sumw2_data = h_data.values(sumw2=True)[()]
    # Add the QCD contribution to the MC
    sumw_mc = h_mc.values()[()] + h_qcd.values * mcscale

    r = sumw_data / sumw_mc
    rerr = np.abs(poisson_interval(r, sumw2_data / sumw_mc**2) - r)

    hep.histplot(
        r,
        h_qcd.edges,
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

    rax.set_xlabel(xlabels[distribution])
    rax.set_ylabel('Data / MC')
    rax.set_ylim(0,2)
    rax.grid(True)
    loc1 = MultipleLocator(0.5)
    loc2 = MultipleLocator(0.25)
    rax.yaxis.set_major_locator(loc1)
    rax.yaxis.set_minor_locator(loc2)

    rax.yaxis.set_ticks_position('both')

    fig.text(1., 1., f'VBF, {lumi(year, mcscale):.1f} fb$^{{-1}}$ ({year})',
                fontsize=14,
                horizontalalignment='right',
                verticalalignment='bottom',
                transform=ax.transAxes
               )

    outdir = f'./output/{outtag}'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outpath = pjoin(outdir, f'datamc_{distribution}.pdf')
    fig.savefig(outpath)
    plt.close(fig)

    print(f'File saved: {outpath}')

def commandline():
    parser = argparse.ArgumentParser(prog='Plotter.')
    parser.add_argument('inpath', type=str, help='Input folder to use.')
    parser.add_argument('--region', type=str, default='sr_vbf', help='Region to plot.')
    parser.add_argument('--distribution', type=str, default='.*', help='Regex specifying the distributions to plot.')
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

    for distribution in distributions:
        if not re.match(args.distribution, distribution):
            continue
        plot_datamc_with_qcd(acc, outtag, year=2017, mcscale=mcscale, distribution=distribution)

if __name__ == "__main__":
    main()