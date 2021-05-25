#!/usr/bin/env python

import os
import sys
import re
import numpy as np
import mplhep as hep

from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
from coffea import hist
from coffea.hist import poisson_interval
from bucoffea.plot.util import merge_datasets, merge_extensions, scale_xs_lumi, fig_ratio
from bucoffea.plot.plotter import legend_labels, colors
from klepto.archives import dir_archive
from pprint import pprint

pjoin = os.path.join

def pretty_eta_label(etaslice):
    return f'${etaslice.start:.2f} < |\\eta_{{j0}}| < {etaslice.stop:.2f}$' 

def preprocess(h, acc, etaslice):
    h = merge_extensions(h, acc)
    scale_xs_lumi(h)
    h = merge_datasets(h)

    # Rebin dhitkpf
    new_bins = [ibin.lo for ibin in h.identifiers('dphi') if ibin.lo < 2] + [3.5]
    
    new_ax = hist.Bin('dphi', r'$\Delta\phi_{TK,PF}$', new_bins)
    h = h.rebin('dphi', new_ax)

    # Integrate out the eta slice
    h = h.integrate('jeteta', etaslice)
    return h

def get_qcd_estimation_for_etaslice(h, outtag, year, etaslice=slice(3, 3.25), fformat='pdf'):
    '''Get QCD estimation for the given leading jet eta slice.'''
    # QCD CR
    region = 'cr_vbf_qcd'
    h = h.integrate('region', region)

    data = f'MET_{year}'
    mc = re.compile(f'(ZJetsToNuNu.*|EW.*|Top_FXFX.*|Diboson.*|DYJetsToLL_M-50_HT_MLM.*|WJetsToLNu.*HT.*).*{year}')

    h_data = h.integrate('dataset', data)
    h_mc = h.integrate('dataset', mc)

    # Data - MC gives the estimation in SR (TF pre-applied)
    h_mc.scale(-1)
    h_data.add(h_mc)

    fig, ax = plt.subplots()
    hist.plot1d(h_data, ax=ax)

    ax.set_yscale('log')
    ax.set_ylim(1e-4,1e4)
    ax.set_ylabel('HF Estimation')

    ax.get_legend().remove()
    ax.yaxis.set_ticks_position('both')

    ax.text(0.,1.,pretty_eta_label(etaslice),
        fontsize=14,
        ha='left',
        va='bottom',
        transform=ax.transAxes
    )

    ax.text(1.,1.,year,
        fontsize=14,
        ha='right',
        va='bottom',
        transform=ax.transAxes
    )

    outdir = f'./output/{outtag}/dphitkpf'
    try:
        os.makedirs(outdir)
    except FileExistsError:
        pass

    etatag = f'{str(etaslice.start).replace(".", "_")}_{str(etaslice.stop).replace(".","_")}'

    outpath = pjoin(outdir, f'hf_estimation_eta_{etatag}_{year}.{fformat}')
    fig.savefig(outpath)
    plt.close(fig)

    print(f'File saved: {outpath}')

    # Return the histogram containing the QCD template
    return h_data

def plot_dphitkpf(acc, outtag, year, region='sr_vbf', distribution='dphitkpf_ak4_eta0', etaslice=slice(3, 3.25), fformat='pdf'):
    '''Plot dphitkpf distribution in data and MC in a stack plot, for the given eta slice for the leading jet.'''
    acc.load(distribution)
    h = preprocess(acc[distribution], acc, etaslice)

    # Get the QCD template
    h_qcd = get_qcd_estimation_for_etaslice(h, outtag, year, etaslice=etaslice, fformat=fformat)

    h = h.integrate('region', region)

    data = f'MET_{year}'
    mc = re.compile(f'(ZJetsToNuNu.*|EW.*|Top_FXFX.*|Diboson.*|DYJetsToLL_M-50_HT_MLM.*|WJetsToLNu.*HT.*).*{year}')
    datasets = list(map(str, h[mc].identifiers('dataset')))
    
    data_err_opts = {
        'linestyle':'none',
        'marker': '.',
        'markersize': 10.,
        'color':'k',
        'elinewidth': 1,
    }
    
    h_data = h.integrate('dataset', data)

    # Stack plot for MC
    plot_info = {
        'label' : [],
        'sumw' : [],
    }

    fig, ax, rax = fig_ratio()
    for dataset in datasets:
        sumw = h[mc].integrate('dataset', dataset).values()[()]
        if (sumw == 0.).all():
            continue
        plot_info['label'].append(dataset)
        plot_info['sumw'].append(sumw)

    # Add the QCD contribution!
    plot_info['label'].insert(6, 'HF Noise Estimate')
    sumw_qcd = h_qcd.values()[()]
    sumw_qcd[sumw_qcd < 0] = 0.
    plot_info['sumw'].insert(6, sumw_qcd)

    xedges = h_data.axis('dphi').edges()

    hist.plot1d(h_data, ax=ax, binwnorm=1, error_opts=data_err_opts)

    hep.histplot(plot_info['sumw'], xedges, 
        ax=ax,
        label=plot_info['label'], 
        histtype='fill',
        binwnorm=1,
        stack=True
        )

    signal = re.compile(f'VBF_HToInvisible.*withDipoleRecoil.*{year}')

    signal_line_opts = {
        'linestyle': '-',
        'color': 'crimson',
    }

    hist.plot1d(
        h[signal],
        ax=ax,
        overlay='dataset',
        line_opts=signal_line_opts,
        binwnorm=1,
        clear=False
    )
    
    ax.set_yscale('log')
    ax.set_ylim(1e-2,1e6)
    ax.yaxis.set_ticks_position('both')

    handles, labels = ax.get_legend_handles_labels()
    for handle, label in zip(handles, labels):
        if label == 'None':
            handle.set_label('Data')
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

    ax.legend(handles=handles, ncol=2)

    ax.text(0.,1.,pretty_eta_label(etaslice),
        fontsize=14,
        ha='left',
        va='bottom',
        transform=ax.transAxes
    )

    ax.text(1.,1.,year,
        fontsize=14,
        ha='right',
        va='bottom',
        transform=ax.transAxes
    )

    # Plot ratio
    h_mc = h[mc].integrate('dataset', mc)

    sumw_data, sumw2_data = h_data.values(sumw2=True)[()]
    sumw_mc = h_mc.values()[()]
    # Add the QCD contribution to the MC
    sumw_mc = sumw_mc + sumw_qcd

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

    rax.set_ylabel('Data / MC')
    rax.set_ylim(0.5,1.5)
    loc1 = MultipleLocator(0.2)
    loc2 = MultipleLocator(0.1)
    rax.yaxis.set_major_locator(loc1)
    rax.yaxis.set_minor_locator(loc2)

    rax.yaxis.set_ticks_position('both')

    sumw_denom, sumw2_denom = h_mc.values(sumw2=True)[()]

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
    
    outdir = f'./output/{outtag}/dphitkpf'
    try:
        os.makedirs(outdir)
    except FileExistsError:
        pass

    etatag = f'{str(etaslice.start).replace(".", "_")}_{str(etaslice.stop).replace(".","_")}'
    
    outpath = pjoin(outdir, f'data_mc_eta_{etatag}_{year}.{fformat}')
    fig.savefig(outpath)
    plt.close(fig)

    print(f'File saved: {outpath}')

def main():
    inpath = sys.argv[1]
    acc = dir_archive(inpath)
    acc.load('sumw')
    acc.load('sumw_pileup')
    acc.load('nevents')

    outtag = re.findall('merged_.*', inpath)[0].replace('/','')

    etaslices = [
        slice(0, 2.5),
        slice(2.5, 3),
        slice(3, 3.25),
        slice(3.25, 5),
    ]
    
    for year in [2017, 2018]:
        for etaslice in etaslices:
            for fformat in ['pdf', 'png']:
                plot_dphitkpf(acc, outtag, year=year, etaslice=etaslice, fformat=fformat)

if __name__ == '__main__':
    main()