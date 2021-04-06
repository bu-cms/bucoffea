#!/usr/bin/env python
import os
import re
import sys
import uproot
import numpy as np

from matplotlib import pyplot as plt
from coffea import hist
from coffea.hist import poisson_interval
from bucoffea.plot.util import merge_datasets, merge_extensions, scale_xs_lumi
from bucoffea.plot.style import matplotlib_rc
from klepto.archives import dir_archive

pjoin = os.path.join

matplotlib_rc()

def derive_weights(acc, outtag, rootf, outdir, year=2017, distribution='ak4_eta0', region='sr_vbf', mcscale=1):
    acc.load(distribution)
    h = acc[distribution]

    h = merge_extensions(h, acc, reweight_pu=False)
    scale_xs_lumi(h)
    h = merge_datasets(h)

    h = h.integrate('region', region)
    
    data = f'MET_{year}'
    mc = re.compile(f'(ZJetsToNuNu.*|EW.*|Top_FXFX.*|Diboson.*|.*DYJetsToLL_M-50_HT_MLM.*|.*WJetsToLNu.*HT.*).*{year}')
    
    h.scale({
        ds : (mcscale  if mc.match(ds) else 1) for ds in map(str,h.axis("dataset").identifiers())
    }, axis='dataset')

    h_data = h.integrate('dataset', data)
    h_mc = h.integrate('dataset', mc)

    data_err_opts = {
        'linestyle':'none',
        'marker': '.',
        'markersize': 10.,
        'color':'k',
    }

    fig, ax = plt.subplots()
    data_sumw, data_sumw2 = h_data.values(sumw2=True)[()]
    mc_sumw = h_mc.values()[()]
    rsumw = h_data.values()[()] / h_mc.values()[()]
    eta_edges = h_data.axis('jeteta').edges()
    eta_ax = h_data.axis('jeteta').centers()

    rsumw_err = np.abs(poisson_interval(rsumw, data_sumw2 / mc_sumw**2) - rsumw)

    jet_in_endcap = (np.abs(eta_ax) > 2.4) & (np.abs(eta_ax) < 3.5)
    ratio = np.where(
        jet_in_endcap,
        rsumw,
        1.0
    )
    
    ratio_err = np.where(
        jet_in_endcap,
        rsumw_err,
        0.
    )

    ax.errorbar(eta_ax, ratio, yerr=ratio_err, **data_err_opts)

    ax.set_ylabel('Data / MC')
    new_xlabels = {
        'ak4_eta0': r'Leading Jet $\eta$',
        'ak4_eta1': r'Trailing Jet $\eta$',
    }
    ax.set_xlabel(new_xlabels[distribution])
    ax.set_ylim(0,2)
    ax.grid(True)

    ax.text(0.,1.,'VBF SR 2017',
        fontsize=14,
        ha='left',
        va='bottom',
        transform=ax.transAxes
    )

    ax.text(1.,1.,'UL Data & ReReco MC',
        fontsize=14,
        ha='right',
        va='bottom',
        transform=ax.transAxes
    )

    outpath = pjoin(outdir, f'{distribution}.pdf')
    fig.savefig(outpath)
    plt.close(fig)
    print(f'File saved: {outpath}')

    rootf[f'mc_weights_{distribution}'] = (ratio, eta_edges)

def main():
    inpath = sys.argv[1]
    acc = dir_archive(inpath)
    acc.load('sumw')
    acc.load('sumw2')

    outtag = re.findall('merged_.*', inpath)[0].replace('/','')

    outdir = f'./output/{outtag}'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    rootpath = pjoin(outdir, 'jet_eta_weights.root')
    rootf = uproot.recreate(rootpath)

    for distribution in ['ak4_eta0', 'ak4_eta1']:
        derive_weights(acc, outtag, rootf, outdir, distribution=distribution, mcscale=0.2)

if __name__ == '__main__':
    main()