#!/usr/bin/env python
from pprint import pprint
import copy
import os
import re
import numpy as np
from matplotlib import pyplot as plt

from coffea import hist

from bucoffea.plot.util import (acc_from_dir, merge_datasets, merge_extensions,
                                scale_xs_lumi)
from bucoffea.helpers.paths import bucoffea_path
import uproot
pjoin = os.path.join

data_err_opts = {
        'linestyle':'none',
        'marker': '.',
        'markersize': 10.,
        'color':'k',
        'elinewidth': 1,
        'emarker': '_'
    }

def get_old_kfac(tag):
    if tag.startswith('w'):
        f = uproot.open(bucoffea_path('data/sf/theory/merged_kfactors_wjets.root'))
    elif tag.startswith('dy'):
        f = uproot.open(bucoffea_path('data/sf/theory/merged_kfactors_wjets.root'))
    return f['kfactor_monojet_qcd']

def sf_1d(acc, tag, regex, outputrootfile):
    outdir = './output/'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    fig, (ax, rax) = plt.subplots(2, 1, figsize=(7,7), gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
    # new_ax = hist.Bin('vpt','LHE V $p_{T}$ (GeV)',list(range(100,500,50)) + list(range(500,1000,100)) + list(range(1000,2000,250)))
    new_ax = hist.Bin('vpt','V $p_{T}$ (GeV)',list(range(100,800,100))+list(range(800,1200,200))+list(range(1200,2800,800)))

    for selection in ['inclusive','monojet','vbf']:
        dist = f'gen_vpt_{selection}'
        h = copy.deepcopy(acc[dist])

        h = h.rebin(h.axis('vpt'), new_ax)

        if selection == 'monojet':
            h = h.integrate(h.axis("jpt"))
        if selection == 'vbf':
            h = h.integrate(h.axis("jpt"))
            h = h.integrate(h.axis("mjj"))
        h = merge_extensions(h, acc, reweight_pu=False)
        scale_xs_lumi(h)
        h = merge_datasets(h)
        h = h[re.compile(regex)]
        hist.plot1d(
            h,
            overlay='dataset',
            overflow='all',
            binwnorm=True,
            ax=ax)
        lo = h[re.compile('.*HT.*')].integrate('dataset')
        nlo = h[re.compile('.*LHE.*')].integrate('dataset')

        hist.plotratio(nlo, lo,
            ax=rax,
            denom_fill_opts={},
            guide_opts={},
            unc='num',
            overflow='all',
            error_opts=data_err_opts,
            label='2017 NLO/LO ratio'
            )
        old = get_old_kfac(tag)
        old_x = 0.5*(old.bins[:,0]+old.bins[:,1])
        rax.plot(old_x, old.values,'ob-', label='2016 QCD k fac')
        rax.plot(old_x, old.values * pdfwgt_sf(old_x),'or-', label='2016 x ad-hoc DY pdfwgt SF')
        ax.set_yscale('log')
        ax.set_ylim(1e-3,1e6)
        rax.set_ylim(0,2)
        rax.legend()


        fig.savefig(pjoin(outdir,f'{tag}_{dist}.pdf'))

        sf_x = lo.axis('vpt').edges()
        sf_y = nlo.values()[()] / lo.values()[()]

        outputrootfile[f'{tag}_{selection}'] = (sf_y,sf_x)


def sf_2d(acc, tag, regex, outputrootfile):
    outdir = './output/2d/'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # fig, (ax, rax) = plt.subplots(2, 1, figsize=(7,7), gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
    fig = plt.gcf()
    ax = plt.gca()
    # new_ax = hist.Bin('vpt','LHE V $p_{T}$ (GeV)',list(range(100,500,50)) + list(range(500,1000,100)) + list(range(1000,2000,250)))
    vpt_ax = hist.Bin('vpt','V $p_{T}$ (GeV)',list(range(250,1000,150))+list(range(1000,1200,200))+list(range(1200,200,800)))
    mjj_ax = hist.Bin('mjj','M(jj) (GeV)',[200]+list(range(500,2500,500)))

    for selection in ['vbf']:
        dist = f'gen_vpt_{selection}'
        h = copy.deepcopy(acc[dist])

        h = h.rebin(h.axis('vpt'), vpt_ax)
        h = h.rebin(h.axis('mjj'), mjj_ax)
        h = h.integrate(h.axis("jpt"))

        h = merge_extensions(h, acc, reweight_pu=False)
        scale_xs_lumi(h)
        h = merge_datasets(h)
        h = h[re.compile(regex)]

        lo = h[re.compile('.*HT.*')].integrate('dataset')
        nlo = h[re.compile('.*LHE.*')].integrate('dataset')

        sumw_lo, sumw2_lo = lo.values(overflow='over', sumw2=True)[()]
        sumw_nlo, sumw2_nlo = nlo.values(overflow='over', sumw2=True)[()]

        print(sumw_nlo)
        sf = sumw_nlo / sumw_lo
        dsf = np.hypot(
            np.sqrt(sumw2_nlo) / sumw_lo,
            sumw_nlo * np.sqrt(sumw2_lo) / (sumw_lo**2)
        )
        xaxis = lo.axes()[0]
        yaxis = lo.axes()[1]

        im = ax.pcolormesh(xaxis.edges(overflow='over'), yaxis.edges(overflow='over'), sf.T)

        x_centers = xaxis.centers(overflow='over')
        y_centers = yaxis.centers(overflow='over')
        for ix in range(len(x_centers)):
            for iy in range(len(y_centers)):
                textcol = 'white' if sf.T[iy, ix] <0.9 else 'black'
                ax.text(
                        x_centers[ix], 
                        y_centers[iy], 
                        f'  {sf.T[iy, ix]:.2f} \n$\\pm$ {dsf.T[iy, ix]:.2f}',
                        ha='center',
                        va='center',
                        color=textcol,
                        fontsize=6
                        )
        # hist.plotratio(nlo, lo,
        #     ax=rax,
        #     denom_fill_opts={},
        #     guide_opts={},
        #     unc='num',
        #     overflow='all',
        #     error_opts=data_err_opts,
        #     label='2017 NLO/LO ratio'
        #     )
        # old = get_old_kfac(tag)
        # old_x = 0.5*(old.bins[:,0]+old.bins[:,1])
        # rax.plot(old_x, old.values,'ob-', label='2016 QCD k fac')
        # rax.plot(old_x, old.values * pdfwgt_sf(old_x),'or-', label='2016 x ad-hoc DY pdfwgt SF')
        # ax.set_yscale('log')
        # ax.set_ylim(1e-3,1e6)
        # rax.set_ylim(0,2)
        # rax.legend()

        ax.set_ylabel('$p_{T}(V)$ (GeV)')
        ax.set_xlabel('M(jj) (GeV)')
        cb = fig.colorbar(im)
        cb.set_label('LO $\\rightarrow$ NLO SF')
        im.set_clim(0.75,1.1)
        fig.savefig(pjoin(outdir,f'2d_{tag}_{dist}.pdf'))

        # sf_x = lo.axis('vpt').edges()
        # sf_y = nlo.values()[()] / lo.values()[()]

        # outputrootfile[f'{tag}_{selection}'] = (sf_y,sf_x)

def pdfwgt_sf(vpt):
    return 1/(1.157 + 2.291e-4 * vpt + 6.0612e-7 * vpt**2)

def main():
    acc = acc_from_dir("./input/2019-10-07_das_lhevpt_v6")
    # outputrootfile = uproot.recreate(f'2017_gen_v_pt_stat1_qcd_sf.root')
    # sf_1d(acc, tag='wjet', regex='W.*',outputrootfile=outputrootfile)
    # sf_1d(acc, tag='dy', regex='.*DY.*',outputrootfile=outputrootfile)

    sf_2d(acc, tag='wjet', regex='W.*',outputrootfile=None)




if __name__ == "__main__":
    main()
