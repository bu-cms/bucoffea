#!/usr/bin/env python
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
def plot_lhe_v_pt(acc, tag, regex, outputrootfile):
    outdir = './output/'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    fig, (ax, rax) = plt.subplots(2, 1, figsize=(7,7), gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
    # new_ax = hist.Bin('vpt','LHE V $p_{T}$ (GeV)',list(range(100,500,50)) + list(range(500,1000,100)) + list(range(1000,2000,250)))
    new_ax = hist.Bin('vpt','LHE V $p_{T}$ (GeV)',list(range(80,800,40))+list(range(800,2000,100)))

    for dist in ['gen_vpt']:
        h = copy.deepcopy(acc[dist])
        h = h.rebin(h.axis('vpt'), new_ax)
        h = merge_extensions(h, acc, reweight_pu=False)
        scale_xs_lumi(h)
        h = merge_datasets(h)
        h = h[re.compile(regex)]
        h = h.integrate('weight_type','nominal')
        h = h.integrate('weight_index',slice(-0.5,0.5))
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

        # try:
        #     f = uproot.create(f'gen_v_pt_qcd_sf.root')
        # except OSError:
        
        outputrootfile[tag] = (sf_y,sf_x)


def pdfwgt_sf(vpt):
    return 1/(1.157 + 2.291e-4 * vpt + 6.0612e-7 * vpt**2)

# def sf_comparison(acc, tag, regex):
#     outdir = './output/'
#     if not os.path.exists(outdir):
#         os.makedirs(outdir)

#     fig = plt.gcf()
#     # new_ax = hist.Bin('vpt','LHE V $p_{T}$ (GeV)',list(range(100,500,50)) + list(range(500,1000,100)) + list(range(1000,2000,250)))
#     new_ax = hist.Bin('vpt','LHE V $p_{T}$ (GeV)',list(range(80,2000,40)))

#     for dist in ['gen_vpt']:
#         fig.clf()
#         h = copy.deepcopy(acc[dist])
#         h = h.rebin(h.axis('vpt'), new_ax)
#         h = merge_extensions(h, acc, reweight_pu=False)
#         scale_xs_lumi(h)
#         h = merge_datasets(h)
#         h = h[re.compile(regex)]
#         h = h.integrate('weight_type','nominal')
#         h = h.integrate('weight_index',slice(-0.5,0.5))

#         lo = h[re.compile('.*HT.*')].integrate('dataset')
#         nlo = h[re.compile('.*LHE.*')].integrate('dataset')

#         x_sf = lo.axis('vpt').identifiers()
#         y_sf = nlo.values () / lo.values()

#         old = get_old_kfac(tag)

#         rax.plot(0.5*(old.bins[:,0]+old.bins[:,1]), old.values,'or-', label='2016 QCD k fac')
#         ax.set_yscale('log')
#         ax.set_ylim(1e-3,1e6)
#         rax.set_ylim(0,2)
#         rax.legend()

#         fig.savefig(pjoin(outdir,f'{tag}_{dist}.pdf'))
def main():
    acc = acc_from_dir("./input/das_lhevpt_v1")
    outputrootfile = uproot.recreate(f'gen_v_pt_qcd_sf.root')
    plot_lhe_v_pt(acc, tag='wjet', regex='W.*',outputrootfile=outputrootfile)
    plot_lhe_v_pt(acc, tag='dy', regex='DY.*',outputrootfile=outputrootfile)




if __name__ == "__main__":
    main()
