#!/usr/bin/env python

import os
import re
import copy
from collections import defaultdict
from pprint import pprint

import matplotlib.ticker
import numpy as np
from coffea import hist
from coffea.util import load
from matplotlib import pyplot as plt

from bucoffea.plot.util import (acc_from_dir, lumi, merge_datasets,
                                merge_extensions, scale_xs_lumi)
from bucoffea.plot.stack_plot import Style, make_plot
from coffea.hist.plot import clopper_pearson_interval
pjoin = os.path.join

def histo_ratio(num, den):
    num_val, num_sw2 = num.values(sumw2=True)
    den_val, den_sw2 = den.values(sumw2=True)

colors ={
    'band' : '#fdd49e',
    'mc' : '#b30000'
}
def ratio_plot(acc, distribution='recoil', regions=['cr_2m_j','cr_1m_j','cr_1e_j','cr_2e_j','cr_g_j'], year=2017, outdir='./output'):
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    # Rebin
    s = Style()
    h = copy.deepcopy(acc[distribution])
    try:
        newax = s.rebin_axes[distribution]
        h = h.rebin(h.axis(newax.name), newax)
    except KeyError:
        pass

    h = merge_extensions(h, acc)

    scale_xs_lumi(h)

    h = merge_datasets(h)

    histograms = {}
    for region in regions:
        histograms[region] = copy.deepcopy(h).integrate(h.axis('region'), region)

    mc = {
        'cr_1m_j' : re.compile(f'(TTJets.*FXFX.*|Diboson.*|ST.*|QCD_HT.*|.*DYJetsToLL_M-50_HT_MLM.*|W.*HT.*).*{year}'),
        'cr_1e_j' : re.compile(f'(TTJets.*FXFX.*|Diboson.*|ST.*|QCD_HT.*|.*DYJetsToLL_M-50_HT_MLM.*|W.*HT.*).*{year}'),
        'cr_2m_j' : re.compile(f'(EW.*|TTJets.*FXFX.*|Diboson.*|ST.*|QCD_HT.*|.*DYJetsToLL_M-50_HT_MLM.*).*{year}'),
        'cr_2e_j' : re.compile(f'(EW.*|TTJets.*FXFX.*|Diboson.*|ST.*|QCD_HT.*|.*DYJetsToLL_M-50_HT_MLM.*).*{year}'),
        'cr_g_j' : re.compile(f'(GJets.*|QCD_HT.*|W.*HT.*).*{year}'),
    }
    data = {
        'cr_1m_j' : f'MET_{year}',
        'cr_2m_j' : f'MET_{year}',
        'cr_1e_j' : f'EGamma_{year}',
        'cr_2e_j' : f'EGamma_{year}',
        'cr_g_j' : f'EGamma_{year}',
    }
    name = {
        'cr_1m_j' : '1$\mu$',
        'cr_2m_j' : '2$\mu$',
        'cr_1e_j' : '1e',
        'cr_2e_j' : '2e',
        'cr_g_j' : '$\gamma$',
    }
    data_err_opts = {
        'linestyle':'none',
        'marker': '.',
        'markersize': 10.,
        'color':'k',
        'elinewidth': 1,
        # 'emarker': '_'
    }

    for i in range(len(regions)):
        for j in range(len(regions)):
            if i==j:
                continue
            fig, (ax, rax) = plt.subplots(2, 1, figsize=(7,7), gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
            h1 = histograms[regions[i]]
            h2 = histograms[regions[j]]

            print(data[regions[i]])
            h1_data = h1[data[regions[i]]].integrate('dataset')
            h1_mc = h1[mc[regions[i]]].integrate('dataset')
            h2_data = h2[data[regions[j]]].integrate('dataset')
            h2_mc = h2[mc[regions[j]]].integrate('dataset')
            # Ratio plot


            def ratio(num, den):
                num_sumw, num_sumw2 = num.values(sumw2=True, overflow='over')[()]
                den_sumw, den_sumw2 = den.values(sumw2=True, overflow='over')[()]
                rsumw_err = np.hypot(
                    np.sqrt(num_sumw2) / den_sumw,
                    num_sumw * np.sqrt(den_sumw2)/ den_sumw ** 2
                )
                rsumw = num_sumw/den_sumw

                return rsumw, rsumw_err


            data_err_opts['color'] = 'k'
            rsumw_data, rsumw_err_data = ratio(h1_data, h2_data)
            ax.errorbar(x=h1_data.axis(distribution).centers(overflow='over'), y=rsumw_data, yerr=rsumw_err_data,label='Data', **data_err_opts)

            # data_err_opts['color'] = 'r'
            rsumw_mc, rsumw_err_mc = ratio(h1_mc, h2_mc)
            edges = h1_mc.axis(distribution).edges(overflow='over')
            ax.step(
                x=edges,
                y=np.r_[rsumw_mc[0], rsumw_mc],
                color=colors['mc'],
                label='MC'
                )

            y1 = np.r_[rsumw_mc - rsumw_err_mc, rsumw_mc[-1] - rsumw_err_mc[-1]]
            y2 =np.r_[rsumw_mc + rsumw_err_mc, rsumw_mc[-1] + rsumw_err_mc[-1]]

            ax.fill_between(edges,
                                    y1 = y1,
                                    y2 = y2,
                                    zorder=-1,
                                    color=colors['band'],
                                    step='post',
                                    label='MC stat. unc'
                                    )

            rrsumw = rsumw_data / rsumw_mc
            rrsumw_err = rsumw_err_data / rsumw_mc
            rax.errorbar(
                         x=h1_data.axis(distribution).centers(overflow='over'),
                         y = rrsumw,
                         yerr =rrsumw_err,
                         **data_err_opts
                         )

            rax.set_ylim(0.75,1.25)

            plt.plot([min(edges), max(edges)],[1,1],color=colors['mc'])


            y1 = np.r_[(rsumw_mc - rsumw_err_mc)/rsumw_mc, (rsumw_mc[-1] - rsumw_err_mc[-1])/rsumw_mc[-1]]
            y2 =np.r_[(rsumw_mc + rsumw_err_mc)/rsumw_mc, (rsumw_mc[-1] + rsumw_err_mc[-1])/rsumw_mc[-1]]

            rax.fill_between(edges,
                                    y1 = y1,
                                    y2 = y2,
                                    zorder=-1,
                                    color=colors['band'],
                                    step='post'
                                    )
            ax.legend(title=f'{name[regions[i]]} over {name[regions[j]]}')
            fig.text(1., 1., f'{lumi(year)} fb$^{{-1}}$ ({year})',
                fontsize=14,
                horizontalalignment='right',
                verticalalignment='bottom',
                transform=ax.transAxes
               )
            fig.text(0., 1., '$\\bf{CMS}$ internal',
                fontsize=14,
                horizontalalignment='left',
                verticalalignment='bottom',
                transform=ax.transAxes
               )

            rax.set_xlabel('Recoil (GeV)',fontsize=14)
            rax.set_ylabel('Data / MC',fontsize=14)
            ax.set_xlabel('Recoil (GeV)',fontsize=14)
            ax.set_ylabel(f'Region ratio: {name[regions[i]]} / {name[regions[j]]} (GeV)',fontsize=14)

            loc1 = matplotlib.ticker.MultipleLocator(base=0.2)
            loc2 = matplotlib.ticker.MultipleLocator(base=0.1)
            rax.yaxis.set_major_locator(loc1)
            rax.yaxis.set_minor_locator(loc2)
            rax.grid(axis='y',which='minor',linestyle='--')
            rax.grid(axis='y',which='major',linestyle='--')
            # Save and close
            fig.savefig(
                pjoin(
                    outdir, f'ratio_{distribution}_{regions[i]}_over_{regions[j]}_{year}.pdf'
                        )
                )
            plt.close(fig)


def main():
    infile=os.path.abspath('../after_new_sfs/input/2019-09-05_all_new_sf_neweletrig/')
    acc = acc_from_dir(infile)
    for year in [2017,2018]:
        ratio_plot(acc, year=year,outdir=f'./output/{os.path.basename(infile)}')

if __name__ == "__main__":
    main()














# for distribution in ['ak4_pt0_eta0']:
#     make_plot(acc, region=region,distribution=distribution, year=year, data=data, mc=mc, integrate=('jeteta',0,1))
#     make_plot(acc, region=region,distribution=distribution, year=year, data=data, mc=mc, integrate=('jeteta',1,3))
#     make_plot(acc, region=region,distribution=distribution, year=year, data=data, mc=mc, integrate=('jetpt',100,200))
#     make_plot(acc, region=region,distribution=distribution, year=year, data=data, mc=mc, integrate=('jetpt',200,500))






# for year in [2017]:
#     data = re.compile(f'MET_{year}')
#     mc = re.compile(f'(TTJets.*|ZZ.*|ST.*|QCD_HT-mg.*|WW.*|WZ.*|W.*HT.*).*{year}')
#     region='cr_1m_j'
#     for distribution in ['ak4_eta0','ak4_btag','ak4_pt','ak4_eta','ak4_ptraw0','muon_pt','met','recoil','ak4_pt0']:
#         make_plot(copy.deepcopy(acc), region=region,distribution=distribution, year=year, data=data, mc=mc)
# for year in [2017]:
#     data = re.compile(f'EGamma.*{year}')
#     mc = re.compile(f'(GJets.*HT|QCD_HT-mg.*).*{year}')
#     region='cr_g_j'
#     for distribution in ['recoil','ak4_pt0','drphotonjet']:
#         make_plot(copy.deepcopy(acc), region=region,distribution=distribution, year=year, data=data, mc=mc)
