#!/usr/bin/env python
import os
import re
from pprint import pprint

import matplotlib.ticker
import numpy as np
import uproot
from coffea.util import load
from matplotlib import pyplot as plt

from bucoffea.plot.stack_plot import Style, make_plot
from bucoffea.plot.util import (acc_from_dir, lumi, merge_datasets,
                                merge_extensions, scale_xs_lumi)

regions = {
'ch1' : 'sr_j',
'ch2' : 'cr_2m_j',
'ch3' : 'cr_1m_j',
'ch4' : 'cr_g_j',
'ch5' : 'cr_2e_j',
'ch6' : 'cr_1e_j'
}
tmp = {}
for k,v in regions.items():
    tmp[v] = k
regions.update(tmp)

axes = {
    'recoil' : 'recoil',
    'ak4_pt0' : 'jetpt',
    'ak4_eta0' : 'jeteta',
    'ak4_eta' : 'jeteta',
    'met' : 'met',
    'muon_pt0' : 'pt',
    'muon_pt1' : 'pt',
    'muon_phi0' : 'phi',
    'muon_phi1' : 'phi',
    'muon_eta0' : 'eta',
    'muon_eta1' : 'eta',
    'dimuon_mass' : 'dilepton_mass'
}
def main():

    indir = "../lo_vs_nlo/input/2019-09-26_checks/"

    acc = acc_from_dir(indir)


    for region in ['cr_2m_j','cr_1m_j','cr_2m_j_noveto_tau','cr_1m_j_noveto_tau','cr_2m_j_noveto_photon','cr_1m_j_noveto_photon']:
        sel_2016 = re.compile(f'MET_2016')
        sel_2017 = re.compile(f'MET_2017')
        sel_2018 = re.compile(f'MET_2018')

        for distribution in axes.keys():
            # fig, ax, rax = make_plot(acc, region=region,distribution=distribution, year=2016, data=data, mc=mc, ylim=(1e-3,1e3), rylim=(0,2),outdir=f'./output/{os.path.basename(indir)}')
            fig, (ax, rax) = plt.subplots(2, 1, figsize=(7,7), gridspec_kw={"height_ratios": (3, 1)}, sharex=True)

            h = acc[distribution]
            h = merge_extensions(h, acc, reweight_pu=('nopu' in distribution))
            # scale_xs_lumi(h)
            h = merge_datasets(h)
            h = h.integrate(h.axis('region'),region)
            s = Style()
            try:
                newax = s.rebin_axes[distribution]
                h = h.rebin(h.axis(newax.name), newax)
            except KeyError:
                pass

            x,y = {}, {}
            h2016 = h[sel_2016].integrate('dataset')
            h2017 = h[sel_2017].integrate('dataset')
            h2018 = h[sel_2018].integrate('dataset')
            try:
                x[2016] = h2016.axis(axes[distribution]).centers()
                y[2016] = h2016.values()[()] / (np.diff(h2016.axis(axes[distribution]).edges())) / lumi(2016)
                x[2017] = h2017.axis(axes[distribution]).centers()
                y[2017] = h2017.values()[()] / (np.diff(h2017.axis(axes[distribution]).edges())) / lumi(2017)
                x[2018] = h2018.axis(axes[distribution]).centers()
                y[2018] = h2018.values()[()] / (np.diff(h2018.axis(axes[distribution]).edges())) / lumi(2018)
            except KeyError:
                continue

            for year in x.keys():
                ax.plot(x[year],y[year],'o-', label=f'{year} / {lumi(year)} fb-1')
                rax.plot(x[year],y[year] / y[2016],'-o')

            ax.legend(title=region)
            ax.set_yscale('log')
            rax.set_xlabel(distribution)
            ax.set_xlabel(distribution)
            rax.set_ylabel('Ratio to 2016')
            ax.set_ylabel('Data cross section / bin')
            loc1 = matplotlib.ticker.MultipleLocator(base=0.2)
            loc2 = matplotlib.ticker.MultipleLocator(base=0.1)
            rax.yaxis.set_major_locator(loc1)
            rax.yaxis.set_minor_locator(loc2)
            rax.grid(axis='y',which='minor',linestyle='--')
            rax.grid(axis='y',which='major',linestyle='--')
            rax.set_ylim(0.5,1.5)
            outname = f'output/bu/{region}_{distribution}.pdf'
            fig.savefig(outname)
            print(f"Saved {outname}")
            plt.close(fig)
        # pprint(dir(h))
    # region='cr_2e_j'
    # mc = re.compile(f'EGamma_2016')
    # data = re.compile(f'EGamma_2017')
    # for distribution in ['recoil', 'ak4_pt0']:
    #     make_plot(acc, region=region,distribution=distribution, year=2016, data=data, mc=mc, ylim=(1e-3,1e3), rylim=(0,2),outdir=f'./output/{os.path.basename(indir)}')
    

    
if __name__ == "__main__":
    main()
