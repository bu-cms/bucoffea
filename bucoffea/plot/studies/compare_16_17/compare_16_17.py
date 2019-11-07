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

def main():

    indir = "../lo_vs_nlo/input/2019-09-23_photon_overlap"

    acc = acc_from_dir(indir)


    for region in ['cr_2m_j','cr_2e_j','cr_1m_j', 'cr_1e_j', 'cr_g_j']:
        if '_1m_' in region or '_2m_' in region:
            sel_2017 = re.compile(f'MET_2017')
            sel_2018 = re.compile(f'MET_2018')
        else:
            sel_2017 = re.compile(f'EGamma_2017')
            sel_2018 = re.compile(f'EGamma_2018')

        for distribution in ['recoil']:
            # fig, ax, rax = make_plot(acc, region=region,distribution=distribution, year=2016, data=data, mc=mc, ylim=(1e-3,1e3), rylim=(0,2),outdir=f'./output/{os.path.basename(indir)}')
            fig, (ax, rax) = plt.subplots(2, 1, figsize=(7,7), gridspec_kw={"height_ratios": (3, 1)}, sharex=True)

            h = acc['recoil']
            h = merge_extensions(h, acc, reweight_pu=('nopu' in distribution))
            # scale_xs_lumi(h)
            h = merge_datasets(h)
            h = h.integrate(h.axis('region'),region)
            s = Style()
            try:
                newax = s.get_binning(distribution, region)
                h = h.rebin(h.axis(newax.name), newax)
            except KeyError:
                pass

            f = uproot.open('fitDiagnostics.root')
            h2016 = f['shapes_prefit'][regions[region]]['data']

            x,y = {}, {}
            # x[2016] = np.r_[0.5*(h2016.bins[:,0] + h2016.bins[:,1]),1500., 1700., 1900.]
            pprint(dir(h2016))
            x[2016] = np.r_[h2016.xvalues,1500., 1700., 1900.]
            y[2016] = np.r_[h2016.yvalues/ lumi(2016),0,0,0]
            h2017 = h[sel_2017].integrate('dataset')
            h2018 = h[sel_2018].integrate('dataset')
            x[2017] = h2017.axis('recoil').centers()
            y[2017] = h2017.values()[()] / (np.diff(h2017.axis('recoil').edges())) / lumi(2017)
            x[2018] = h2018.axis('recoil').centers()
            y[2018] = h2018.values()[()] / (np.diff(h2018.axis('recoil').edges())) / lumi(2018)


            for year in x.keys():
                ax.plot(x[year],y[year],'o-', label=f'{year}')
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
            fig.savefig(f'output/{region}.pdf')
        # pprint(dir(h))
    # region='cr_2e_j'
    # mc = re.compile(f'EGamma_2016')
    # data = re.compile(f'EGamma_2017')
    # for distribution in ['recoil', 'ak4_pt0']:
    #     make_plot(acc, region=region,distribution=distribution, year=2016, data=data, mc=mc, ylim=(1e-3,1e3), rylim=(0,2),outdir=f'./output/{os.path.basename(indir)}')
    

    
if __name__ == "__main__":
    main()
