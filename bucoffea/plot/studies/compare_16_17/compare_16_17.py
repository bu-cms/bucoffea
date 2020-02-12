#!/usr/bin/env python
import os
import re
from pprint import pprint

import matplotlib.ticker
import numpy as np
import uproot
from coffea import hist
from coffea.util import load
from matplotlib import pyplot as plt

from bucoffea.plot.stack_plot import Style, make_plot
from bucoffea.plot.util import (acc_from_dir, lumi, merge_datasets,
                                merge_extensions, scale_xs_lumi)
from klepto.archives import dir_archive


def main():

    for mode in ['monojet', 'zpt']:
        if mode == 'monojet':
            regions = {
                    # 'ch1' : 'sr_j',
                    'ch2' : 'cr_2m_j',
                    'ch3' : 'cr_1m_j',
                    # # 'ch4' : 'cr_g_j',
                    # # 'ch5' : 'cr_2e_j',
                    'ch6' : 'cr_1e_j'
                        }
            fitfile = 'fitDiagnostics.root'
            bins = [ 250,  280,  310,  340,  370,  400,  430,  470,  510, 550,  590,  640,  690,  740,  790,  840,  900,  960, 1020, 1090, 1160, 1250, 1400]
        elif mode == 'zpt':
            regions = {
            # 'ch1' : 'sr_j',
            'monojet_singlemu' : 'cr_1m_j',
            # 'ch4' : 'cr_g_j',
            # 'ch5' : 'cr_2e_j',
            'monojet_singleel' : 'cr_1e_j'
            }
            fitfile = 'fitDiagnostics_zpt.root'
            bins = [250,275,300,350,400,450,500,650,800,1150,1400]  
        tmp = {}
        for k,v in regions.items():
            tmp[v] = k
        regions.update(tmp)

        inpath = "input/2020-01-07_fine_recoil_bins_v4"

        acc = dir_archive(
                        inpath,
                        serialized=True,
                        compression=0,
                        memsize=1e3,
                        )
        acc.load('recoil')
        acc.load('sumw')
        acc.load('sumw_pileup')
        acc.load('nevents')
        for region in [x for x in regions.keys() if 'r_' in x]:
            if '_1m_' in region or '_2m_' in region:
                sel_2016 = re.compile(f'MET_2016')
                sel_2017 = re.compile(f'MET_2017')
                sel_2018 = re.compile(f'MET_2018')
            else:
                sel_2016 = re.compile(f'EGamma_2016')
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

                # Rebin
                newax = hist.Bin('recoil','Recoil (GeV)', bins)
                h = h.rebin(h.axis(newax.name), newax)

                # Retrieve input from fit diagnostics
                f = uproot.open(fitfile)
                tg2016 = f['shapes_prefit'][regions[region]]['data']

                # Prepare x and y values for easy plotting
                x,y = {}, {}
                # x[2016] = np.r_[0.5*(h2016.bins[:,0] + h2016.bins[:,1]),1500., 1700., 1900.]
                x['2016orig'] = np.r_[tg2016.xvalues]
                y['2016orig'] = np.r_[tg2016.yvalues/ lumi(2016)]
                h2016 = h[sel_2016].integrate('dataset')
                h2017 = h[sel_2017].integrate('dataset')
                h2018 = h[sel_2018].integrate('dataset')

                x[2016] = h2016.axis('recoil').centers()
                y[2016] = h2016.values()[()] / (np.diff(h2016.axis('recoil').edges())) / lumi(2016)
                x[2017] = h2017.axis('recoil').centers()
                y[2017] = h2017.values()[()] / (np.diff(h2017.axis('recoil').edges())) / lumi(2017)
                x[2018] = h2018.axis('recoil').centers()
                y[2018] = h2018.values()[()] / (np.diff(h2018.axis('recoil').edges())) / lumi(2018)


                # Actual plotting
                for year in x.keys():
                    print(x[year])
                    ax.plot(x[year],y[year],'o-', label=f'{year} ({sum(y[year]*np.diff(bins)):.2e})')
                    rax.plot(x[year],y[year] / y['2016orig'],'-o')

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
                fig.savefig(f'output/{mode}_{region}.pdf')

    
if __name__ == "__main__":
    main()
