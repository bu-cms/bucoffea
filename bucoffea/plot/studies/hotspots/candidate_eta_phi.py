#!/usr/bin/env python

import os
import re
import sys
from pprint import pprint
from coffea.hist.plot import plot2d
from bucoffea.plot.util import acc_from_dir, merge_extensions, scale_xs_lumi, merge_datasets, lumi
from bucoffea.plot.cr_ratio_plot import cr_ratio_plot
from bucoffea.plot.style import plot_settings
import copy
from collections import defaultdict
from klepto.archives import dir_archive

pjoin = os.path.join
def eta_phi_plot(inpath):
    indir=os.path.abspath(inpath)

    acc = dir_archive(
        indir,
        serialized=True,
        compression=0,
        memsize=1e3
        )
                      
    acc.load('sumw')
    acc.load('sumw_pileup')
    acc.load('sumw2')
    acc.load('nevents')

    outdir = pjoin('./output/', os.path.basename(indir))
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    for year in [2017,2018]:
        data = {
                'cr_1m_j' : f'MET_{year}',
                'cr_2m_j' : f'MET_{year}',
                'cr_1e_j' : f'EGamma_{year}',
                'cr_2e_j' : f'EGamma_{year}',
                # 'cr_g_j' : f'EGamma_{year}',
            }
        for region, datare in data.items():
            distributions = ['ak4_eta_phi']

            if 'e_' in region:
                distributions.append('electron_eta_phi')
            elif 'm_' in region:
                distributions.append('muon_eta_phi')
            for distribution in distributions:
                acc.load(distribution)
                h = copy.deepcopy(acc[distribution])
                h = merge_extensions(h, acc, reweight_pu=('nopu' in distribution))
                scale_xs_lumi(h)
                h = merge_datasets(h)

                h = h.integrate('dataset', datare)
                h = h.integrate(h.axis('region'),region)
                fig, ax, _ = plot2d(h, xaxis='eta')


                ax.text(0., 1., region,
                            fontsize=10,
                            horizontalalignment='left',
                            verticalalignment='top',
                            color='white',
                            transform=ax.transAxes
                        )
                ax.text(1., 0., distribution,
                            fontsize=10,
                            horizontalalignment='right',
                            verticalalignment='bottom',
                            transform=ax.transAxes
                        )
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
                outname = pjoin(outdir,f'{region}_{distribution}_{year}.pdf')
                fig.savefig( outname)
                print(f'Created file {outname}')
    
def eta_phi_plot_photon(inpath):
    '''Create 2D eta-phi plot for photons in VBF photon CR.'''
    indir=os.path.abspath(inpath)

    acc = dir_archive(
        indir,
        serialized=True,
        compression=0,
        memsize=1e3
        )

    acc.load('sumw')
    acc.load('sumw_pileup')
    acc.load('sumw2')
    acc.load('nevents') 

    outdir = pjoin('./output/', os.path.basename(indir))
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    for year in [2017,2018]:
        data = {
                 'cr_g_vbf' : f'EGamma_{year}'
            }
        for region, datare in data.items():
            distributions = ['photon_eta_phi']

            for distribution in distributions:
                acc.load(distribution)
                h = copy.deepcopy(acc[distribution])
                h = merge_extensions(h, acc, reweight_pu=('nopu' in distribution))
                scale_xs_lumi(h)
                h = merge_datasets(h)

                h = h.integrate('dataset', datare)
                h = h.integrate(h.axis('region'),region)
                fig, ax, _ = plot2d(h, xaxis='eta')


                ax.text(0., 1., region,
                            fontsize=10,
                            horizontalalignment='left',
                            verticalalignment='top',
                            color='white',
                            transform=ax.transAxes
                        )
                ax.text(1., 0., distribution,
                            fontsize=10,
                            horizontalalignment='right',
                            verticalalignment='bottom',
                            transform=ax.transAxes
                        )
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
                outname = pjoin(outdir,f'{region}_{distribution}_{year}.pdf')
                fig.savefig( outname)
                print(f'Created file {outname}')
def main():
    inpath = sys.argv[1]
    job_type = sys.argv[2]
    
    if job_type == 'monojet':
        eta_phi_plot(inpath)

    elif job_type == 'photon':
        eta_phi_plot_photon(inpath)


if __name__ == "__main__":
    main()
