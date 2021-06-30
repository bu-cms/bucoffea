#!/usr/bin/env python

import os
import sys
import re
import uproot
import numpy as np

from matplotlib import pyplot as plt
from coffea import hist
from bucoffea.helpers.paths import bucoffea_path
from bucoffea.plot.util import merge_datasets, merge_extensions, scale_xs_lumi, fig_ratio
from klepto.archives import dir_archive
from pprint import pprint

pjoin = os.path.join

pretty_dataset_label = {
    'cr_2e_vbf': {
        'WJetsToLNu' : r'$W(e\nu)$ {year}',
        'DYJetsToLL' : r'$Z(ee)$ {year}',
    },
    'cr_2m_vbf': {
        'WJetsToLNu' : r'$W(\mu\nu)$ {year}',
        'DYJetsToLL' : r'$Z(\mu\mu)$ {year}',
    },
    'sr_vbf_no_veto_all': {
        'WJetsToLNu' : r'$W(\ell\nu)$ {year}',
        'DYJetsToLL' : r'$Z(\ell\ell)$ {year}',
    },
}

COLOR_CODE_FOR_VETOWEIGHTS = {
    'ele_id.*': 1,
    'ele_reco.*' : 2,
    'muon_id.*' : 3,
    'muon_iso.*' : 4,
    'tau_id.*' : 5 
}

def plot_lepton_sf_uncertainties(acc, outtag, sftag, region, dataset, outrootfile):
    '''Plot lepton ID/ISO/RECO SF uncertainties as a function of mjj.'''
    distribution = f'mjj_{sftag}'
    acc.load(distribution)
    h = acc[distribution]

    h = merge_extensions(h, acc)
    scale_xs_lumi(h)
    h = merge_datasets(h)

    if 'mjj' in distribution:
        mjj_bins = [200., 400., 600., 900., 1200., 1500., 2000., 2750., 3500., 5000.]
        mjj_ax = hist.Bin('mjj', r'$M_{jj} \ (GeV)$', mjj_bins)
        h = h.rebin('mjj', mjj_ax)            

    h = h.integrate('region', region)

    outdir = f'./output/{outtag}'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    data_err_opts = {
        'linestyle':'none',
        'marker': '.',
        'markersize': 10.,
    }

    for year in [2017, 2018]:
        _h = h.integrate('dataset', re.compile(f'{dataset}.*{year}'))
        fig, ax, rax = fig_ratio()
        hist.plot1d(_h,ax=ax,overlay='variation')
        
        ax.text(0., 1., pretty_dataset_label[region][dataset].format(year=year),
            fontsize=14,
            ha='left',
            va='bottom',
            transform=ax.transAxes
        )
        
        ax.text(1., 1., sftag,
            fontsize=14,
            ha='right',
            va='bottom',
            transform=ax.transAxes
        )

        ax.set_yscale('log')
        ax.set_ylim(1e-2,1e6)

        hist.plotratio(
            _h.integrate('variation', 'up'),
            _h.integrate('variation', 'nom'),
            ax=rax,
            unc='num',
            label='Up',
            error_opts=data_err_opts
        )

        hist.plotratio(
            _h.integrate('variation', 'down'),
            _h.integrate('variation', 'nom'),
            ax=rax,
            unc='num',
            label='Down',
            error_opts=data_err_opts,
            clear=False
        )

        rax.legend(ncol=2)
        rax.grid(True)
        if 'DY' in dataset:
            rax.set_ylim(0.9,1.1)
        else:
            rax.set_ylim(0.95,1.05)
        rax.set_ylabel('Variation')

        outpath = pjoin(outdir, f'{dataset}_{sftag}_uncs_{year}.pdf')        
        fig.savefig(outpath)
        plt.close(fig)

        print(f'File saved: {outpath}')

        # Save to output ROOT file
        edges = _h.axis('mjj').edges()

        if re.match('cr_(\d)m_vbf', region):
            if 'DY' in dataset:
                proctag = 'Zmumu'
            else:
                proctag = 'Wmunu'
        elif re.match('cr_(\d)e_vbf', region):
            if 'DY' in dataset:
                proctag = 'Zee'
            else:
                proctag = 'Wenu'
        else:
            raise RuntimeError(f'Could not determine process tag for region: {region}')

        r_up = _h.integrate('variation', 'up').values()[()] / _h.integrate('variation', 'nom').values()[()]
        r_down = _h.integrate('variation', 'down').values()[()] / _h.integrate('variation', 'nom').values()[()]
        outrootfile[f'{proctag}_{sftag}_{year}_up'] = (r_up, edges)
        outrootfile[f'{proctag}_{sftag}_{year}_down'] = (r_down, edges)

def plot_lepton_veto_weight_uncs(acc, outtag, region, dataset, outrootfile):
    distribution = f'mjj_veto_weight'
    acc.load(distribution)
    h = acc[distribution]

    h = merge_extensions(h, acc)
    scale_xs_lumi(h)
    h = merge_datasets(h)

    if 'mjj' in distribution:
        mjj_bins = [200., 400., 600., 900., 1200., 1500., 2000., 2750., 3500., 5000.]
        mjj_ax = hist.Bin('mjj', r'$M_{jj} \ (GeV)$', mjj_bins)
        h = h.rebin('mjj', mjj_ax)            

    h = h.integrate('region', region)

    outdir = f'./output/{outtag}'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    data_err_opts = {
        'linestyle':'none',
        'marker': '.',
        'markersize': 10.,
    }

    for year in [2017, 2018]:
        _h = h.integrate('dataset', re.compile(f'{dataset}.*{year}'))
        fig, ax, rax = fig_ratio()
        hist.plot1d(_h,ax=ax,overlay='variation')

        ax.set_yscale('log')
        ax.set_ylim(1e-2,1e6)
        ax.legend(ncol=2, title='Veto Weight Variations')

        ax.text(0., 1., pretty_dataset_label[region][dataset].format(year=year),
            fontsize=14,
            ha='left',
            va='bottom',
            transform=ax.transAxes
        )

        variations = _h.axis('variation').identifiers()
        for v in variations:
            if v.name == 'nominal':
                continue
            opts = data_err_opts
            opts['color'] = 'C0'
            for name, code in COLOR_CODE_FOR_VETOWEIGHTS.items():
                if re.match(name, str(v)):
                    opts['color'] = f'C{code}'
                    break
            hist.plotratio(
                _h.integrate('variation', v),
                _h.integrate('variation', 'nominal'),
                ax=rax,
                unc='num',
                label=v,
                error_opts=opts,
                clear=False
            )

        rax.legend(ncol=2)
        rax.grid(True)
        rax.set_ylim(0.97,1.03)
        rax.set_ylabel('Variation')

        outpath = pjoin(outdir, f'{dataset}_vetoweight_uncs_{year}.pdf')        
        fig.savefig(outpath)
        plt.close(fig)

        print(f'File saved: {outpath}')

        # Save to output ROOT file
        edges = _h.axis('mjj').edges()

        if re.match('cr_(\d)m_vbf', region):
            if 'DY' in dataset:
                proctag = 'Zmumu'
            else:
                proctag = 'Wmunu'
        elif re.match('cr_(\d)e_vbf', region):
            if 'DY' in dataset:
                proctag = 'Zee'
            else:
                proctag = 'Wenu'
        elif re.match('.*no_veto.*', region):
            if 'DY' in dataset:
                proctag = 'Zll'
            else:
                proctag = 'Wlnu'
        else:
            raise RuntimeError(f'Could not determine process tag for region: {region}')

        for v in variations:
            r = _h.integrate('variation', v).values()[()] / _h.integrate('variation','nominal').values()[()]
            outrootfile[f'{proctag}_vetow_{year}_{v}'] = (r, edges)

def main():
    inpath = sys.argv[1]
    acc = dir_archive(inpath)
    acc.load('sumw')
    acc.load('sumw_pileup')
    acc.load('nevents')

    outtag = re.findall('merged_.*', inpath)[0].replace('/','')

    # Lepton ID / ISO SF
    lepton_sf_uncs = {
        'ele_id' : [
            {'dataset' : 'WJetsToLNu', 'region': 'cr_1e_vbf'},
            {'dataset' : 'DYJetsToLL', 'region': 'cr_2e_vbf'},
        ], 
        'ele_reco' : [
            {'dataset' : 'WJetsToLNu', 'region': 'cr_1e_vbf'},
            {'dataset' : 'DYJetsToLL', 'region': 'cr_2e_vbf'},
        ], 
        'muon_id' : [
            {'dataset' : 'WJetsToLNu', 'region': 'cr_1m_vbf'},
            {'dataset' : 'DYJetsToLL', 'region': 'cr_2m_vbf'},
        ], 
        'muon_iso' : [
            {'dataset' : 'WJetsToLNu', 'region': 'cr_1m_vbf'},
            {'dataset' : 'DYJetsToLL', 'region': 'cr_2m_vbf'},
        ], 

    }

    outdir = f'./output/{outtag}'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    outrootfile = uproot.recreate(pjoin(outdir,  'lepton_sf_uncs.root'))
    for unc, data in lepton_sf_uncs.items():
        for entry in data:
            try:
                plot_lepton_sf_uncertainties(acc, outtag, 
                    sftag=unc,
                    region=entry['region'],
                    dataset=entry['dataset'],
                    outrootfile=outrootfile
                    )
            except (KeyError, AssertionError) as e:
                print(f'Skipping: {unc}')
                continue

    # Lepton veto uncertainties
    outrootfile_veto = uproot.recreate(pjoin(outdir,  'lepton_vetow_uncs.root'))
    datasets_for_veto = [
        {'dataset' : 'DYJetsToLL', 'region' : 'sr_vbf_no_veto_all'},
        {'dataset' : 'WJetsToLNu', 'region' : 'sr_vbf_no_veto_all'},
    ]

    for entry in datasets_for_veto:
        try:
            plot_lepton_veto_weight_uncs(acc, outtag,
                region=entry['region'],
                dataset=entry['dataset'],
                outrootfile=outrootfile_veto
                )
        except KeyError:
            print(f'Skipping: {entry["dataset"]}')
            continue


if __name__ == '__main__':
    main()