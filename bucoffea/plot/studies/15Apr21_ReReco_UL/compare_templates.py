#!/usr/bin/env python
import os
import re
import sys
import argparse
import datetime
import numpy as np

from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
from coffea import hist
from bucoffea.plot.util import merge_datasets, merge_extensions, scale_xs_lumi, fig_ratio
from klepto.archives import dir_archive
from pprint import pprint

pjoin = os.path.join

def parse_cli():
    parser = argparse.ArgumentParser()
    parser.add_argument('--inpath1', help='Path to merged acc 1.')
    parser.add_argument('--inpath2', help='Path to merged acc 2.')
    parser.add_argument('--tag1', help='Accumulator tag for acc 1.')
    parser.add_argument('--tag2', help='Accumulator tag for acc 2.')
    args = parser.parse_args()
    return args

def preprocess(h, acc, distribution, region, dataset_regex):
    h = merge_extensions(h, acc)
    scale_xs_lumi(h)
    h = merge_datasets(h)

    if distribution ==  'mjj':
        mjj_ax = hist.Bin('mjj', r'$M_{jj} \ (GeV)$', [200., 400., 600., 900., 1200., 1500., 2000., 2750., 3500., 5000.])
        h = h.rebin('mjj', mjj_ax)

    h = h.integrate('region', region).integrate('dataset', dataset_regex)
    return h

def compare_templates(acc_dict, outtag, dataset_tag, dataset_regex, distribution='mjj', region='sr_vbf_no_veto_all'):
    '''Compare the two templates.'''
    histos = {}
    for key, acc in acc_dict.items():
        acc.load(distribution)
        histos[key] = preprocess(
            acc[distribution],
            acc,
            distribution=distribution,
            region=region,
            dataset_regex=dataset_regex
        )

    fig, ax, rax = fig_ratio()
    keys = list(histos.keys())
    hist.plot1d(histos[keys[0]], ax=ax)
    hist.plot1d(histos[keys[1]], ax=ax, clear=False)

    ax.set_yscale('log')
    ax.set_ylim(1e0,1e6)

    ax.legend(labels=keys)

    ax.text(0.,1.,'VBF SR',
        fontsize=14,
        ha='left',
        va='bottom',
        transform=ax.transAxes
    )

    year = dataset_tag.split('_')[-1]
    assert int(year) in [2017, 2018]

    ax.text(1.,1.,year,
        fontsize=14,
        ha='right',
        va='bottom',
        transform=ax.transAxes
    )

    data_err_opts = {
        'linestyle':'none',
        'marker': '.',
        'markersize': 10.,
        'color':'k',
    }

    hist.plotratio(
        histos[keys[0]],
        histos[keys[1]],
        ax=rax,
        unc='num',
        error_opts=data_err_opts
    )

    rax.grid(True)
    rax.set_ylim(0.9,1.1)
    rax.set_ylabel(f'{keys[0]} / {keys[1]}')

    loc = MultipleLocator(0.05)
    rax.yaxis.set_major_locator(loc)

    outdir = f'./output/{outtag}'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outpath = pjoin(outdir, f'{dataset_tag}_{distribution}.pdf')
    fig.savefig(outpath)
    plt.close(fig)

    print(f'File saved: {outpath}')

def main():
    args = parse_cli()
    inpaths = {}

    inpaths[args.tag1] = args.inpath1
    inpaths[args.tag2] = args.inpath2

    acc_dict = {}
    for tag, inpath in inpaths.items():
        acc_dict[tag] = dir_archive(inpath)
        acc_dict[tag].load('sumw')
        acc_dict[tag].load('sumw_pileup')
        acc_dict[tag].load('nevents')

    outtag = f'{datetime.date.today().strftime("%m_%d_%y")}_{args.tag1}_{args.tag2}'

    for year in [2017, 2018]:
        dataset_regex = re.compile(f'(ZJetsToNuNu|WJetsToLNu|EWK|Top|Diboson).*{year}')
        dataset_tag = f'mc_all_{year}'

        compare_templates(acc_dict, outtag,
            dataset_tag=dataset_tag,
            dataset_regex=dataset_regex,
        )

if __name__ == '__main__':
    main()