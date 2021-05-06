#!/usr/bin/env python
import os
import re
import sys
import tabulate
import numpy as np
import mplhep as hep

from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
from coffea import hist
from coffea.hist import poisson_interval
from bucoffea.plot.util import merge_datasets, merge_extensions, scale_xs_lumi, fig_ratio, lumi
from bucoffea.helpers.paths import bucoffea_path
from klepto.archives import dir_archive
from pprint import pprint

pjoin = os.path.join

HEADERS = [
    'Mjj Bin Range (GeV)',
    '200-400',
    '400-600',
    '600-900',
    '900-1200',
    '1200-1500',
    '1500-2000',
    '2000-2750',
    '2750-3500',
    '>3500',
]

def tabulate_bkg_predictions(acc,outtag,region='sr_vbf_no_veto_all', isUL=True):
    '''Tabulate bkg predictions from MC in the SR.'''
    distribution='mjj'
    acc.load(distribution)
    h = acc[distribution]

    h = merge_extensions(h, acc)
    scale_xs_lumi(h)
    h = merge_datasets(h)

    # Update mjj binning
    mjj_ax = hist.Bin('mjj', r'$M_{jj} \ (GeV)$', [200., 400., 600., 900., 1200., 1500., 2000., 2750., 3500.])
    h = h.rebin('mjj', mjj_ax)

    h = h.integrate('region', region)

    outdir = f'./output/prediction_tables'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    for year in [2017, 2018]:
        datasets = [
            ('QCD Z(vv)',   re.compile(f'ZJetsToNuNu.*{year}')),
            ('EWK Z(vv)',   re.compile(f'EWKZ2Jets.*ZToNuNu.*{year}')),
            ('QCD W(lv)',  re.compile(f'WJetsToLNu_HT.*{year}')),
            ('EWK W(lv)',  re.compile(f'EWKW2Jets.*{year}')),
            ('QCD Z(ll)', re.compile(f'DYJetsToLL.*{year}')),
            ('EWK Z(ll)', re.compile(f'EWKZ2Jets.*ZToLL.*{year}')),
            ('Top', re.compile(f'Top.*{year}')),
            ('Diboson', re.compile(f'Diboson.*{year}')),
        ]

        table = []

        for tag, regex in datasets:
            _h = h.integrate('dataset', regex)
            sumw, sumw2 = _h.values(overflow='over', sumw2=True)[()]

            # Calculate the statistical uncertainty on the MC samples
            err = np.abs(poisson_interval(sumw, sumw2) - sumw)
            # Get the average error
            avgerr = 0.5 * (err[0] + err[1])

            data = [tag]

            # Per mjj bin table
            for _sumw, _err in zip(sumw, avgerr):
                data.append(f'{_sumw if _sumw > 0 else 0.:.1f} +- {_err if not (np.isnan(_err) or _sumw < 0) else 0.:.1f}')

            table.append(data)

        # Write the table to an output file
        outpath = pjoin(outdir, f'table_{year}_{"ul" if isUL else "eoy"}.txt')
        with open(outpath, 'w+') as f:
            f.write(f'Job tag: {outtag}')
            f.write('\n')
            f.write(f'{"UL" if isUL else "EOY"}')
            f.write('\n')
            
            f.write(
                tabulate.tabulate(table, headers=HEADERS, floatfmt=".1f")
            )

        print(f'File saved: {outpath}')

def main():
    inpath = sys.argv[1]
    acc = dir_archive(inpath)
    acc.load('sumw')
    acc.load('sumw_pileup')
    acc.load('nevents')

    outtag = re.findall('merged_.*', inpath)[0].replace('/','')

    tabulate_bkg_predictions(acc, outtag, isUL='ULv8' in outtag)

if __name__ == '__main__':
    main()