#!/usr/bin/env python

import os
import sys
import re
import numpy as np

from coffea import hist
from scipy.stats import distributions
from bucoffea.plot.util import merge_datasets, merge_extensions, scale_xs_lumi
from matplotlib import pyplot as plt
from klepto.archives import dir_archive
from pprint import pprint

pjoin = os.path.join

def plot_strip_size(acc, outtag, distribution, region, dataset='MET_2017', etaslice=None):
    acc.load(distribution)
    h = acc[distribution]

    h = merge_extensions(h, acc, reweight_pu=False)
    scale_xs_lumi(h)
    h = merge_datasets(h)

    h = h.integrate('dataset', dataset).integrate('region', region)

    # If the histogram is also split by jet eta, take the slice we're interested in
    # Typically: 2.9 < |eta| < 3.25 and 3.25 < |eta| < 5.0
    if etaslice is not None:
        h = h.integrate('jeta', etaslice)
    
    fig, ax = plt.subplots()
    hist.plot2d(h, ax=ax, xaxis='centraletastripsize')

    if re.match('cr_2m_vbf.*', region):
        regiontag = r'$Z(\mu\mu)$'
    elif region == 'cr_vbf_qcd':
        regiontag = 'QCD CR'
    else:
        raise RuntimeError(f'Check region: {region}')

    ax.text(0., 1., f'{regiontag}, Data',
        fontsize=14,
        ha='left',
        va='bottom',
        transform=ax.transAxes
    )

    ax.text(1., 1., r'HF Jets, $p_T > 80 \ GeV$',
        fontsize=14,
        ha='right',
        va='bottom',
        transform=ax.transAxes
    )

    if distribution == 'ak4_hfcentral_adjacent_etastripsize0':
        ax.set_xlabel('Leading Jet Central Strip Size')
        ax.set_ylabel('Leading Jet Adjacent Strip Size')
    elif distribution == 'ak4_hfcentral_adjacent_etastripsize1':
        ax.set_xlabel('Trailing Jet Central Strip Size')
        ax.set_ylabel('Trailing Jet Adjacent Strip Size')

    if etaslice in [slice(2.9, 3.25), slice(3.25, 5.0)]:
        tagforslice = {
            (2.9, 3.25): r'$2.9 < |\eta| < 3.25$',
            (3.25, 5.0): r'$3.25 < |\eta| < 5.0$',
        }
        
        for (lo, hi), _tag in tagforslice.items():
            if (etaslice.start, etaslice.stop) == (lo, hi):
                _tagforslice = _tag
                break

        ax.text(1., 0., _tagforslice,
            fontsize=14,
            ha='right',
            va='bottom',
            transform=ax.transAxes
        )

    # Save figure
    outdir = f'./output/{outtag}/2d'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    if etaslice in [slice(2.9, 3.25), slice(3.25, 5.0)]:
        etaslicetag = f'_eta_{str(etaslice.start).replace(".", "p")}_{str(etaslice.stop).replace(".", "p")}'
    elif etaslice == slice(2.9, 5.0):
        etaslicetag= '_eta_combined'
    else:
        etaslicetag = ''

    outpath = pjoin(outdir, f'MET_2017_{region}_{distribution}{etaslicetag}.pdf')
    fig.savefig(outpath)
    plt.close(fig)
    print(f'File saved: {outpath}')

def main():
    inpath = sys.argv[1]
    acc = dir_archive(inpath)
    acc.load('sumw')
    acc.load('sumw2')

    outtag = re.findall('merged_.*', inpath)[0].replace('/','')

    distributions = [
        'ak4_hfcentral_adjacent_etastripsize',
        'ak4_hfcentral_adjacent_etastripsize0',
        'ak4_hfcentral_adjacent_etastripsize1',
    ]

    regions = [
        # 'cr_2m_vbf_relaxed_sel',
        'cr_vbf_qcd',
    ]

    etaslices = [
        slice(2.9, 3.25),
        slice(3.25, 5.0),
        slice(2.9, 5.0),
    ]

    for region in regions:
        for distribution in distributions:
            for etaslice in etaslices:
                plot_strip_size(acc, outtag,
                        distribution=distribution,
                        region=region,
                        etaslice=etaslice
                    )

if __name__ == '__main__':
    main()