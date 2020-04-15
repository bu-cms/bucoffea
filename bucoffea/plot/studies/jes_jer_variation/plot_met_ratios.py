#!/usr/bin/env python

import os
import sys
import re
import numpy as np
from bucoffea.plot.util import merge_datasets, merge_extensions, scale_xs_lumi
from klepto.archives import dir_archive
from matplotlib import pyplot as plt
from pprint import pprint
import mplhep as hep

pjoin = os.path.join

pretty_var_label = {
    '_jerup' : 'JER up',
    '_jerdown' : 'JER down',
    '_jesup' : 'JES up',
    '_jesdown' : 'JES down'
}

pretty_process_label = {
    'cr_1e_vbf' : r'$W(e\nu)$',
    'cr_1m_vbf' : r'$W(\mu\nu)$',
    'cr_2e_vbf' : r'$Z(ee)$',
    'cr_2m_vbf' : r'$Z(\mu\mu)$'
}

def plot_met_ratios(acc, tag, dataset_regex, regions, variation, outtag):
    acc.load('met_ratio')
    h = acc['met_ratio']

    h = merge_extensions(h, acc, reweight_pu=False)
    scale_xs_lumi(h)
    h = merge_datasets(h)

    # Pick the regions and datasets
    h1 = h.integrate('region', f'{regions[0]+variation}').integrate('dataset', re.compile(dataset_regex))
    h2 = h.integrate('region', f'{regions[1]+variation}').integrate('dataset', re.compile(dataset_regex))

    edges = h1.axes()[0].edges()

    h1_vals = h1.values()
    h2_vals = h2.values()

    # h1 and h2 contain histograms for different MET pt bins
    for metbin in h1_vals.keys():
        arr1, arr2 = h1_vals[metbin], h2_vals[metbin]
        # Store the maximum in each array
        max1, max2 = max(arr1), max(arr2)
        fig, ax = plt.subplots()
        hep.histplot(arr1, edges, ax=ax, histtype='step', label=pretty_process_label[regions[0]])
        hep.histplot(arr2, edges, ax=ax, histtype='step', label=pretty_process_label[regions[1]])

        ax.legend(title=f'Variation: {pretty_var_label[variation]}')
        ax.set_xlabel(r'Varied MET $p_T$ / Nominal MET $p_T$ - 1')
        ax.set_ylabel('Counts')
        ax.set_ylim(0, np.maximum(max1,max2)*1.2)
        ax.set_title(metbin[0])

        # Save figure
        outpath = f'./output/{outtag}/met_ratios'
        if not os.path.exists(outpath):
            os.makedirs(outpath)
        
        outfile = pjoin(outpath, f'{tag}{variation}_{metbin[0]}.pdf')

        fig.savefig(outfile)
        print(f'Figure saved: {outfile}')

def main():
    inpath = sys.argv[1]
    acc = dir_archive(
        inpath,
        memsize=1e3,
        compression=0,
        serialized=True
    )
    acc.load('sumw')
    acc.load('sumw2')

    if inpath.endswith('/'):
        outtag = inpath.split('/')[-2]
    else:
        outtag = inpath.split('/')[-1]

    variations = ['_jerup', '_jerdown', '_jesup', '_jesdown']
    tag_regex_regions = {
        'wjets17' : {'regex' : 'WJetsToLNu.*2017', 'regions' : ('cr_1e_vbf', 'cr_1m_vbf')},
        'wjets18' : {'regex' : 'WJetsToLNu.*2018', 'regions' : ('cr_1e_vbf', 'cr_1m_vbf')},
        'zjets17' : {'regex' : 'DYJetsToLL.*2017', 'regions' : ('cr_2e_vbf', 'cr_2m_vbf')},
        'zjets18' : {'regex' : 'DYJetsToLL.*2018', 'regions' : ('cr_2e_vbf', 'cr_2m_vbf')},
    }
    
    for tag, info in tag_regex_regions.items():
        regex, regions = info.values()
        for var in variations:
            plot_met_ratios(acc, tag=tag, dataset_regex=regex, regions=regions, variation=var, outtag=outtag)

    plt.close()

if __name__ == '__main__':
    main()