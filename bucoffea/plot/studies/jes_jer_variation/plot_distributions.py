#!/usr/bin/env python

import os
import sys
import re
from bucoffea.plot.util import merge_datasets, merge_extensions, scale_xs_lumi
from coffea import hist
from matplotlib import pyplot as plt
from klepto.archives import dir_archive
from pprint import pprint

pjoin = os.path.join

# Plot aesthetics
colors = plt.rcParams['axes.prop_cycle'].by_key()['color'][:5]

recoil_bins_2016 = [250,  280,  310,  340,  370,  400,  430,  470,  510, 550,  590,  640,  690,  740,  790,  840,  900]

REBIN = {
    'mjj' : hist.Bin('mjj', r'$M_{jj}$ (GeV)', list(range(200,800,300)) + list(range(800,2000,400)) + [2000, 2750, 3500]),
    'recoil' : hist.Bin('recoil','Recoil (GeV)', recoil_bins_2016),
    'ak4_pt0' : hist.Bin('jetpt',r'Leading AK4 jet $p_{T}$ (GeV)',list(range(80,700,20)) ),
    'ak4_pt1' : hist.Bin('jetpt',r'Trailing AK4 jet $p_{T}$ (GeV)',list(range(40,400,20)) ),
    'ak4_pt' : hist.Bin('jetpt',r'All AK4 jet $p_{T}$ (GeV)',list(range(100,700,20)) )
}

def plot_dist(acc, tag, param, regex, region, outtag):
    '''
    Given the input accumulator, plot distribution of param for ALL variations for given dataset and region.
    =================
    PARAMETERS:
    =================
    acc:    Input accumulator containing all histograms.
    tag:    Tag representing the process.
    param:  The parameter to be plotted.
    regex:  Regular expression matching the dataset name.
    region: The region from which data will be taken.
    outtag: Tag for naming output directory.
    '''
    acc.load(param)
    h = acc[param]

    # Rebin, if neccessary
    if param in REBIN.keys():
        new_bin = REBIN[param]
        axis_name = 'jetpt' if 'ak4_pt' in param else param
        h = h.rebin(h.axis(axis_name), new_bin)
    
    # Merge extensions/datasets, scale w.r.t. xs and lumi
    h = merge_extensions(h, acc, reweight_pu=False)
    scale_xs_lumi(h)
    h = merge_datasets(h)

    # Pick the relevant region and dataset
    h = h[re.compile(regex)].integrate('dataset')[re.compile(f'{region}.*')]

    var_to_label = {
        '' : 'Nominal',
        '_jerup' : 'JER up',
        '_jerdown' : 'JER down',
        '_jesup' : 'JES up',
        '_jesdown' : 'JES down',
    }

    fig, (ax, rax) = plt.subplots(2, 1, figsize=(7,7), gridspec_kw={"height_ratios": (2,1)}, sharex=True)
    
    # Plot for each variation
    hist.plot1d(h, ax=ax, overlay='region')
    ax.set_ylabel('Counts')        

    # Calculate and plot ratios
    h_nom = h.integrate('region', region).values()[()]
    centers = h.axis(axis_name).centers()
    for idx, (var, label) in enumerate(var_to_label.items() ):
        if var == '':
            continue
        ratio = h.integrate('region', f'{region}{var}').values()[()] / h_nom
        rax.plot(centers, ratio, 'o', color=colors[idx], label=label)

    # Display dataset name and type (old or new) as title 
    dataset_name = regex.replace('.*', '_HT_')
    if '06JanInput' in outtag:
        dataset_tag = '06Jan20'
    else:        
        dataset_tag = '19Feb20'
    title = f'{dataset_name}: {dataset_tag}'
    ax.set_title(title)
    ax.set_xlabel('')

    rax.grid(True)
    rax.legend(ncol=2)
    rax.set_ylim(0.7, 1.3)
    rax.set_xlabel(h.axis(axis_name).label)
    rax.set_ylabel('Var / Nom')

    outdir = f'./output/{outtag}/distributions'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    outpath = pjoin(outdir, f'{tag}_{param}.pdf')
    fig.savefig(outpath)
    plt.close()

    print(f'File saved: {outpath}')

def main():
    inpath = sys.argv[1]
    
    # Handle optional parameter argument
    param_to_plot = None
    if (len(sys.argv) > 2):
        param_to_plot = sys.argv[2]

    # Get outtag for output directory naming
    if inpath.endswith('/'):
        outtag = inpath.split('/')[-2]
    else:
        outtag = inpath.split('/')[-1]

    acc = dir_archive(
        inpath,
        memsize=1e3,
        compression=0,
        serialized=True
    )

    acc.load('sumw')
    acc.load('sumw2')

    params = ['mjj', 'recoil', 'ak4_pt', 'ak4_pt0', 'ak4_pt1']
    dataset_info = {
        'wjets17' : {'regex': 'WJetsToLNu.*2017', 'region': 'sr_vbf'},
        'wjets18' : {'regex': 'WJetsToLNu.*2018', 'region': 'sr_vbf'},
        'zjets17' : {'regex': 'ZJetsToNuNu.*2017', 'region': 'sr_vbf'},
        'zjets18' : {'regex': 'ZJetsToNuNu.*2018', 'region': 'sr_vbf'},
        'gjets17' : {'regex': 'GJets_DR-0p4.*2017', 'region': 'cr_g_vbf'},
    }

    for param in params:
        if param_to_plot:
            if param != param_to_plot:
                continue
        for tag, info in dataset_info.items():
            regex = info['regex']
            region = info['region']
            plot_dist(acc, tag=tag, param=param, regex=regex, region=region, outtag=outtag)

if __name__ == '__main__':
    main()

