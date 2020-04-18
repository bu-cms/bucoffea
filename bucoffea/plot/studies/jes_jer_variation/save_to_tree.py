#!/usr/bin/env python

import os 
import sys
import re
import uproot
from bucoffea.plot.util import merge_datasets, merge_extensions, scale_xs_lumi
from coffea import hist
from klepto.archives import dir_archive
from pprint import pprint

pjoin = os.path.join

# List of variables to be saved and the regular expression
# for dataset used for each region 
variable_dict = {
    'cr_1m_j' : {
        'dataset_regex' : 'WJetsToLNu.*2017',
        'variables' : [
            'recoil',
            'muon_eta',
            'muon_pt',
            'muon_phi'
        ]
    },
    'cr_1e_j' : {
        'dataset_regex' : 'WJetsToLNu.*2017',
        'variables' : [
            'recoil',
            'electron_eta',
            'electron_pt',
            'electron_phi'
        ]
    },
    'cr_2m_j' : {
        'dataset_regex' : 'DYJetsToLL.*400to600.*2017',
        'variables' : [
            'recoil',
            'muon_eta0', 'muon_pt0', 'muon_phi0',
            'muon_eta1', 'muon_pt1', 'muon_phi1',
            'dimuon_pt', 'dimuon_eta', 'dimuon_mass'
        ]
    },
    'cr_2e_j' : {
        'dataset_regex' : 'DYJetsToLL.*400to600.*2017',
        'variables' : [
            'recoil',
            'electron_eta0', 'electron_pt0', 'electron_phi0',
            'electron_eta1', 'electron_pt1', 'electron_phi1',
            'dielectron_pt', 'dielectron_eta', 'dielectron_mass'
        ]
    }

}

def save_to_tree(acc, dataset_regex, region, variables, outputrootfile):
    '''Save variables for selected dataset (in the selected region) to a ROOT tree.'''
    
    # Load variables from cache
    for var in variables:
        acc.load(var)
        h = acc[var]

        # Merge extensions + rescale 
        h = merge_extensions(h, acc, reweight_pu=False)
        scale_xs_lumi(h)

        # Pick relevant dataset and region
        h = h.integrate('region', region).integrate('dataset', re.compile(dataset_regex))

        # Export to 1D ROOT histogram
        outputrootfile[var] = hist.export1d(h)
        print(f'Variable saved: {var}, region: {region}')

def main():
    inpath = sys.argv[1]

    if inpath.endswith('/'):
        outtag = inpath.split('/')[-2]
    else:
        outtag = inpath.split('/')[-1]
    
    acc = dir_archive(
        inpath,
        serialized=True,
        compression=0,
        memsize=1e3
    )
    acc.load('sumw')
    acc.load('sumw2')

    # Loop over the regions
    for region, info in variable_dict.items():
        # Create output ROOT file to save the histograms
        outdir = f'./output/{outtag}/comparison_with_uw'
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        
        outfile = pjoin(outdir, f'{region}_values.root')
        outputrootfile = uproot.recreate(outfile)

        dataset_regex, variables = info.values()

        save_to_tree(acc, dataset_regex=dataset_regex, region=region, variables=variables, outputrootfile=outputrootfile)

if __name__ == '__main__':
    main()