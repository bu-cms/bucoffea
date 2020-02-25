#!/usr/bin/env python

import os
import sys
import re
import numpy as np
from bucoffea.plot.util import merge_datasets, merge_extensions, scale_xs_lumi
from klepto.archives import dir_archive
from tabulate import tabulate
from pprint import pprint

pjoin = os.path.join

region_headers = {
	'sr_vbf' : 'Signal Region',
	'cr_1m_vbf' : 'Single Muon CR',
	'cr_2m_vbf' : 'Double Muon CR',
	'cr_1e_vbf' : 'Single Electron CR',
	'cr_2e_vbf' : 'Double Electron CR',
	'cr_g_vbf' : 'Photon CR'
}	

def dump_cutflows(acc, region, dataset, year):
	'''Dump cutflows for all variations in a given region and dataset.'''
	cutflow_dict = {}
	variations = ['', '_jesup', '_jesdown', '_jerup', '_jerdown']

	ncuts = None
	dataset_name = None

	for var in variations:
		tag = f'cutflow_{region}{var}'
		acc.load(tag)
		c = acc[tag]

		nevents = []

		for data in c.keys():
			if data.startswith(dataset) and year in data:
				cutflow = c[data]
				if not ncuts:
					ncuts = len(cutflow.values())
				if not dataset_name:
					dataset_name = data
				nevents.append( list(cutflow.values()) )
			
		sum_nevents = np.sum(nevents, axis=0)
		cutflow_dict[f'{region}{var}'] = sum_nevents

	# Calculate % differences w.r.t. nominal
	diffs = {}
	for var in variations: 
		if var == '': continue
		key = f'{region}{var}'
		abs_diff = np.abs( np.array( list(cutflow_dict[key]) ) - np.array( list(cutflow_dict[region]) ) )
		percent_diff = abs_diff*100 / np.array( list(cutflow_dict[region]) )
		diffs[f'{var.strip("_")}_nom'] = percent_diff

	# List of cut names
	cut_names = acc['cutflow_sr_vbf'][dataset_name].keys()

	headers = ['Nominal', 'JES up', 'JES down', 'JER up', 'JER down', 'JES Up-Nom (%)', 'JES Down-Nom (%)', 'JER Up-Nom (%)', 'JER Down-Nom (%)']
	
	# Merge the two dicts
	cutflow_dict.update(diffs)

	# Create table
	table = tabulate(cutflow_dict, headers=headers, showindex=cut_names, floatfmt=['.0f']*6 + ['.3f']*4, numalign='right')
	
	# Dump as output
	outdir = './output/cutflow_comparisons'
	if not os.path.exists(outdir):
		os.makedirs(outdir)
	
	outfile = pjoin(outdir, f'{dataset}_{region}_{year}_cutflow_comp.txt')
	with open(outfile, 'w+') as f:
		f.write(f'---{region_headers[region]}---\n\n')
		f.write(f'Dataset: {dataset}_{year}\n\n')
		f.write(table)

	print(f'File saved: {outfile}')

	return cutflow_dict, cut_names

def dump_cutflows_ratio(cutflows, datasets, year, region, tag, cut_names):
	'''Dump cutflow as a ratio of two datasets (e.g. Z/W)'''
	variations = ['', '_jesup', '_jesdown', '_jerup', '_jerdown']
	ratios = {}

	for var in variations:
		cutflow1 = cutflows[f'{datasets[0]}_{year}'][f'{region}{var}']
		cutflow2 = cutflows[f'{datasets[1]}_{year}'][f'{region}{var}']
		# Get the ratio of two cutflows
		ratio = cutflow1 / cutflow2
		ratios[var] = ratio
	
	# Calculate % differences w.r.t. nominal
	diffs = {}
	for var in variations: 
		if var == '': continue
		key = f'{region}{var}'
		abs_diff = np.abs(ratios[var] - ratios[''])
		percent_diff = abs_diff*100 / ratios['']
		diffs[f'{var.strip("_")}_nom'] = percent_diff
	
	# Merge the two dicts
	ratios.update(diffs)

	headers = ['Nominal', 'JES up', 'JES down', 'JER up', 'JER down', 'JES Up-Nom (%)', 'JES Down-Nom (%)', 'JER Up-Nom (%)', 'JER Down-Nom (%)']

	# Create table
	table = tabulate(ratios, headers=headers, showindex=cut_names, floatfmt='.3f', numalign='right')

	outdir='./output/cutflow_comparisons/ratios'
	if not os.path.exists(outdir):
		os.makedirs(outdir)
	
	outfile = pjoin(outdir, f'{datasets[0]}_{datasets[1]}_ratio_{year}_cutflow.txt')
	with open(outfile, 'w+') as f:
		f.write(f'---{region_headers[region]}---\n\n')
		f.write(f'Dataset 1: {datasets[0]}_{year}\n')
		f.write(f'Dataset 2: {datasets[1]}_{year}\n\n')
		f.write(table)

	print(f'File saved: {outfile}')
	

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

	datasets = ['ZJetsToNuNu', 'WJetsToLNu']
	years = ['2017', '2018']

	cutflows = {}

	for dataset in datasets:
		for year in years:
			cutflows[f'{dataset}_{year}'], cut_names = dump_cutflows(acc, region='sr_vbf', dataset=dataset, year=year)


	dump_cutflows_ratio(cutflows, datasets, year=2017, region='sr_vbf', tag='zoverw', cut_names=cut_names)


if __name__ == '__main__':
	main()
