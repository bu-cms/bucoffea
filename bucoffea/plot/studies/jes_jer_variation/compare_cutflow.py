#!/usr/bin/env python

import os
import sys
import re
import numpy as np
from klepto.archives import dir_archive
from tabulate import tabulate
from pprint import pprint
from bucoffea.plot.util import load_xs, lumi

pjoin = os.path.join

region_headers = {
	'sr_vbf' : 'Signal Region',
	'cr_1m_vbf' : 'Single Muon CR',
	'cr_2m_vbf' : 'Double Muon CR',
	'cr_1e_vbf' : 'Single Electron CR',
	'cr_2e_vbf' : 'Double Electron CR',
	'cr_g_vbf' : 'Photon CR'
}	

def get_weighted_cutflow(acc, dataset, region, var, year):
	'''For a specified scale variation (var), get the cutflow from 
	the input coffea files (acc).
	Cutflow for the given dataset, region and year will be considered.
	RETURNS:
	List containing the cut names and an array containing the cutflow.
	'''
	tag = f'cutflow_{region}{var}'
	acc.load(tag)
	c = acc[tag]
	r = re.compile(f'{dataset}.*{year}')
	# Get the relevant dataset names
	datasets = list(filter(r.match, c.keys()))

	# Get cut names
	cutnames = c[datasets[0]].keys()

	def get_cutflow(dataset):
		return np.array(list(c[dataset].values() ), dtype=float)
	
	def get_sumw(dataset):
		return acc['sumw'][dataset]

	mapped = {d: {'cutflow' : get_cutflow(d), 'sumw' : get_sumw(d)} for d in datasets} 

	# Merge extensions
	exts = [
		'.*(_ext\d+).*',
		'.*(_new_+pmx).*',
		'.*(_PSweights).*'
	]
	
	def merge_exts(_mapped, _exts):
		for dataset in datasets:
			for regex in _exts:
				m = re.match(regex, dataset)
				if m:
					dataset_new = dataset.replace(m.groups()[0], '')
					to_be_merged = _mapped.pop(dataset) 
					_mapped[dataset_new]['cutflow'] += to_be_merged['cutflow'] 
					_mapped[dataset_new]['sumw'] += to_be_merged['sumw'] 

		return _mapped

	merged_map = merge_exts(mapped, exts)

	# Reweight w.r.t. x-sec and lumi
	xs = load_xs()

	for d in merged_map.keys():
		reweight = 1e3*xs[d]*lumi(year) / merged_map[d]['sumw']
		merged_map[d]['cutflow'] *= reweight
	
	# Merge bins and get the final cutflow
	cutflow_list = [] 
	for d in merged_map.keys():
		cutflow_list.append(merged_map[d]['cutflow'])

	cutflow_arr = np.sum(cutflow_list, axis=0)

	return cutnames, cutflow_arr

def dump_cutflows(acc, region, dataset, year, out_tag):
	'''Dump cutflows for all variations in a given region and dataset.'''
	cutflow_dict = {}
	variations = ['', '_jesup', '_jesdown', '_jerup', '_jerdown']

	# Get cutflows for each variation
	for var in variations:
		cutnames, cutflow_dict[f'{region}{var}'] = get_weighted_cutflow(acc,
														   dataset=dataset,
														   region=region,
														   var=var,
														   year=year)

	# Calculate % differences w.r.t. nominal
	diffs = {}
	for var in variations: 
		if var == '': continue
		key = f'{region}{var}'
		abs_diff = np.abs( np.array( list(cutflow_dict[key]) ) - np.array( list(cutflow_dict[region]) ) )
		percent_diff = abs_diff*100 / np.array( list(cutflow_dict[region]) )
		diffs[f'{var.strip("_")}_nom'] = percent_diff

	headers = ['Nominal', 'JES up', 'JES down', 'JER up', 'JER down', 'JES Up-Nom (%)', 'JES Down-Nom (%)', 'JER Up-Nom (%)', 'JER Down-Nom (%)']
	
	# Merge the two dicts
	cutflow_dict.update(diffs)

	# Clean cut names 
	def clean_cutname(cutname):
		clean_regex = '.*_(jes|jer)(up|down)'
		m = re.match(clean_regex, cutname)
		if m:
			toreplace = '_' + ''.join(m.groups())
			cutname = cutname.replace(toreplace, '')

		return cutname

	cutnames = list(map(clean_cutname, cutnames))

	# Create table
	table = tabulate(cutflow_dict, headers=headers, showindex=cutnames, floatfmt=['.0f']*6 + ['.3f']*4, numalign='right')

	# Dump as output
	outdir = f'./output/{out_tag}/cutflow_comparisons'
	if not os.path.exists(outdir):
		os.makedirs(outdir)
	
	outfile = pjoin(outdir, f'{dataset}_{region}_{year}_cutflow_comp.txt')
	with open(outfile, 'w+') as f:
		f.write(f'---{region_headers[region]}---\n\n')
		f.write(f'Dataset: {dataset}_{year}\n\n')
		f.write(table)

	print(f'File saved: {outfile}')

	return cutflow_dict, cutnames

def dump_cutflows_ratio(cutflows, datasets, year, region, tag, cut_names, out_tag):
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

	outdir=f'./output/{out_tag}/cutflow_comparisons/ratios'
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

	if inpath.endswith('/'):
		out_tag = inpath.split('/')[-2]
	else:
		out_tag = inpath.split('/')[-1]

	datasets = ['ZJetsToNuNu', 'WJetsToLNu']
	years = [2017, 2018]

	cutflows = {}

	for dataset in datasets:
		for year in years:
			cutflows[f'{dataset}_{year}'], cut_names = dump_cutflows(acc, region='sr_vbf', dataset=dataset, year=year, out_tag=out_tag)

	dump_cutflows_ratio(cutflows, datasets, year=2017, region='sr_vbf', tag='zoverw', cut_names=cut_names, out_tag=out_tag)
	dump_cutflows_ratio(cutflows, datasets, year=2018, region='sr_vbf', tag='zoverw', cut_names=cut_names, out_tag=out_tag)


if __name__ == '__main__':
	main()
