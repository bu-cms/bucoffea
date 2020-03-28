#!/usr/bin/env python

import os
import sys
import re
import numpy as np
from matplotlib import pyplot as plt
from coffea import hist
from bucoffea.plot.util import merge_extensions, merge_datasets, scale_xs_lumi
from klepto.archives import dir_archive
from pprint import pprint

pjoin = os.path.join

param_to_label = {
	'mjj' : r'Varied $M_{jj}$ / Nominal $M_{jj}$ - 1',
	'recoil' : r'Varied recoil / Nominal recoil - 1',
	'detajj' : r'Varied $\Delta\eta_{jj}$ / Nominal $\Delta\eta_{jj}$ - 1',
	'dphijj' : r'Varied $\Delta\phi_{jj}$ / Nominal $\Delta\phi_{jj}$ - 1',
}

var_to_label = {
	'_jesup' : 'JES up',
	'_jesdown' : 'JES down',
	'_jerup' : 'JER up',
	'_jerdown' : 'JER down',
}

def plot_comparison(acc, param, vars, dataregex, tag, outtag):
	'''
	Plot the distribution of the quantity:

	(Varied param / Nominal param) - 1
	
	for specified param and for data  
	and for all JES/JER variations listed in vars.

	Returns a dictionary containing:
	- Bin centers for the histogram
	- The distribution for all variations specified in vars list
	'''
	dist = f'{param}_varovernom'
	acc.load(dist)
	h = acc[dist]

	h = merge_extensions(h, acc, reweight_pu=False)
	scale_xs_lumi(h)
	h = merge_datasets(h)

	dataset_name = None

	# Find dataset name
	for dataset in h.identifiers('dataset'):
		if re.match(dataregex, dataset.name):
			dataset_name = dataset.name
		
	histo = h.integrate('dataset', re.compile(dataregex) )[re.compile('sr_vbf.*')]

	# Plot the distributions
	fig, (ax1, ax2) = plt.subplots(1,2,figsize=(13,6))

	centers = histo.axes()[1].centers()

	# Store the centers and distributions in a dict
	histdict = {
		'centers' : centers,
	}
	# Plot JES variations on left-hand side plot 
	# and JER variations on right-hand side plot
	for var in vars:
		histdict[f'sumw{var}'] = histo.integrate('region', f'sr_vbf{var}').values()[()]
		if 'jes' in var:
			ax1.step(centers, histdict[f'sumw{var}'], label=var_to_label[var])
		elif 'jer' in var:
			ax2.step(centers, histdict[f'sumw{var}'], label=var_to_label[var])

	for ax in [ax1, ax2]:
		ax.legend()
		ax.set_xlabel(param_to_label[param])
		ax.set_ylabel('Counts')
		ax.set_ylim(bottom=0)

		ylim = ax.get_ylim()
		ax.plot([0,0], ylim, 'r--')
		ax.set_ylim(ylim)

		ax.text(1,1,
				dataset_name,
				horizontalalignment='right',
				verticalalignment='bottom',
				fontsize=12,
				transform=ax.transAxes
				)
	
	# Save the figure
	outdir = f'./output/{outtag}/variation_comparisons'
	if not os.path.exists(outdir):
		os.makedirs(outdir)
	
	outpath = pjoin(outdir, f'{tag}_{param}.pdf')
	fig.savefig(outpath)

	print(f'File saved: {outpath}')

	return histdict

def plot_comparison_for_ratio(histinfos, param, var1, var2, data1, data2, tag, outtag):
	'''Given the dictionary containing dictionaries for each histogram (gotten from plot_comparison func),
	plots the distribution for the given parameter (param) for two variations (var1, var2) 
	for the ratio of two datasets.'''
	# Get the information related to two datasets 
	# that are being considered
	histinfo1 = histinfos[data1]
	histinfo2 = histinfos[data2]

	centers = histinfo1['centers']

	# Calculate the ratio of two datasets for each variation 
	# (variations are set in plot_comparison func)
	ratios = {
		var1 : histinfo1[f'sumw{var1}'] / histinfo2[f'sumw{var1}'],
		var2 : histinfo1[f'sumw{var2}'] / histinfo2[f'sumw{var2}']
	}

	# y-axis label according to tag
	tag_to_ylabel = {
		'zoverw17' : r'$Z\rightarrow \ell \ell$ / $W\rightarrow \ell \nu$ (2017)'
	}

	# Plot the ratios 
	fig, ax = plt.subplots(1,1)
	ax.step(centers, ratios[var1], label=var_to_label[var1])
	ax.step(centers, ratios[var2], label=var_to_label[var2])
	ax.legend()
	ax.set_xlabel(param_to_label[param])
	ax.set_ylabel(tag_to_ylabel[tag])
	
	# Save the plot
	outpath = f'./output/{outtag}/variation_comparisons/ratios'
	if not os.path.exists(outpath):
		os.makedirs(outpath)
	outfile = pjoin(outpath, f'{tag}_{param}{var1}{var2}.pdf')
	fig.savefig(outfile)

	print(f'File saved: {outfile}')

def main():
	inpath = sys.argv[1]
	acc = dir_archive(
			inpath,
			serialized=True,
			memsize=1e3,
			compression=0
			)
	
	acc.load('sumw')
	acc.load('sumw2')

	if inpath.endswith('/'):
		outtag = inpath.split('/')[-2]
	else:
		outtag = inpath.split('/')[-1]

	params = ['mjj', 'recoil']

	# Store dictionaries for each dataset in a dictionary
	histinfos = {}

	# List of JES/JER variations
	vars = ['_jesup', '_jesdown', '_jerup', '_jerdown']

	for param in params:
		# Create individual comparison plots for each dataset
		histinfos['ZJetsToNuNu_2017'] = plot_comparison(acc, param=param, 
										vars=vars,
										dataregex='ZJetsToNuNu.*2017', 
										tag='zjets', 
										outtag=outtag)

		histinfos['WJetsToLNu_2017'] = plot_comparison(acc, param=param, 
										vars=vars, 
										dataregex='WJetsToLNu.*2017', 
										tag='wjets', 
										outtag=outtag)
		
		# Create comparison plot for ratios of the two
		plot_comparison_for_ratio(histinfos, 
			param=param,
			var1='_jesup',
			var2='_jesdown',
			data1='ZJetsToNuNu_2017',
			data2='WJetsToLNu_2017',
			tag='zoverw17',
			outtag=outtag
			)

if __name__ == '__main__':
	main()

