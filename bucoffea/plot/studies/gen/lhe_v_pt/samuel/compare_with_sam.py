#!/usr/bin/env python

import os
import sys
import re
import uproot
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import ticker
from pprint import pprint

pjoin = os.path.join
plotcolors = plt.rcParams['axes.prop_cycle'].by_key()['color'][:4]

def get_sams_variations(sam_rootfile):
	'''Get scale variations from Samuel's root file.'''
	f = uproot.open(sam_rootfile)
	histdir = f['kfactors_shape']
	histnames = [histname.decode('utf-8') for histname in histdir]

	# Store the values of TH1 histograms in numpy arrays
	histvals = [histdir[histname].numpy() for histname in histnames]
	# Main container 
	histos = list(zip(histnames, histvals)) 

	return histos

def get_bu_variations(bu_rootfile):
	'''Get scale variations from BU root file.'''
	f = uproot.open(bu_rootfile)
	# Get the histograms in the same order as Sam's file
	order = ['renScaleDown', 'facScaleDown', 'facScaleUp', 'renScaleUp']
	histnames = []
	histvals = []
	for histtype in order:
		for histname in f:
			histname = histname.decode('utf-8')
			if histtype in histname:
				histnames.append(histname)
				histvals.append(f[histname].numpy())

	buhistos = list(zip(histnames, histvals))

	return buhistos

def get_combine_inputs(sam_rootfile, bu_rootfile, vary='all'):
	'''Get and combine the inputs from IC and BU for a specific variation.
	If vary = 'all' : The final container will have combined variations of the ratio.
	If vary = 'w' : The final container will have W variations of the ratio.
	If vary = 'z' : The final container will have Z variations of the ratio.
	'''
	f_sam = uproot.open(sam_rootfile)
	f_bu = uproot.open(bu_rootfile)

	sam_histdir = f_sam['kfactors_shape']

	sam_histnames = [histname.decode('utf-8') for histname in sam_histdir]
	bu_histnames = [histname.decode('utf-8') for histname in f_bu]

	# Get the combined variations if vary is set to "all"
	if vary == 'all':
		histnames_sam = [histname for histname in sam_histnames if 'Uncorr' in histname]
		histnames_bu = [histname for histname in bu_histnames if histname.startswith('w_over_z_')]
	# Get the W variations in ratio if vary is set to "w"
	elif vary == 'w':
		histnames_sam = [histname for histname in sam_histnames if 'Wup' in histname]
		histnames_bu = [histname for histname in bu_histnames if histname.startswith('wvar_over_z')]
	# Get the Z variations in ratio if vary is set to "z"
	elif vary == 'z':
		histnames_sam = [histname for histname in sam_histnames if 'Zup' in histname]
		histnames_bu = [histname for histname in bu_histnames if histname.startswith('w_over_zvar')]

	# Store the histogram values and bin edges as numpy arrays in the same order
	variations_sam = ['Renorm-Down', 'Renorm-Up', 'Fact-Down', 'Fact-Up']
	variations_bu = ['renScaleDown', 'renScaleUp', 'facScaleDown', 'facScaleUp']
	variations = zip(variations_sam, variations_bu)

	histvals_sam = []
	histvals_bu = []

	for var_sam, var_bu in variations:
		for histname in histnames_sam:
			if var_sam in histname:
				histvals_sam.append( sam_histdir[histname].numpy() )

		for histname in histnames_bu:
			if var_bu in histname:
				histvals_bu.append( f_bu[histname].numpy() )

	# Final object containing both histograms for all variations
	histvals_all = list(zip(histvals_bu, histvals_sam))

	return histvals_all

def plot_comparison(histvals_all, vartag):
	'''Given the container item containing two scale variations (one from IC, one from BU),
	   create a comparison plot between the two. Save the plot in ./output directory.'''
	edges = histvals_all[0][1][1]
	# Get centers
	centers = ( (edges + np.roll(edges, -1))/2)[:-1]

	# Labels compatible with histogram ordering
	labels = [r'$\mu_R$ down', r'$\mu_R$ up', r'$\mu_F$ down', r'$\mu_F$ up']
	# Create the comparison plot
	fig, (ax, rax) = plt.subplots(2,1, figsize=(7,7), gridspec_kw={'height_ratios' : (3,1)}, sharex=True)
	for idx, histvals in enumerate(histvals_all):
		legend_label = labels[idx]
		bu_val, sam_val = histvals
		# Plot Sam's points
		ax.plot(centers, sam_val[0], label=f'{legend_label} (IC)', marker='*')
		# Plot BU points
		ax.plot(centers, bu_val[0], 'o', label=f'{legend_label} (BU)', c=plotcolors[idx])
		# Calculate and plot the ratio between the two
		ratio_opts = {
			'linestyle' : 'none',
			'marker' : '.',
			'markersize' : 10.
		}
		ratio = sam_val[0] / bu_val[0]
		rax.plot(centers, ratio, 'o', color=plotcolors[idx], **ratio_opts)

	ax.legend(ncol=2)
	ax.set_ylim(0.7,1.5)
	ax.set_xlim(200,1000)
	ax.set_ylabel(r'($\frac{W}{Z}$ varied) / ($\frac{W}{Z}$ nominal)')
	if vartag == 'allvar':
		title = 'Uncorrelated Scale Variations'
		rax.set_ylim(0.96, 1.04)
		tickbase = 0.02
	elif vartag == 'wvar':
		title = r'$W \rightarrow \ell \nu$ Variations'
		rax.set_ylim(0.98, 1.02)
		tickbase = 0.01
	elif vartag == 'zvar':
		title = r'$Z \rightarrow \ell \ell$ Variations'
		rax.set_ylim(0.98, 1.02)
		tickbase = 0.01

	ax.set_title(title)

	loc = ticker.MultipleLocator(base=tickbase)

	rax.set_xlabel(r'$p_T (V)\ (GeV)$')
	rax.set_ylabel('IC / BU')
	rax.yaxis.set_major_locator(loc)
	rax.grid(True)

	# Save the output 
	outpath = f'./output'
	if not os.path.exists(outpath):
		os.makedirs(outpath)

	outfile = pjoin(outpath, f'scale_var_comp_{vartag}.pdf')
	fig.savefig(outfile)

	print(f'File saved: {outfile}')

def main():
	# Input root files from both groups
	sam_rootfile = './wz-ratio-uncertainty.root'
	bu_rootfile = './bu_woverz_scale_unc.root'
	
	histvals_all_wvar = get_combine_inputs(sam_rootfile, bu_rootfile, vary='w')
	histvals_all_zvar = get_combine_inputs(sam_rootfile, bu_rootfile, vary='z')
	histvals_all_combvar = get_combine_inputs(sam_rootfile, bu_rootfile, vary='all')

	plot_comparison(histvals_all_wvar, vartag='wvar')
	plot_comparison(histvals_all_zvar, vartag='zvar')
	plot_comparison(histvals_all_combvar, vartag='allvar')

if __name__ == '__main__':
	main()
