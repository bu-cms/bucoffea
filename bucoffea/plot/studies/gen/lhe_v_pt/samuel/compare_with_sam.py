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
	f = uproot.open(sam_rootfile)
	histdir = f['kfactors_shape']
	histnames_all = [histname.decode('utf-8') for histname in histdir]
	# Take the correlated histograms only
	histnames = [histname for histname in histnames_all if 'Uncorr' in histname]
	
	# Store the values of TH1 histograms in numpy arrays
	histvals = [histdir[histname].numpy() for histname in histnames]
	# Main container 
	histos = list(zip(histnames, histvals)) 

	return histos

def get_bu_variations(bu_rootfile):
	f = uproot.open(bu_rootfile)
	# Get the histograms in the same order as Sam's file
	order = ['renScaleDown', 'facScaleDown', 'facScaleUp', 'renScaleUp']
	histnames = []
	histvals = []
	for histtype in order:
		for histname in f:
			if histtype in histname.decode('utf-8'):
				histnames.append(histname)
				histvals.append(f[histname].numpy())

	buhistos = list(zip(histnames, histvals))

	return buhistos

def plot_comparison(sam_var, bu_var):
	'''Given two containers containing two scale variations (one from IC, one from BU),
	   create a comparison plot between the two. Save the plot in ./output directory.'''
	edges = sam_var[0][1][1]
	# Get centers
	centers = ( (edges + np.roll(edges, -1))/2)[:-1]

	sam_histnames = [item[0] for item in sam_var] 
	sam_histvals = [item[1][0] for item in sam_var]

	bu_histnames = [item[0] for item in bu_var] 
	bu_histvals = [item[1][0] for item in bu_var]

	# Create the comparison plot
	fig, (ax, rax) = plt.subplots(2,1, figsize=(7,7), gridspec_kw={'height_ratios' : (3,1)}, sharex=True)
	for idx in range(len(sam_histvals)):
		labeltokens = re.findall('(Renorm|Fact)-(Up|Down)', sam_histnames[idx])[0]
		label = ' '.join(labeltokens)
		sam_val = sam_histvals[idx]
		bu_val = bu_histvals[idx]
		# Plot Sam's points
		ax.plot(centers, sam_val, label=f'{label} (IC)', marker='*')
		# Plot BU points
		ax.plot(centers, bu_val, 'o', label=f'{label} (BU)', c=plotcolors[idx])
		# Calculate and plot the ratio between the two
		ratio_opts = {
			'linestyle' : 'none',
			'marker' : '.',
			'markersize' : 10.
		}
		ratio = sam_val / bu_val
		rax.plot(centers, ratio, 'o', color=plotcolors[idx], **ratio_opts)

	ax.legend(ncol=2)
	ax.set_ylim(0.7,1.5)
	ax.set_xlim(200,1000)
	ax.set_ylabel(r'($\frac{W}{Z}$ varied) / ($\frac{W}{Z}$ nominal)')
	ax.set_title('Uncorrelated Scale Variations')

	loc = ticker.MultipleLocator(base=0.02)

	rax.set_xlabel(r'$p_T (V)\ (GeV)$')
	rax.set_ylabel('IC / BU')
	rax.set_ylim(0.96, 1.04)
	rax.yaxis.set_major_locator(loc)
	rax.grid(True)

	# Save the output 
	outpath = f'./output'
	if not os.path.exists(outpath):
		os.makedirs(outpath)
	outfile = pjoin(outpath, f'scale_var_comp.pdf')
	fig.savefig(outfile)


def main():
	sam_rootfile = './wz-ratio-uncertainty.root'
	sam_var = get_sams_variations(sam_rootfile)

	bu_rootfile = './bu_woverz_scale_unc.root'
	bu_var = get_bu_variations(bu_rootfile)

	plot_comparison(sam_var, bu_var)

if __name__ == '__main__':
	main()
