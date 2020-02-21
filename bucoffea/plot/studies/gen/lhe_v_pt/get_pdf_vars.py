#!/usr/bin/env python

import os
import sys
import re
import uproot
import numpy as np

from bucoffea.plot.util import merge_datasets, merge_extensions, scale_xs_lumi
from bucoffea.helpers.paths import bucoffea_path
from coffea import hist
from klepto.archives import dir_archive
from matplotlib import pyplot as plt

pjoin = os.path.join

def plot_variations(nom, var, tag):
    '''For the nominal weights and each set of
       variational weights, plot the ratio nom/var
       as a histogram.'''
    var_over_nom = []
    for variation in var.values():
        for idx in range(len(variation)):
            ratio = variation[idx]/nom[idx]
            var_over_nom.append(ratio)
    
    # Plot the results
    fig, ax = plt.subplots(1,1)
    tag_to_title = {
        'dy'    : r'PDF variations: $Z\rightarrow \ell \ell$',
        'wjet'  : r'PDF variations: $W\rightarrow \ell \nu$',
        'gjets' : r'PDF variations: $\gamma$ + jets'
    }
    if tag == 'd':
        bins = np.linspace(0.2, 1.8)
    else:
        bins = np.linspace(0.9, 1.1, 20)
    ax.hist(var_over_nom, bins=bins)
    ax.set_xlabel('Var / Nom')
    ax.set_ylabel('Counts')
    title = tag_to_title[tag]
    ax.set_title(title)

    # Save figure
    outdir = './output/theory_variations/pdf'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outpath = pjoin(outdir, f'{tag}_var_nom_weights.pdf')
    fig.savefig(outpath)
    print(f'Histogram saved in: {outpath}')

def hessian_unc(nom, var):
    '''Calculate PDF uncertainty for a Hessian set.'''
    unc=np.zeros_like(nom) 
    for idx,variation in enumerate(var.values()):
        unc += (nom - variation)**2
    return np.sqrt(unc)

def mc_unc(nom, var):
    '''Calculate PDF uncertainty for a MC set.'''
    # Calculate the average of all variations
    var_avg = np.zeros_like(nom)
    for variation in var.values():
        var_avg += variation
    var_avg /= len(var)
    # Calculate the MC uncertainty
    unc = np.zeros_like(nom)
    for variation in var.values():
        unc += (variation-var_avg)**2
    return np.sqrt(unc/(len(var)-1))

def calculate_pdf_unc(nom, var, tag):
    '''Given the nominal and varied weight content,
       calculate the PDF uncertainty.'''
    # Use PDF uncertainties for Hessian sets 
    # if samples is a DY or W sample
    if tag in ['wjet', 'dy']:
        unc = hessian_unc(nom, var)
    elif tag == 'gjets':
        unc = mc_unc(nom, var) 
    # Return uncertainty and percent uncertainty
    return unc, unc/nom 

def get_pdf_uncertainty(acc, regex, tag, outputrootfile, nominal='pdf_0'):
    '''Given the input accumulator, calculate the
       PDF uncertainty from all PDF variations.'''
    # Define rebinning
    vpt_ax_fine = list(range(0,400,40)) + list(range(400,1200,80))
    if tag in ['wjet', 'dy']:
        vpt_ax = hist.Bin('vpt','V $p_{T}$ (GeV)', vpt_ax_fine)
        mjj_ax = hist.Bin('mjj','M(jj) (GeV)',[0,200]+list(range(500,2500,500)))
    elif tag in ['gjets']:
        vpt_ax = hist.Bin('vpt','V $p_{T}$ (GeV)', vpt_ax_fine)
        mjj_ax = hist.Bin('mjj','M(jj) (GeV)',[0,200,500,1000,1500,2000])

    # Set the correct pt type
    pt_tag = 'combined' if tag != 'gjets' else 'stat1'
    acc.load(f'gen_vpt_vbf_{pt_tag}')
    h = acc[f'gen_vpt_vbf_{pt_tag}']

    h = h.rebin('vpt', vpt_ax)

    h = merge_extensions(h, acc, reweight_pu=False)
    scale_xs_lumi(h)
    h = merge_datasets(h)
    h = h[re.compile(regex)]

    # Integrate out mjj to get 1D variations 
    # as a function of V-pt
    mjj_slice = slice(200,7500)
    h = h.integrate('mjj', mjj_slice, overflow='over')

    # Get NLO distribution
    nlo = h[re.compile('.*(LHE|amcat).*')].integrate('dataset')

    # Nominal NLO weights, as specified in arguments
    # By defualt, use first PDF variation as nominal
    nlo_nom = nlo.integrate('var', nominal).values(overflow='over')[()]

    # NLO with PDF variations
    # Use a dict to collect NLO contents with all PDF variations
    nlo_var = {}

    for var in nlo.identifiers('var'):
        var_name = var.name
        if 'pdf' not in var_name: 
            continue
        nlo_var[var_name] = nlo.integrate('var', var_name).values(overflow='over')[()]

    unc, percent_unc = calculate_pdf_unc(nlo_nom, nlo_var, tag)
    print(percent_unc)

    plot_variations(nlo_nom, nlo_var, tag)

    # Plot the % uncertainty as a function of V-pt
    fig, ax = plt.subplots(1,1)
    vpt_edges = vpt_ax.edges(overflow='over')
    vpt_centers = ((vpt_edges + np.roll(vpt_edges, -1))/2)[:-1]
    ax.plot(vpt_centers, percent_unc, 'o')
    ax.set_xlabel(r'$p_T(V) \ (GeV)$')
    ax.set_ylabel(r'$\sigma_{pdf}$ / Nominal Counts')
    tag_to_title = {
        'dy'    : r'$Z\rightarrow \ell \ell$',
        'wjet'  : r'$W\rightarrow \ell \nu$',
        'gjets' : r'$\gamma$ + jets'
    }
    title = tag_to_title[tag]
    ax.set_title(title)
    ax.grid(True)
    ax.plot([200, 200], [0, 0.07], 'r')
    ax.set_ylim(0, 0.07)
           

    # Save the figure
    outdir = './output/theory_variations/pdf'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    outpath = pjoin(outdir, f'{tag}_pdf_unc.pdf')
    fig.savefig(outpath)

    # Save the uncertainties as a ROOT histogram
    outputrootfile[f'{tag}_pdf_unc'] = (percent_unc, vpt_ax.edges(overflow='over'))
    
    # Return nominal weights and uncertainty
    return nlo_nom, unc, vpt_edges, vpt_centers

def plot_ratio(noms, uncs, tag, vpt_edges, vpt_centers, outputrootfile):
    '''Plot the ratio of two processes for nominal case and two variations:
       Variation 1: (nom1+unc1)/(nom2+unc2)
       Variation 2: (nom1-unc1)/(nom2-unc2)
       Nominal: nom1/nom2
       '''
    nom1, nom2 = noms
    unc1, unc2 = uncs

    # Labels for y-axis
    tag_to_label = {
        'z_over_w' : r'$Z \rightarrow \ell \ell$ / $W \rightarrow \ell \nu$',
        'g_over_z' : r'$\gamma$ + jets / $Z \rightarrow \ell \ell$',
    }

    ratio_nom = nom1/nom2
    ratio_up = (nom1+unc1) / (nom2+unc2)
    ratio_down = (nom1-unc1) / (nom2-unc2)

    # Plot the ratios
    fig, ax = plt.subplots(1,1)
    ax.plot(vpt_centers, ratio_nom, 'o', label='Nominal')    
    ax.plot(vpt_centers, ratio_up, 'o', label='PDF up')    
    ax.plot(vpt_centers, ratio_down, 'o', label='PDF down')    
    ax.set_xlabel(r'$p_{T} (V)\ (GeV)$')
    ax.set_ylabel(tag_to_label[tag])
    ax.grid(True)
    ax.legend()

    ax.set_ylim(0,10)
    ax.plot([200,200], [0,10], 'r')

    # Save output
    outdir = './output/theory_variations/pdf/ratioplots'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    filepath = pjoin(outdir, f'{tag}_pdfunc_ratio.pdf')
    fig.savefig(filepath)

    ###########################
    # Plot the ratio of ratios!
    # Ratio(var) / Ratio(nom)
    ###########################
    dratio_up = ratio_up/ratio_nom
    dratio_down = ratio_down/ratio_nom
    
    # Plot the ratio of ratios
    plt.close('all')
    fig, ax = plt.subplots(1,1)
    ax.plot(vpt_centers, dratio_up, marker='o', label='PDF up / Nominal')
    ax.plot(vpt_centers, dratio_down, marker='o', label='PDF down / Nominal')
    
    ax.set_xlabel(r'$p_T (V)\ (GeV)$')
    ax.set_ylabel(f'{tag_to_label[tag]} (Var / Nom)')
    ax.legend()
    ax.grid(True)

    if tag == 'z_over_w':
        ax.set_ylim(0.99,1.01)
    elif tag == 'g_over_z':
        ax.set_ylim(0.95,1.05)
    
    ax.plot([200, 200], [ax.get_ylim()[0], ax.get_ylim()[1]], 'r')

    # Save output
    filepath = pjoin(outdir, f'{tag}_pdfunc_doubleratio.pdf')
    fig.savefig(filepath)

    # Save into root file
    outputrootfile[f'{tag}_var_over_nom_pdfup'] = (dratio_up, vpt_edges)
    outputrootfile[f'{tag}_var_over_nom_pdfdown'] = (dratio_down, vpt_edges)

def main():
    inpath = sys.argv[1]

    acc = dir_archive(
                       inpath,
                       serialized=True,
                       compression=0,
                       memsize=1e3
                     )

    acc.load('sumw')
    acc.load('sumw2')

    # Create the output ROOT file to save the 
    # PDF uncertainties as a function of v-pt
    outputrootfile = uproot.recreate('vbf_pdf_var.root')

    w_nom, w_unc, vpt_edges, vpt_centers = get_pdf_uncertainty(acc, regex='WNJetsToLNu.*', tag='wjet', outputrootfile=outputrootfile)
    dy_nom, dy_unc, vpt_edges, vpt_centers = get_pdf_uncertainty(acc, regex='DYNJetsToLL.*', tag='dy', outputrootfile=outputrootfile)
    gjets_nom, gjets_unc, vpt_edges, vpt_centers = get_pdf_uncertainty(acc, regex='G1Jet.*', tag='gjets', outputrootfile=outputrootfile)

    data_for_ratio = {
        'z_over_w' : {'noms' : (dy_nom, w_nom), 'uncs' : (dy_unc, w_unc)},
        'g_over_z' : {'noms' : (gjets_nom, dy_nom), 'uncs' : (gjets_unc, dy_unc)},
    }
    
    for tag, entry in data_for_ratio.items():
        noms = entry['noms']
        uncs = entry['uncs']
        plot_ratio(noms=noms, 
                   uncs=uncs,
                   tag=tag,
                   vpt_edges=vpt_edges,
                   vpt_centers=vpt_centers,
                   outputrootfile=outputrootfile)

if __name__ == '__main__':
    main()

