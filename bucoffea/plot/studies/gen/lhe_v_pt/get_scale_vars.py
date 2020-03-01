#!/usr/bin/env python

import os
import sys
import re
import uproot
import numpy as np
import pickle
from bucoffea.plot.util import merge_datasets, merge_extensions, scale_xs_lumi
from bucoffea.helpers.paths import bucoffea_path
from coffea import hist
from klepto.archives import dir_archive
from matplotlib import pyplot as plt

pjoin = os.path.join

def get_old_kfac(tag):
    '''Given the dataset tag, get the nominal 2D VBF k-factors.'''
    f = uproot.open(bucoffea_path('data/sf/theory/2017_gen_v_pt_qcd_sf.root'))
    return f[f'2d_{tag}_vbf'].values

def get_scale_variations(acc, regex, tag, scale_var):
    '''With the given scale variation, get the ratio between the varied k-factors
    and nominal k-factors, for the specified dataset. 
    Returns two tuples:
    1. A tuple containing this ratio as a function of V-pt,
    and an array containing V-pt histogram edges. 
    2. A tuple containing sum of NLO weights for varied case and nominal case.
    ==============
    PARAMETERS:
    ==============
    acc   : The coffea accumulator containing all the input histograms.
    regex : The regular expression matching the dataset name to be considered.
    tag   : The tag for the dataset being used.
    scale_var      :  Tag for which scale variation is going to be considered. 
                      Used to get the relevant scale variated inputs from the coffea file.
                      These are in the from of "scale_idx" (e.g. scale_01)
    '''

    print(f'Working on: {tag}, {scale_var}')

    # Define rebinning
    if tag in ['wjet', 'dy']:
        vpt_ax_coarse = [0, 40, 80, 120, 160, 200, 240, 280, 320, 400, 520, 640, 760, 880, 1100]
        vpt_ax_fine = list(range(0,400,40)) + list(range(400,1200,80))
        vpt_ax = hist.Bin('vpt','V $p_{T}$ (GeV)', vpt_ax_coarse)
        mjj_ax = hist.Bin('mjj','M(jj) (GeV)', [0,200] + list(range(500,2500,500)))
    elif tag in ['gjets']:
        vpt_ax_coarse = [0, 40, 80, 120, 160, 200, 240, 280, 320, 400, 520, 640, 760, 880, 1100]
        vpt_ax_fine = list(range(0,400,40)) + list(range(400,1200,80)) 
        vpt_ax = hist.Bin('vpt','V $p_{T}$ (GeV)', vpt_ax_coarse)
        mjj_ax = hist.Bin('mjj','M(jj) (GeV)',[0,200,500,1000,1500])

    # Get the correct pt type from coffea input
    pt_tag = 'combined' if tag != 'gjets' else 'stat1'
    acc.load(f'gen_vpt_vbf_{pt_tag}')
    h = acc[f'gen_vpt_vbf_{pt_tag}']

    # Rebin
    h = h.rebin('vpt', vpt_ax)
    h = h.rebin('mjj', mjj_ax)

    # Merging extensions/datasets, scaling w.r.t xs and lumi 
    h = merge_extensions(h, acc, reweight_pu=False)
    scale_xs_lumi(h)
    h = merge_datasets(h)
    h = h[re.compile(regex)]

    # Get LO and NLO inputs, to calculate the scale factors later
    lo = h[re.compile('.*HT.*')].integrate('dataset')
    nlo = h[re.compile('.*(LHE|amcat).*')].integrate('dataset')
    
    # Print choose the relevant scale variation (relevant to NLO only)
    # For LO, choose the nominal (i.e. no variation)
    lo = lo.integrate('var', 'nominal')
    nlo_var = nlo.integrate('var', scale_var)
    nlo_nom = nlo.integrate('var', 'nominal')

    # Get 1D LO and NLO weights to calculate the variation
    if tag in ['wjet', 'dy']:
        mjj_slice = slice(200,2000)
    elif tag in ['gjets']:
        mjj_slice = slice(200,1500)
    lo_1d = lo.integrate('mjj', mjj_slice)
    nlo_var_1d = nlo_var.integrate('mjj', mjj_slice)
    nlo_nom_1d = nlo_nom.integrate('mjj', mjj_slice)
    
    sumw_lo_1d = lo_1d.values()[()]
    sumw_nlo_var_1d = nlo_var_1d.values()[()]
    sumw_nlo_nom_1d = nlo_nom_1d.values()[()]

    # Calculate 1D scale factors, nominal and varied
    # as a function of V-pt
    sf_nom_1d = sumw_nlo_nom_1d / sumw_lo_1d
    sf_var_1d = sumw_nlo_var_1d / sumw_lo_1d

    # Calculate 1D variation ratio, as a function of V-pt
    var_ratio = sf_var_1d / sf_nom_1d
    
    tup1 = (var_ratio, h.axis('vpt').edges() )
    tup2 = (sumw_nlo_var_1d, sumw_nlo_nom_1d)
    # Return tuple containing the SF ratios and
    # NLO weights with and without variation
    return tup1, tup2

def plot_ratio_vpt(tup, var, tag, outtag):
    '''Given the tuple contatining the SF ratio (variational/nominal) and 
    bin edges, plot ratios in each mjj bin as a function of v-pt.
    Plot is made ONLY for one scale variation.'''
    ratio, x_edges = tup
    vpt_centers = ((x_edges + np.roll(x_edges,-1))/2)[:-1]
    fig, ax = plt.subplots(1,1)
    
    # Figure out the variation and the relevant title
    var_title = {
        'gjets' : {
            'scale_1' : r'$\gamma$ + jets: $\mu_R = 0.5$, $\mu_F = 1.0$',
            'scale_3' : r'$\gamma$ + jets: $\mu_R = 1.0$, $\mu_F = 0.5$',
            'scale_5' : r'$\gamma$ + jets: $\mu_R = 1.0$, $\mu_F = 2.0$',
            'scale_7' : r'$\gamma$ + jets: $\mu_R = 2.0$, $\mu_F = 1.0$'
        },
        'wjet' : {
            'scale_1' : r'$W\rightarrow l \nu$: $\mu_R = 0.5$, $\mu_F = 1.0$',
            'scale_3' : r'$W\rightarrow l \nu$: $\mu_R = 1.0$, $\mu_F = 0.5$',
            'scale_4' : r'$W\rightarrow l \nu$: $\mu_R = 1.0$, $\mu_F = 2.0$',
            'scale_6' : r'$W\rightarrow l \nu$: $\mu_R = 2.0$, $\mu_F = 1.0$'
        },
        'dy' : {
            'scale_1' : r'$Z\rightarrow ll$: $\mu_R = 0.5$, $\mu_F = 1.0$',
            'scale_3' : r'$Z\rightarrow ll$: $\mu_R = 1.0$, $\mu_F = 0.5$',
            'scale_4' : r'$Z\rightarrow ll$: $\mu_R = 1.0$, $\mu_F = 2.0$',
            'scale_6' : r'$Z\rightarrow ll$: $\mu_R = 2.0$, $\mu_F = 1.0$'
        },
    }

    fig_title = var_title[tag][var] 

    ax.set_xlim(-50.0, vpt_centers[-1]+50.0)
    ax.set_ylim(0.8, 1.2)
    ax.set_xlabel(r'$p_T(V) \ (GeV)$')
    ax.set_ylabel('Varied / Nominal SF')
    ax.set_title(fig_title)
    ax.plot([-50.0, vpt_centers[-1]+50.0], [1, 1], 'r')
    ax.plot(vpt_centers, ratio, marker='o')
    ax.grid(True)

    # Save the figure
    outpath = f'./output/theory_variations/{outtag}/scale/individual'
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    outfile = pjoin(outpath, f'{tag}_kfac_ratio_{var}.pdf')
    fig.savefig(outfile)

def plot_ratio_vpt_combined(tup_combined, tag, outtag):
    '''Given the tup_combined contatining the SF ratio (variational/nominal) and 
    bin edges for all variations for a given process, 
    plot ratios in each mjj bin as a function of v-pt.
    
    Plot is made for ALL scale variations that are present in tup_combined.'''
    ratio_combined, x_edges = tup_combined[:,0], tup_combined[0,1]
    vpt_centers = ((x_edges + np.roll(x_edges,-1))/2)[:-1]
    fig, ax = plt.subplots(1,1)

    # Figure out the variation and the relevant title, and label for the legend
    var_title_label = {
        'gjets' : {
            'scale_1' : (r'$\gamma$ + jets', r'$\mu_R = 0.5$, $\mu_F = 1.0$'),
            'scale_3' : (r'$\gamma$ + jets', r'$\mu_R = 1.0$, $\mu_F = 0.5$'),
            'scale_5' : (r'$\gamma$ + jets', r'$\mu_R = 1.0$, $\mu_F = 2.0$'),
            'scale_7' : (r'$\gamma$ + jets', r'$\mu_R = 2.0$, $\mu_F = 1.0$')
        },
        'wjet' : {
            'scale_1' : (r'$W\rightarrow \ell \nu$', r'$\mu_R = 0.5$, $\mu_F = 1.0$'),
            'scale_3' : (r'$W\rightarrow \ell \nu$', r'$\mu_R = 1.0$, $\mu_F = 0.5$'),
            'scale_4' : (r'$W\rightarrow \ell \nu$', r'$\mu_R = 1.0$, $\mu_F = 2.0$'),
            'scale_6' : (r'$W\rightarrow \ell \nu$', r'$\mu_R = 2.0$, $\mu_F = 1.0$')
        },
        'dy' : {
            'scale_1' : (r'$Z\rightarrow \ell \ell$', r'$\mu_R = 0.5$, $\mu_F = 1.0$'),
            'scale_3' : (r'$Z\rightarrow \ell \ell$', r'$\mu_R = 1.0$, $\mu_F = 0.5$'),
            'scale_4' : (r'$Z\rightarrow \ell \ell$', r'$\mu_R = 1.0$, $\mu_F = 2.0$'),
            'scale_6' : (r'$Z\rightarrow \ell \ell$', r'$\mu_R = 2.0$, $\mu_F = 1.0$')
        },
    }

    variations = var_title_label[tag]
    for idx, var_info in enumerate(variations.values()):
        ax.plot([-50.0, vpt_centers[-1]+50.0], [1, 1], 'r')
        ax.plot(vpt_centers, ratio_combined[idx], marker='o', label=var_info[1])

    ax.set_xlim(-50.0, vpt_centers[-1]+50.0)
    ax.set_ylim(0.8, 1.2)
    ax.set_xlabel(r'$p_T(V) \ (GeV)$')
    ax.set_ylabel('Varied / Nominal SF')
    fig_title = variations['scale_1'][0]
    ax.set_title(fig_title)
    ax.grid(True)
    plt.legend()

    # Save the figure
    outpath = f'./output/theory_variations/{outtag}/scale/combinedplot'
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    outfile = pjoin(outpath, f'{tag}_kfac_ratio_combined.pdf')
    fig.savefig(outfile)

def get_ratios(sumw_var, tag, var):
    '''Get ratio for two physics processes, for a given scale variation.'''
    # Figure out the processes
    if tag in ['zvar_over_w', 'z_over_wvar']:
        tag1 = 'dy'
        tag2 = 'wjet'
    elif tag in ['wvar_over_z', 'w_over_zvar']:
        tag1 = 'wjet'
        tag2 = 'dy'
    elif tag in ['gvar_over_z', 'g_over_zvar']:
        tag1 = 'gjets'
        tag2 = 'dy'

    # Get varied and nominal weights for
    # the given scale variation
    sumw1_var, sumw1_nom = sumw_var[tag1][var]    
    sumw2_var, sumw2_nom = sumw_var[tag2][var]
    
    ratio_var1 = sumw1_var / sumw2_nom
    ratio_var2 = sumw1_nom / sumw2_var
    ratio_nom = sumw1_nom / sumw2_nom
    
    # Return the varied and nominal ratios
    return (ratio_var1, ratio_var2), ratio_nom 

def get_ratio_arbitrary(sumw_var, tag, var1, var2):
    '''Get ratio of two processes, such that:
       - Process 1 is varied with var1
       - Process 2 is varied with var2'''
    if re.match('z.*_over_w.*', tag):
        tag1 = 'dy'
        tag2 = 'wjet'
    elif re.match('g.*_over_z.*', tag):
        tag1 = 'gjets'
        tag2 = 'dy'

    sumw1_var, sumw1_nom = sumw_var[tag1][var1]
    sumw2_var, sumw2_nom = sumw_var[tag2][var2]

    ratio_var = sumw1_var / sumw2_var
    ratio_nom = sumw1_nom / sumw2_nom

    return ratio_var, ratio_nom
    
def plot_ratio(sumw_var, tag, xedges, outtag, varied='num'):
    '''Plot ratio for two processes, for all variations.
       Specify which process is to be varied (num or denom).'''
    # List of all variations
    varlist = sumw_var['wjet'].keys()

    xcenters = ((xedges + np.roll(xedges,-1))/2)[:-1]

    tag_to_ylabel = {
        'zvar_over_w' : r'$Z\rightarrow \ell \ell$ / $W\rightarrow \ell \nu$ (Var / Nom)',
        'z_over_wvar' : r'$Z\rightarrow \ell \ell$ / $W\rightarrow \ell \nu$ (Var / Nom)',
        'wvar_over_z' : r'$W\rightarrow \ell \nu$ / $Z\rightarrow \ell \ell$ (Var / Nom)',
        'w_over_zvar' : r'$W\rightarrow \ell \nu$ / $Z\rightarrow \ell \ell$ (Var / Nom)',
        'gvar_over_z' : r'$\gamma$ + jets / $Z\rightarrow \ell \ell$ (Var / Nom)',
        'g_over_zvar' : r'$\gamma$ + jets / $Z\rightarrow \ell \ell$ (Var / Nom)',
    }

    var_to_label = {
        'mu_r_down' : r'$\mu_R = 0.5$, $\mu_F = 1.0$',
        'mu_r_up' : r'$\mu_R = 2.0$, $\mu_F = 1.0$',
        'mu_f_down' : r'$\mu_R = 1.0$, $\mu_F = 0.5$',
        'mu_f_up' : r'$\mu_R = 1.0$, $\mu_F = 2.0$',
    }

    # Construct the plot
    fig, ax = plt.subplots(1,1)

    # Store double ratios in a dict
    dratios = {}

    for var in varlist:
        ratios_var, ratio_nom = get_ratios(sumw_var, tag, var)
        # Get the ratio with relevant variations
        if varied == 'num':
            ratio_var = ratios_var[0]
        elif varied == 'denom':
            ratio_var = ratios_var[1]

        # Calculate the ratio of ratios!
        dratio = ratio_var / ratio_nom

        ax.plot(xcenters, dratio, marker='o', label=var_to_label[var])
        
        dratios[var] = dratio
    
    ax.set_xlabel(r'$p_T (V)\ (GeV)$')
    ax.set_ylabel(tag_to_ylabel[tag])
    ax.set_ylim(0.9, 1.1)
    if 'zvar' in tag:
        legend_title = r'$Z\rightarrow \ell \ell$ variations'
    elif 'wvar' in tag:
        legend_title = r'$W\rightarrow \ell \nu$ variations'
    elif 'gvar' in tag:
        legend_title = r'$\gamma$ + jets variations'
    ax.legend(title=legend_title)
    ax.grid(True)

    ax.plot([200,200], [0.9,1.1], 'k')

    # Save figure
    outdir = f'./output/theory_variations/{outtag}/scale/ratioplots'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    outpath = pjoin(outdir, f'{tag}_scaleunc_ratio.pdf')
    fig.savefig(outpath)

    print(f'Saved file in {outpath}')

    return dratios

def plot_ratio_anticorr(sumw_var, tag, xedges):
    '''Plot ratio for two processes, for all anti-correlated variations.'''
    xcenters = ((xedges + np.roll(xedges,-1))/2)[:-1]
    
    tag_to_ylabel = {
        'z_over_w' : r'$Z\rightarrow \ell \ell$ / $W\rightarrow \ell \nu$',
        'g_over_z' : r'$\gamma$ + jets / $Z\rightarrow \ell \ell$',
    }

    var_to_label = {
        ('mu_r_up', 'mu_r_down') : r'$\mu_R$ up (Z), $\mu_R$ down (W)',
        ('mu_r_down', 'mu_r_up') : r'$\mu_R$ down (Z), $\mu_R$ up (W)',
        ('mu_f_up', 'mu_f_down') : r'$\mu_F$ up (Z), $\mu_F$ down (W)',
        ('mu_f_down', 'mu_f_up') : r'$\mu_F$ down (Z), $\mu_F$ up (W)'
    }

    # Construct the plot
    fig, ax = plt.subplots(1,1)

    # Store double ratios in a dict
    dratios = {}

    for var1, var2 in var_to_label.keys():
        ratio_var, ratio_nom = get_ratio_arbitrary(sumw_var, tag, var1, var2)

        # Calculate the ratio of ratios!
        dratio = ratio_var / ratio_nom

        ax.plot(xcenters, dratio, marker='o', label=var_to_label[(var1, var2)])
        
        dratios[(var1, var2)] = dratio
    
    ax.set_xlabel(r'$p_T (V)\ (GeV)$')
    ax.set_ylabel(tag_to_ylabel[tag])
    ax.set_ylim(0.8, 1.2)
    ax.legend()
    ax.grid(True)

    ax.plot([200,200], [0.8,1.2], 'k')

    # Save figure
    outdir = './output/theory_variations/scale/ratioplots/anticorr'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    outpath = pjoin(outdir, f'{tag}_scaleunc_ratio.pdf')
    fig.savefig(outpath)

    print(f'Saved file in {outpath}')

    return dratios

def plot_combined_scale_uncs(dratio_dict, xedges, outputrootfiles, outtag):
    '''Plot combined scale variations for Z/W and photon/W ratios.
       Combined var: (var_num**2 + var_denom**2)**0.5'''
    to_be_combined = [('zvar_over_w', 'z_over_wvar'), ('wvar_over_z', 'w_over_zvar'), ('gvar_over_z', 'g_over_zvar')]
    tags = ['z_over_w', 'w_over_z', 'g_over_z']
    var_pairs = [
        ('mu_r_down', 'mu_r_up'),
        ('mu_r_up', 'mu_r_down'),
        ('mu_f_down', 'mu_f_up'),
        ('mu_f_up', 'mu_f_down')
    ]
    var_to_label = {
        'mu_r_down' : r'$\mu_{R,1}$',
        'mu_r_up' : r'$\mu_{R,2}$',
        'mu_f_down' : r'$\mu_{F,1}$',
        'mu_f_up' : r'$\mu_{F,2}$',
    }
    var_to_roothistname = {
        'mu_r_down' : 'renScaleDown',
        'mu_r_up' : 'renScaleUp',
        'mu_f_down' : 'facScaleDown',
        'mu_f_up' : 'facScaleUp',
    }
    tag_to_ylabel = {
        'z_over_w' : r'Combined Scale Unc: $Z\rightarrow \ell \ell$ / $W\rightarrow \ell \nu$',
        'w_over_z' : r'Combined Scale Unc: $W\rightarrow \ell \nu$ / $Z\rightarrow \ell \ell$',
        'g_over_z' : r'Combined Scale Unc: $\gamma$ + jets / $Z\rightarrow \ell \ell$',
    }
    xcenters = ((xedges + np.roll(xedges,-1))/2)[:-1]
    for idx, pair in enumerate(to_be_combined):
        dratios_num = dratio_dict[pair[0]]
        dratios_denom = dratio_dict[pair[1]]
        fig, ax = plt.subplots(1,1)
        # Combine the variations
        for var_pair in var_pairs:
            combined_dratio_unc = 1 + np.sign(1-dratios_num[var_pair[0]]) * np.sqrt( (1-dratios_num[var_pair[0]] )**2 + ( 1-dratios_denom[var_pair[1]] )**2)
            ax.plot(xcenters, combined_dratio_unc, marker='o', label=var_to_label[var_pair[0]])
            # Save to output ROOT file
            tup = (combined_dratio_unc, xedges)
            outputrootfile = outputrootfiles[tags[idx]]
            outputrootfile[f'{tags[idx]}_{var_to_roothistname[var_pair[0]]}'] = tup
        ax.legend()
        ax.set_xlabel(r'$p_T(V)\ (GeV)$')
        ax.set_ylabel(tag_to_ylabel[tags[idx]])
        ax.grid(True)
        ax.set_xlim(-50,1200)
        ax.set_ylim(0.8,1.2)
        ax.plot([200,200], [0.8,1.2], 'k')
        ax.plot([-50,1200], [1,1], 'k--')
        outpath = f'./output/theory_variations/{outtag}/scale/ratioplots/combinedunc'
        if not os.path.exists(outpath):
            os.makedirs(outpath)
        outfile = pjoin(outpath, f'{tags[idx]}_combined_scaleuncs.pdf')
        fig.savefig(outfile)


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

    if inpath.endswith('/'):
        outtag = inpath.split('/')[-2]
    else:
        outtag = inpath.split('/')[-1]
        
    # Create the output ROOT file to save the 
    # PDF uncertainties as a function of v-pt
    outputrootpath = f'./output/theory_variations/{outtag}/rootfiles'
    if not os.path.exists(outputrootpath):
        os.makedirs(outputrootpath)
    
    outputrootfile_z_over_w = uproot.recreate( pjoin(outputrootpath, 'zoverw_scale_unc.root') )
    outputrootfile_w_over_z = uproot.recreate( pjoin(outputrootpath, 'woverz_scale_unc.root') )
    outputrootfile_g_over_z = uproot.recreate( pjoin(outputrootpath, 'goverz_scale_unc.root') )

    outputrootfiles = {
        'z_over_w' : outputrootfile_z_over_w,
        'w_over_z' : outputrootfile_w_over_z,
        'g_over_z' : outputrootfile_g_over_z
    }

    scale_var_dict = {
        'gjets' : [
            ('scale_1', 'mu_r_down'),
            ('scale_3', 'mu_f_down'),
            ('scale_5', 'mu_f_up'),
            ('scale_7', 'mu_r_up')
                ],
        'wjet/dy' : [
            ('scale_1', 'mu_r_down'),
            ('scale_3', 'mu_f_down'),
            ('scale_4', 'mu_f_up'),
            ('scale_6', 'mu_r_up')
                ]
                
    }

    tag_regex = {
        'wjet'  : 'WN?JetsToLNu.*',
        'dy'    : 'DYN?JetsToLL.*',
        'gjets' : 'G\d?Jet.*' 
    }

    pkl_dir = f'./output/theory_variations/{outtag}/scale'
    if not os.path.exists(pkl_dir):
        os.makedirs(pkl_dir)

    pkl_file = pjoin(pkl_dir, 'vbf_scale_sumwvar.pkl')

    # Store double ratios in a dict
    dratio_dict = {}

    # Load the results from pkl file if it exists
    # Otherwise, run the code to get scale variations
    if os.path.exists(pkl_file): 
        with open(pkl_file, 'rb') as f:
            tup = pickle.load(f)
            sumw_var = pickle.load(f)

        xedges = tup[1]
        dratio_dict['zvar_over_w'] = plot_ratio(sumw_var, tag='zvar_over_w', xedges=xedges, varied='num', outtag=outtag)
        dratio_dict['z_over_wvar'] = plot_ratio(sumw_var, tag='z_over_wvar', xedges=xedges, varied='denom', outtag=outtag)
        dratio_dict['wvar_over_z'] = plot_ratio(sumw_var, tag='wvar_over_z', xedges=xedges, varied='num', outtag=outtag)
        dratio_dict['w_over_zvar'] = plot_ratio(sumw_var, tag='w_over_zvar', xedges=xedges, varied='denom', outtag=outtag)
        dratio_dict['gvar_over_z'] = plot_ratio(sumw_var, tag='gvar_over_z', xedges=xedges, varied='num', outtag=outtag)
        dratio_dict['g_over_zvar'] = plot_ratio(sumw_var, tag='g_over_zvar', xedges=xedges, varied='denom', outtag=outtag)

        #dratio_dict['z_over_w'] = plot_ratio_anticorr(sumw_var, tag='z_over_w', xedges=xedges)

        plot_combined_scale_uncs(dratio_dict, xedges, outputrootfiles=outputrootfiles, outtag=outtag)

    else:
        sumw_var = {}
        for tag,regex in tag_regex.items():
            scale_var_list = scale_var_dict['gjets'] if tag == 'gjets' else scale_var_dict['wjet/dy']
            
            sumw_var[tag] = {}
            tup_combined = []

            for scale_var, scale_var_type in scale_var_list:
                tup, sumw_var[tag][scale_var_type] = get_scale_variations( acc=acc,
                                                                           regex=regex,
                                                                           tag=tag,
                                                                           scale_var=scale_var,
                                                                          )        

                tup_combined.append(tup)

                plot_ratio_vpt(tup, var=scale_var, tag=tag, outtag=outtag)
            
            tup_combined = np.array(tup_combined)
            plot_ratio_vpt_combined(tup_combined, tag=tag, outtag=outtag)
        
        # Dump contents of sumw_var into pickle file
        with open(pkl_file, 'wb+') as f:
            pickle.dump(tup, f)
            pickle.dump(sumw_var, f)


if __name__ == '__main__':
    main()
