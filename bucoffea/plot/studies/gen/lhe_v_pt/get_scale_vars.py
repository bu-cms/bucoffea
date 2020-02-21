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

def get_scale_variations(acc, regex, tag, scale_var, scale_var_type, outputrootfile):
    '''Calculate the new k-factors with a scale weight variation.
       Dumps the ratio: 
       --- New k-factors with variation / Old (nominal) k-factors
       into the given output ROOT file.'''

    print(f'Working on: {tag}, {scale_var}')

    # Define rebinning
    if tag in ['wjet', 'dy']:
        vpt_ax_coarse = [0, 40, 80, 120, 160, 200, 240, 280, 320, 400, 520, 640, 760, 880,1200]
        vpt_ax_fine = list(range(0,400,40)) + list(range(400,1200,80))
        vpt_ax = hist.Bin('vpt','V $p_{T}$ (GeV)', vpt_ax_fine)
        mjj_ax = hist.Bin('mjj','M(jj) (GeV)', [0,200] + list(range(500,2500,500)))
    elif tag in ['gjets']:
        vpt_ax_coarse = [0, 40, 80, 120, 160, 200, 240, 280, 320, 400, 520, 640]
        vpt_ax_fine = list(range(0,400,40)) + list(range(400,1200,80)) 
        vpt_ax = hist.Bin('vpt','V $p_{T}$ (GeV)', vpt_ax_fine)
        mjj_ax = hist.Bin('mjj','M(jj) (GeV)',[0,200,500,1000,1500])

    # Set the correct pt type
    pt_tag = 'combined' if tag != 'gjets' else 'stat1'
    acc.load(f'gen_vpt_vbf_{pt_tag}')
    h = acc[f'gen_vpt_vbf_{pt_tag}']

    h = h.rebin('vpt', vpt_ax)
    h = h.rebin('mjj', mjj_ax)

    h = merge_extensions(h, acc, reweight_pu=False)
    scale_xs_lumi(h)
    h = merge_datasets(h)
    h = h[re.compile(regex)]

    lo = h[re.compile('.*HT.*')].integrate('dataset')
    nlo = h[re.compile('.*(LHE|amcat).*')].integrate('dataset')
    
    xaxis = lo.axes()[0]
    yaxis = lo.axes()[1]

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
    lo_1d = lo.integrate('mjj', mjj_slice, overflow='over')
    nlo_var_1d = nlo_var.integrate('mjj', mjj_slice, overflow='over')
    nlo_nom_1d = nlo_nom.integrate('mjj', mjj_slice, overflow='over')
    
    sumw_lo_1d = lo_1d.values(overflow='over')[()]
    sumw_nlo_var_1d = nlo_var_1d.values(overflow='over')[()]
    sumw_nlo_nom_1d = nlo_nom_1d.values(overflow='over')[()]

    # Calculate 1D scale factors, nominal and varied
    # as a function of V-pt
    sf_nom_1d = sumw_nlo_nom_1d / sumw_lo_1d
    sf_var_1d = sumw_nlo_var_1d / sumw_lo_1d

    # Calculate 1D variation ratio, as a function of V-pt
    var_ratio = sf_var_1d / sf_nom_1d
    
    # Calculate nominal 2D scale factor 
    sumw_lo = lo.values(overflow='over')[()]
    sumw_nlo_nom = nlo_nom.values(overflow='over')[()]

    sf_nom = sumw_nlo_nom / sumw_lo 

    tup = (var_ratio, yaxis.edges(overflow='over'))
 
    # Save to the ROOT file
    outputrootfile[f'{tag}_vbf_{scale_var_type}'] = tup

    # Return tuple containing the SF ratios and
    # NLO weights with and without variation
    return tup, (sumw_nlo_var_1d, sumw_nlo_nom_1d)

def plot_ratio_vpt(tup, var, tag):
    '''Given the tuple contatining the SF ratio (variational/nominal) and 
       bin edges, plot ratios in each mjj bin as a function of v-pt.'''
    ratio, x_edges = tup
    vpt_centers = ((x_edges + np.roll(x_edges,0))/2)[:-1]
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
    ax.plot(vpt_centers, ratio, 'o')
    ax.grid(True)

    # Save the figure
    outpath = './output/theory_variations'
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    outfile = pjoin(outpath, f'{tag}_kfac_ratio_{var}.png')
    fig.savefig(outfile)

def plot_ratio_vpt_combined(tup_combined, tag):
    '''Given the tup_combined contatining the SF ratio (variational/nominal) and 
       bin edges for all variations for a given process, 
       plot ratios in each mjj bin as a function of v-pt.'''
    ratio_combined, x_edges = tup_combined[:,0], tup_combined[0,1]
    vpt_centers = ((x_edges + np.roll(x_edges,0))/2)[:-1]
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
        ax.plot(vpt_centers, ratio_combined[idx], 'o', label=var_info[1])

    ax.set_xlim(-50.0, vpt_centers[-1]+50.0)
    ax.set_ylim(0.8, 1.2)
    ax.set_xlabel(r'$p_T(V) \ (GeV)$')
    ax.set_ylabel('Varied / Nominal SF')
    fig_title = variations['scale_1'][0]
    ax.set_title(fig_title)
    ax.grid(True)
    plt.legend()

    # Save the figure
    outpath = './output/theory_variations'
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    outfile = pjoin(outpath, f'{tag}_kfac_ratio_combined.png')
    fig.savefig(outfile)

def get_ratios(sumw_var, tag, var):
    '''Get ratio for two physics processes, for a given scale variation.'''
    # Figure out the processes
    if tag in ['zvar_over_w', 'z_over_wvar']:
        tag1 = 'dy'
        tag2 = 'wjet'
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

def plot_ratio(sumw_var, tag, xedges, outputrootfile, varied='num'):
    '''Plot ratio for two processes, for all variations.
       Specify which process is to be varied (num or denom).'''
    # List of all variations
    varlist = sumw_var['wjet'].keys()

    xcenters = ((xedges + np.roll(xedges,0))/2)[:-1]
    
    tag_to_ylabel = {
        'zvar_over_w' : r'$Z\rightarrow \ell \ell$ / $W\rightarrow \ell \nu$ (Var / Nom)',
        'z_over_wvar' : r'$Z\rightarrow \ell \ell$ / $W\rightarrow \ell \nu$ (Var / Nom)',
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
        else:
            raise ValueError('Invalid value for varied argument: Use num or denom.')

        # Calculate the ratio of ratios!
        dratio = ratio_var / ratio_nom

        ax.plot(xcenters, dratio, marker='o', label=var_to_label[var])
        
        dratios[var] = dratio

        # Save the double ratios in ROOT file
        tup = (dratio, xedges)
        outputrootfile[f'{tag}_{var}'] = tup 
    
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
    outdir = './output/theory_variations/scale/ratioplots'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    outpath = pjoin(outdir, f'{tag}_scaleunc_ratio.pdf')
    fig.savefig(outpath)

    print(f'Saved file in {outpath}')

    return dratios

def plot_combined_scale_uncs(dratio_dict, xedges):
    '''Plot combined scale variations for Z/W and photon/W ratios.
       Combined var: (var_num**2 + var_denom**2)**0.5'''
    to_be_combined = [('zvar_over_w', 'z_over_wvar'), ('gvar_over_z', 'g_over_zvar')]
    tags = ['z_over_w', 'g_over_z']
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
    tag_to_ylabel = {
        'z_over_w' : r'Combined Scale Unc: $Z\rightarrow \ell \ell$ / $W\rightarrow \ell \nu$',
        'g_over_z' : r'Combined Scale Unc: $\gamma$ + jets / $Z\rightarrow \ell \ell$',
    }
    xcenters = ((xedges + np.roll(xedges,0))/2)[:-1]
    for idx, pair in enumerate(to_be_combined):
        dratios_num = dratio_dict[pair[0]]
        dratios_denom = dratio_dict[pair[1]]
        fig, ax = plt.subplots(1,1)
        # Combine the variations
        for var_pair in var_pairs:
            combined_dratio_unc = 1 + np.sign(1-dratios_num[var_pair[0]]) * np.sqrt( (1-dratios_num[var_pair[0]] )**2 + ( 1-dratios_denom[var_pair[1]] )**2)
            ax.plot(xcenters, combined_dratio_unc, marker='o', label=var_to_label[var_pair[0]])
        ax.legend()
        ax.set_xlabel(r'$p_T(V)\ (GeV)$')
        ax.set_ylabel(tag_to_ylabel[tags[idx]])
        ax.grid(True)
        ax.set_xlim(-50,1200)
        ax.set_ylim(0.8,1.2)
        ax.plot([200,200], [0.8,1.2], 'k')
        ax.plot([-50,1200], [1,1], 'k--')
        outpath = './output/theory_variations/scale/ratioplots'
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
    
    outputrootfile = uproot.recreate('vbf_scale_var.root')

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

    pkl_dir = './output/theory_variations/scale'
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
        dratio_dict['zvar_over_w'] = plot_ratio(sumw_var, tag='zvar_over_w', xedges=xedges, varied='num', outputrootfile=outputrootfile)
        dratio_dict['z_over_wvar'] = plot_ratio(sumw_var, tag='z_over_wvar', xedges=xedges, varied='denom', outputrootfile=outputrootfile)
        dratio_dict['gvar_over_z'] = plot_ratio(sumw_var, tag='gvar_over_z', xedges=xedges, varied='num', outputrootfile=outputrootfile)
        dratio_dict['g_over_zvar'] = plot_ratio(sumw_var, tag='g_over_zvar', xedges=xedges, varied='denom', outputrootfile=outputrootfile)

        plot_combined_scale_uncs(dratio_dict, xedges)

    else:
        sumw_var = {}
        for tag,regex in tag_regex.items():
            scale_var_list = scale_var_dict['gjets'] if tag == 'gjets' else scale_var_dict['wjet/dy']
            
            sumw_var[tag] = {}
            tup_combined = []

            for scale_var, scale_var_type in scale_var_list:
                tup, sumw_var[tag][scale_var_type] = get_scale_variations( acc=acc,
                                                                           regex=regex,
                                                                           tag=tag    ,
                                                                           scale_var=scale_var,
                                                                           scale_var_type=scale_var_type,
                                                                           outputrootfile=outputrootfile
                                                                          )        

                tup_combined.append(tup)

                plot_ratio_vpt(tup, var=scale_var, tag=tag)
            
            tup_combined = np.array(tup_combined)
            plot_ratio_vpt_combined(tup_combined, tag=tag)
        
        # Dump contents of sumw_var into pickle file
        with open(pkl_file, 'wb+') as f:
            pickle.dump(tup, f)
            pickle.dump(sumw_var, f)


if __name__ == '__main__':
    main()
