#!/usr/bin/env python

import uproot
from matplotlib import pyplot as plt


f = uproot.open('2017_gen_v_pt_stat1_qcd_sf.root')


for process in ['wjet','dy']:
    x = {}
    y = {}

    for selection in ['inclusive','monojet','vbf']:
        h = f[f'{process}_{selection}']
        
        x[selection] = 0.5*(h.bins[:,0]+h.bins[:,1])
        y[selection] = h.values
    
    fig, (ax, rax) = plt.subplots(2, 1, figsize=(7,7), gridspec_kw={"height_ratios": (3, 1)}, sharex=True)

    for sel in x.keys():
        ax.plot(x[sel],y[sel],'-o',label=sel)
        rax.plot(x[sel],y[sel] / y['inclusive'],'-o')

    rax.set_ylim(0.9,1.1)
    ax.legend()
    fig.savefig(f'output/selection_comparison/selection_comparison_{process}.pdf')
    plt.close(fig)