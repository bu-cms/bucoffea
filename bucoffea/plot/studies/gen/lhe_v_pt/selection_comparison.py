#!/usr/bin/env python

import uproot
from matplotlib import pyplot as plt
import matplotlib.ticker

f = uproot.open('2017_gen_v_pt_qcd_sf.root')


for process in ['wjet','dy','gjets']:
    x = {}
    y = {}

    for selection in ['inclusive','monojet','vbf']:
        if process in ['wjet','dy']:
            h = f[f'{process}_dress_{selection}']
        elif process in ['gjets']:
            h = f[f'{process}_stat1_{selection}']
        
        x[selection] = 0.5*(h.bins[:,0]+h.bins[:,1])
        y[selection] = h.values
    
    fig, (ax, rax) = plt.subplots(2, 1, figsize=(7,7), gridspec_kw={"height_ratios": (3, 1)}, sharex=True)

    for sel in x.keys():
        ax.plot(x[sel],y[sel],'-o',label=sel)
        rax.plot(x[sel],y[sel] / y['inclusive'],'-o')

    rax.set_ylim(0.95,1.15)
    ax.legend(title='Selection', fontsize=14)
    ax.set_ylabel('LO -> NLO SF', fontsize=14)
    ax.set_xlabel('$p_{T}$ (V) (GeV)', fontsize=14)
    rax.set_xlabel('$p_{T}$ (V) (GeV)', fontsize=14)
    rax.set_ylabel('Ratio to inclusive', fontsize=14)


    loc1 = matplotlib.ticker.MultipleLocator(base=0.02)
    rax.yaxis.set_major_locator(loc1)
    rax.grid(axis='y',which='major',linestyle='--')
    fig.savefig(f'output/selection_comparison/selection_comparison_{process}.pdf')
    plt.close(fig)