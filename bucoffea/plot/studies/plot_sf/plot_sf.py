#!/usr/bin/env python
import os
import uproot
from matplotlib import pyplot as plt
import matplotlib
pjoin = os.path.join

font = {'family' : 'normal',
        'size'   : 14}

matplotlib.rc('font', **font)

name = {
    'gjets' : '$\gamma$ + jets',
    'wjets' : 'W + jets',
    'wjet' : 'W + jets',
    'dy' : 'DY',
    'zjets' : 'Z + jets',
}
def plot_nlo_ewk():
    outdir = './output'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    fig = plt.gcf()
    fig.clf()
    for tag in ['gjets','zjets','wjets']:
        f = uproot.open(f'../../../data/sf/theory/merged_kfactors_{tag}.root')
        h = f['kfactor_monojet_ewk']
        plt.plot(0.5*(h.bins[:,0]+h.bins[:,1]), h.values,'o-', label=name[tag])
    plt.ylabel('LO -> NLO EWK SF')
    plt.xlabel('Boson $p_{T}$ (GeV)')

    plt.legend()
    fig.savefig(pjoin(outdir, f'nlo_ewk.pdf'))

def plot_nlo_qcd():
    outdir = './output'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    f = uproot.open(f'../../../data/sf/theory/2017_gen_v_pt_qcd_sf.root')
    for selection in ['monojet','vbf']:
        fig = plt.gcf()
        fig.clf()
        for tag in ['wjet','dy']:
            h = f[f'{tag}_dress_{selection}']
            plt.plot(0.5*(h.bins[:,0]+h.bins[:,1]), h.values,'o-', label=name[tag])
        plt.ylabel('LO -> NLO QCD SF')
        plt.xlabel('Boson $p_{T}$ (GeV)')
        plt.grid(linestyle='--')
        plt.legend()
        fig.savefig(pjoin(outdir, f'nlo_qcd_{selection}.pdf'))

import numpy as np
def plot_nlo_qcd_forewk():
    outdir = './output'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    infiles = {
        'W' : 'kFactor_WToLNu_pT_Mjj.root',
        'DY' : 'kFactor_ZToNuNu_pT_Mjj.root'
    }

    for proc, infile in infiles.items():
        f = uproot.open(f'../../../data/sf/theory/{infile}')
        fig = plt.gcf()
        fig.clf()
        h = f['TH2F_kFactor']
        im = plt.gca().pcolormesh(h.edges[0], h.edges[1], h.values.T)
        plt.ylabel('LO -> NLO QCD SF')
        plt.xlabel('Boson $p_{T}$ (GeV)')
        cb = fig.colorbar(im)
        cb.set_label('LO $\\rightarrow$ NLO SF')
        im.set_clim(0.85,1.1)

        for ix in range(len(h.bins[0])):
            for iy in range(len(h.bins[1])):
                textcol = 'white' if h.values.T[iy, ix] < 0.95 else 'black'
                plt.gca().text(
                        np.mean(h.bins[0],axis=1)[ix],
                        np.mean(h.bins[1],axis=1)[iy],
                        f'  {h.values.T[iy, ix]:.3f}',
                        ha='center',
                        va='center',
                        color=textcol,
                        fontsize=9
                        )
        fig.tight_layout()
        fig.savefig(pjoin(outdir, f'nlo_qcd_for_ewk_{proc.lower()}.pdf'))

def plot_consistency():
    outdir = './output'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    fig = plt.gcf()
    fig.clf()
    for tag in ['zjets','wjets']:
        f = uproot.open(f'../../../data/sf/theory/merged_kfactors_{tag}.root')
        qcd = f['kfactor_monojet_qcd']
        ewk = f['kfactor_monojet_ewk']
        both = f['kfactor_monojet_qcd_ewk']

        plt.plot(0.5*(both.bins[:,0]+both.bins[:,1]), both.values,'o-', label=name[tag])
        plt.plot(0.5*(both.bins[:,0]+both.bins[:,1]), qcd.values*ewk.values,'o-', label=name[tag])
    plt.ylabel('LO -> NLO EWK SF')
    plt.xlabel('Boson $p_{T}$ (GeV)')

    plt.legend()
    fig.savefig(pjoin(outdir, f'consistency.pdf'))

# plot_consistency()
# plot_nlo_ewk()
# plot_nlo_qcd()
plot_nlo_qcd_forewk()