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

def plot_nnlo_qcd():
    outdir = './output'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    f = uproot.open(f'../../../data/sf/theory/lindert_qcd_nnlo_sf.root')
    processes = {
        'aj' : r'$\gamma$',
        'eej' : 'DY',
        'evj' : 'W',
        'vvj' : r'Z($\nu\nu$)',
    }
    fig = plt.gcf()
    fig.clf()
    for proc, name in processes.items():
        h = f[proc]
        plt.plot(0.5*(h.bins[:,0]+h.bins[:,1]), h.values,'o-', label=name,fillstyle='none')
    plt.ylabel('NLO -> NNLO QCD SF')
    plt.xlabel('Boson $p_{T}$ (GeV)')
    plt.grid(linestyle='--')
    plt.legend()
    plt.gca().set_xscale('log')
    fig.savefig(pjoin(outdir, f'nnlo_qcd.pdf'))

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
plot_nnlo_qcd()