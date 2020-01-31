
import os
from matplotlib import pyplot as plt
import numpy as np
from scipy.optimize import curve_fit, minimize
pjoin = os.path.join
from bucoffea.plot.util import fig_ratio, lumi
from bucoffea.helpers import sigmoid

def load(tag, dataset, year):
    fname = f'output/gamma/table_g_{tag}_photon_pt0_{dataset}_{year}.txt'
    data = np.loadtxt(fname, skiprows=2)
    x = data[:,2]
    eff = data[:,5]
    eff_up = data[:,6]
    eff_down = data[:,7]

    return x, eff, eff_up, eff_down


colors = {
    'mc' : 'mediumblue',
    'data' : 'darkorange',
}
def fit(tag, year):
    outdir = './output/gamma/fit/'
    try:
        os.makedirs(outdir)
    except FileExistsError:
        pass

    x, eff, eff_up, eff_down = {}, {}, {}, {}
    x['data'], eff['data'], eff_up['data'], eff_down['data'] = load(tag, 'JetHT', year)
    x['mc'], eff['mc'], eff_up['mc'], eff_down['mc'] = load(tag, 'GJets_HT_MLM', year)

    pars = {}
    cross = {}
    for key in ['data','mc']:
        pars[key], _ = curve_fit(sigmoid, x[key], eff[key], sigma=0.5*(eff_up[key]-eff_down[key]), p0=[1e-3,200,0.1,1])    
        cross[key] = minimize( lambda x: np.abs(sigmoid(x, *pars[key]) - 0.95),x0=230)

    fig, ax, rax = fig_ratio()

    xinterp = np.linspace(min(x['data']), max(x['data']), 1000)

    handles = []
    ax.errorbar(x['data'], eff['data'], 0.5*(eff_up['data']-eff_down['data']),fmt='o',label='Data',color=colors["data"])
    ax.errorbar(x['mc'], eff['mc'], 0.5*(eff_up['mc']-eff_down['mc']),fmt='s',label='MC',fillstyle='none',color=colors["mc"])
    ax.plot(xinterp, sigmoid(xinterp, *pars['data']), label='Data fit',color=colors["data"],zorder=-1)
    ax.plot(xinterp, sigmoid(xinterp, *pars['mc']), label='MC fit',color=colors["mc"],zorder=-1,linestyle='--')

    ax.set_ylim(0.,1.1)
    ax.set_xlim(100,1100)
    ax.legend()


    ax.text(350,.4,'f(x) = c + (d-c) / (1 + exp(-a * (x-b)))')
    ax.text(
            300,
            0.1,
            '\n'.join([
                f"a = {pars['data'][0]:.3f} / GeV",
                f"b = {pars['data'][1]:.2f} GeV",
                f"c = {pars['data'][2]:.3f}",
                f"d = {pars['data'][3]:.3f}"
            ]),
            color=colors['data']
            )
    ax.text(
            600,
            0.1,
            '\n'.join([
                f"a = {pars['mc'][0]:.3f} / GeV",
                f"b = {pars['mc'][1]:.2f} GeV",
                f"c = {pars['mc'][2]:.3f}",
                f"d = {pars['mc'][3]:.3f}"
            ]),
            color=colors['mc']
            )

    ax.text(700,0.8, 
                    "\n".join([
                              f"Data > 95% @ {cross['data'].x[0]:.0f} GeV",
                              f"MC > 95% @ {cross['mc'].x[0]:.0f} GeV",
                    ])
                              )
    ax.text(1., 1., r"%.1f fb$^{-1}$ (13 TeV)" % lumi(year),
                fontsize=16,
                horizontalalignment='right',
                verticalalignment='bottom',
                transform=ax.transAxes
                )
    rax.set_ylim(0.95,1.1)
    rax.grid()
    # rax.plot([cross['data'].x[0],cross['data'].x[0]], [0.8,1.05],color='k',)
    # ax.plot([cross['data'].x[0],cross['data'].x[0]], [0.,1.05],color='k',linestyle='--')
    # ax.plot([cross['mc'].x[0],cross['mc'].x[0]], [0.9,1.05],color='r')


    # ax.plot([215,215],[0.9,1.05])

    rax.errorbar(x["data"], eff["data"] / eff["mc"], 0.5*(eff_up["data"] - eff_down["data"]) / eff["mc"],fmt='o',label='Data / MC',color=colors['data'])
    rxinterp = np.linspace(cross['data'].x[0], max(xinterp),1000)

    
    rax.plot(rxinterp, sigmoid(rxinterp, *pars['data']) / sigmoid(rxinterp, *pars['mc']),label=f"Data / MC fit ratio, plateau at {100*(sigmoid(rxinterp, *pars['data']) / sigmoid(rxinterp, *pars['mc']))[-1]:.1f} %",color='k')
    rax.plot(rxinterp, 0.99*sigmoid(rxinterp, *pars['data']) / sigmoid(rxinterp, *pars['mc']), label='1% uncertainty on fit',linestyle='--',color='gray')
    rax.plot(rxinterp, 1.01*sigmoid(rxinterp, *pars['data']) / sigmoid(rxinterp, *pars['mc']),linestyle='--',color='gray')
    # rxinterp2 = np.linspace(min(xinterp),cross['data'].x[0], 1000)
    # rax.plot(rxinterp2, sigmoid(rxinterp2, *pars['data']) / sigmoid(rxinterp2, *pars['mc']),color='k',linestyle=':')

    rax.legend()
    rax.set_ylabel("Ratio")
    ax.set_ylabel("Trigger efficiency")
    ax.set_xlabel("Photon $p_{T}$ (GeV)")
    rax.set_xlabel("Photon $p_{T}$ (GeV)")
    ax.figure.savefig(pjoin(outdir, f'fit_{tag}_{year}.pdf'))
    ax.figure.clf()



tags = ['HLT_PFHT1050_ht1500']

for tag in tags:
        for year in [2017,2018]:
            fit(tag, year)

