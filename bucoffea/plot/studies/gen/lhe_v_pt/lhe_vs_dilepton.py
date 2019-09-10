#!/usr/bin/env python
import re
from coffea.util import load
from coffea import hist
from matplotlib import pyplot as plt
from bucoffea.plot.util import acc_from_dir
import copy
from bucoffea.plot.util import merge_extensions, merge_datasets, scale_xs_lumi
acc = acc_from_dir('input/das_lhevpt_v2')

data_err_opts = {
        'linestyle':'none',
        'marker': '.',
        'markersize': 10.,
        'color':'k',
        'elinewidth': 1,
        'emarker': '_'
    }

h = copy.deepcopy(acc['gen_vpt'])

h = merge_extensions(h, acc, reweight_pu=False)
scale_xs_lumi(h)
h = merge_datasets(h)

h = h.integrate('weight_type','nominal')
h = h.integrate('weight_index',slice(-0.5,0.5))
h = h[re.compile('.*DY.*HT.*')].integrate('dataset')

new_ax = hist.Bin('vpt','Gen V $p_{T}$ (GeV)',list(range(80,800,40))+list(range(800,2000,100)))
h = h.rebin(h.axis('vpt'), new_ax)
print(h)

fig, (ax, rax) = plt.subplots(2, 1, figsize=(7,7), gridspec_kw={"height_ratios": (3, 1)}, sharex=True)

hist.plot1d(
    h,
    overlay='type',
    overflow='all',
    ax = ax,
    clear=False,
    binwnorm=True)
ax.legend()

hist.plotratio(h['dilepton'].integrate('type'), h['nano'].integrate('type'),
    ax=rax,
    denom_fill_opts={},
    guide_opts={},
    unc='num',
    overflow='all',
    error_opts=data_err_opts,
    label='status 1 dilepton / LHE',
    clear=False
    )

# rax.set_ylim(0.9,1.1)
fig.savefig('test_lhe.pdf')