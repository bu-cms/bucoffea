from coffea.util import load
from coffea import hist
from matplotlib import pyplot as plt
# Load input
output = load('monojet.coffea')


h=output['genvpt_check'].project("dataset")
oldax = h.axis("vpt")
newax = hist.Bin("vpt",r"$p_{T}^{V}$ (GeV)", 25, 0, 2000)

h = h.rebin(oldax, newax)
fig, (ax, rax) = plt.subplots(2, 1, figsize=(7,7), gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
fig, ax, _ = hist.plot1d(
    h,
    overlay='type',
    ax=ax)


# Ratio plot
data_err_opts = {
    'linestyle':'none',
    'marker': '.',
    'markersize': 10.,
    'color':'k',
    'elinewidth': 1,
    'emarker': '_'
}

a=h.axis('type')
h1=h.project(a,"BU")
h2=h.project(a,"Nano")
hist.plotratio(h2, h1, 
               ax=rax,
               denom_fill_opts={},
               guide_opts={},
               unc='num',
               error_opts=data_err_opts
              )

# Aesthetics
ax.set_yscale("log")
ax.set_ylim(1,1e6)
rax.set_ylim(0.9,1.1)
fig.savefig("genvpt.pdf")