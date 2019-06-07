import lz4.frame as lz4f
import cloudpickle
from coffea import hist

with lz4f.open("hists.cpkl.lz4", mode="r", compression_level=5) as fin:
    hists = cloudpickle.load(fin)



fig, ax, _ = hist.plot1d(hists["sr_met"],overlay="dataset")
# ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylim(0.1, 1e5)
fig.savefig("test.pdf")

# print(hists)