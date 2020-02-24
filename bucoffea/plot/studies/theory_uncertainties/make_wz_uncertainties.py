#!/usr/bin/env python

from coffea.util import load
from coffea.hist.plot import plot1d
from matplotlib import pyplot as plt
from klepto.archives import dir_archive
from bucoffea.plot.util import scale_xs_lumi, merge_extensions, merge_datasets
import uproot
import re
import sys
from coffea.hist.export import export1d
from coffea import hist
import ROOT as r
def from_coffea(inpath, outfile):

    acc = dir_archive(
                        inpath,
                        serialized=True,
                        compression=0,
                        memsize=1e3,
                        )

    # Merging, scaling, etc
    acc.load('sumw')
    acc.load('sumw_pileup')
    acc.load('nevents')
    mjj_ax = hist.Bin('mjj', r'$M_{jj}$ (GeV)', [200, 400, 600, 900, 1200, 1500, 2000, 2750, 3500, 5000])
    for distribution in ['mjj','mjj_unc']:
        acc.load(distribution)
        acc[distribution] = merge_extensions(
                                            acc[distribution],
                                            acc, 
                                            reweight_pu=not ('nopu' in distribution)
                                            )
        scale_xs_lumi(acc[distribution])
        acc[distribution] = merge_datasets(acc[distribution])
        acc[distribution] = acc[distribution].integrate(acc[distribution].axis('region'),'sr_vbf')
        acc[distribution] = acc[distribution].rebin(acc[distribution].axis('mjj'), mjj_ax)

    print(acc[distribution].axis('dataset').identifiers())
    histos = {}
    f = uproot.recreate(outfile)
    for year in [2017,2018]:
        
        # QCD V
        h_z = acc['mjj'][re.compile(f'ZJ.*HT.*{year}')].integrate('dataset')
        f[f'z_qcd_mjj_nominal_{year}'] = export1d(h_z)

        h_w = acc['mjj'][re.compile(f'W.*HT.*{year}')].integrate('dataset')
        f[f'w_qcd_mjj_nominal_{year}'] = export1d(h_w)

        # QCD Variations
        h_z_unc = acc['mjj_unc'][re.compile(f'ZJ.*HT.*{year}')].integrate('dataset')
        for unc in map(str, h_z_unc.axis('uncertainty').identifiers()):
            h = h_z_unc.integrate(h_z_unc.axis('uncertainty'), unc)
            f[f'z_qcd_mjj_{unc}_{year}'] = export1d(h)

        # EWK V
        h_z = acc['mjj'][re.compile(f'.*EWKZ.*{year}')].integrate('dataset')
        f[f'z_ewk_mjj_nominal_{year}'] = export1d(h_z)

        h_w = acc['mjj'][re.compile(f'.*EWKW.*{year}')].integrate('dataset')
        f[f'w_ewk_mjj_nominal_{year}'] = export1d(h_w)

        # EWK Variations
        h_z_unc = acc['mjj_unc'][re.compile(f'.*EWKZ.*{year}')].integrate('dataset')
        for unc in map(str, h_z_unc.axis('uncertainty').identifiers()):
            h = h_z_unc.integrate(h_z_unc.axis('uncertainty'), unc)
            f[f'z_ewk_mjj_{unc}_{year}'] = export1d(h)

def make_ratios(infile):
    f = r.TFile(infile)
    of = r.TFile(infile.replace('.root','_ratio.root'),'RECREATE')
    of.cd()

    for source in ['ewk','qcd']:
        for year in [2017,2018]:
            denominator = f.Get(f'w_{source}_mjj_nominal_{year}')
            for name in map(lambda x:x.GetName(), f.GetListOfKeys()):
                if not name.startswith(f'z_{source}'):
                    continue
                if not f"{year}" in name:
                    continue
                ratio = f.Get(name).Clone(f'ratio_{name}')
                ratio.Divide(denominator)
                ratio.SetDirectory(of)
                ratio.Write()
    of.Close()
    return str(of.GetName())
def make_uncertainties(infile):
    f = r.TFile(infile)
    of = r.TFile(infile.replace('_ratio','_ratio_unc'),'RECREATE')
    of.cd()
    for source in ['ewk','qcd']:
        for year in [2017,2018]:
            nominal = f.Get(f'ratio_z_{source}_mjj_nominal_{year}')
            for name in map(lambda x:x.GetName(), f.GetListOfKeys()):
                m = re.match(f'.*z_{source}_mjj_unc_(.*)_{year}', name)
                if not m:
                    continue
                variation_name = m.groups()[0]
                ratio = f.Get(name)
                variation = ratio.Clone(f'uncertainty_{name}')

                # New = variation / nominal -1
                variation.Divide(nominal)
                for i in range(1,variation.GetNbinsX()+1):
                    content = variation.GetBinContent(i)
                    variation.SetBinContent(i, content-1)

                variation.SetDirectory(of)
                variation.Write()
                
                ratio.SetDirectory(of)
                ratio.Write()
    of.Close()
    return str(of.GetName())
import os
pjoin = os.path.join
def main():
    outdir = './output'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    inpath = sys.argv[1] 
    outfile = pjoin(outdir, f'vbf_z_w_theory_unc.root')
    from_coffea(inpath, outfile)
    outfile = make_ratios(outfile)
    make_uncertainties(outfile)
if __name__ == "__main__":
    main()
# print(h)
# fig, ax, _ = plot1d(
#                     h,
#                     overlay='uncertainty'
# )
# # plt.yscale('log')
# plt.ylim(2e3,3e3)
# plt.xlim(0,5e2)
# fig.savefig('test.pdf')
