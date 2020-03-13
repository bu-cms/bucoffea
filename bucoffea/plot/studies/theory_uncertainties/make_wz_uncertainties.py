#!/usr/bin/env python

from coffea.util import load
from coffea.hist.plot import plot1d
from matplotlib import pyplot as plt
from klepto.archives import dir_archive
from bucoffea.plot.util import scale_xs_lumi, merge_extensions, merge_datasets
import uproot
import re
import sys
import numpy as np
from coffea.hist.export import export1d
from coffea import hist
import ROOT as r
from pprint import pprint
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
    for distribution in ['mjj','mjj_unc', 'mjj_noewk']:
        acc.load(distribution)
        acc[distribution] = merge_extensions(
                                            acc[distribution],
                                            acc, 
                                            reweight_pu=not ('nopu' in distribution)
                                            )
        scale_xs_lumi(acc[distribution])
        acc[distribution] = merge_datasets(acc[distribution])
        acc[distribution] = acc[distribution].rebin(acc[distribution].axis('mjj'), mjj_ax)

    pprint(acc[distribution].axis('dataset').identifiers())
    f = uproot.recreate(outfile)
    for year in [2017,2018]:
        # QCD V
        h_z = acc['mjj'][re.compile(f'ZJetsToNuNu.*HT.*{year}')].integrate('region', 'sr_vbf').integrate('dataset')
        f[f'z_qcd_mjj_nominal_{year}'] = export1d(h_z)

        h_w = acc['mjj'][re.compile(f'WJetsToLNu.*HT.*{year}')].integrate('region', 'sr_vbf').integrate('dataset')
        f[f'w_qcd_mjj_nominal_{year}'] = export1d(h_w)

        h_ph = acc['mjj'][re.compile(f'GJets_DR-0p4.*HT.*{year}')].integrate('region', 'cr_g_vbf').integrate('dataset')
        f[f'gjets_qcd_mjj_nominal_{year}'] = export1d(h_ph)

        # Scale + PDF variations for QCD Z 
        h_z_unc = acc['mjj_unc'][re.compile(f'ZJ.*HT.*{year}')].integrate('region', 'sr_vbf').integrate('dataset')
        for unc in map(str, h_z_unc.axis('uncertainty').identifiers()):
            if 'goverz' in unc or 'ewkcorr' in unc:
                continue
            h = h_z_unc.integrate(h_z_unc.axis('uncertainty'), unc)
            f[f'z_qcd_mjj_{unc}_{year}'] = export1d(h)

        # EWK variations for QCD Z
        # Get EWK down variation first
        h_z_unc_ewk = acc['mjj_noewk'][re.compile(f'ZJetsToNuNu.*HT.*{year}')].integrate('region', 'sr_vbf').integrate('dataset')
        f[f'z_qcd_mjj_unc_w_ewkcorr_overz_common_down_{year}'] = export1d(h_z_unc_ewk)

        # Get EWK up variation
        h_z_unc_ewk.scale(-1)
        h_z_diff = h_z.copy().add(h_z_unc_ewk)
        h_z_unc_ewk_down = h_z.add(h_z_diff) 
        f[f'z_qcd_mjj_unc_w_ewkcorr_overz_common_up_{year}'] = export1d(h_z_unc_ewk_down)

        # EWK variations for QCD W
        # Get EWK down variation first
        h_w_unc_ewk = acc['mjj_noewk'][re.compile(f'WJetsToLNu.*HT.*{year}')].integrate('region', 'sr_vbf').integrate('dataset')
        f[f'w_qcd_mjj_unc_w_ewkcorr_overz_common_down_{year}'] = export1d(h_w_unc_ewk)

        # Get EWK up variation
        h_w_unc_ewk.scale(-1)
        h_w_diff = h_w.copy().add(h_w_unc_ewk)
        h_w_unc_ewk_down = h_w.add(h_w_diff)
        f[f'w_qcd_mjj_unc_w_ewkcorr_overz_common_up_{year}'] = export1d(h_w_unc_ewk_down)

        # Scale + PDF variations for QCD photons
        h_ph_unc = acc['mjj_unc'][re.compile(f'GJets_DR-0p4.*HT.*{year}')].integrate('region', 'cr_g_vbf').integrate('dataset')
        for unc in map(str, h_ph_unc.axis('uncertainty').identifiers()):
            if 'zoverw' in unc or 'ewkcorr' in unc:
                continue
            h = h_ph_unc.integrate(h_ph_unc.axis('uncertainty'), unc)
            f[f'gjets_qcd_mjj_{unc}_{year}'] = export1d(h)

        # EWK variations for QCD photons
        # Get EWK down variation first
        h_ph_unc_ewk = acc['mjj_noewk'][re.compile(f'GJets_DR-0p4.*HT.*{year}')].integrate('region', 'cr_g_vbf').integrate('dataset')
        f[f'gjets_qcd_mjj_unc_w_ewkcorr_overz_common_down_{year}'] = export1d(h_ph_unc_ewk)

        # Get EWK up variation
        h_ph_unc_ewk.scale(-1)
        h_ph_diff = h_ph.copy().add(h_ph_unc_ewk)
        h_ph_unc_ewk_down = h_ph.add(h_ph_diff)
        f[f'gjets_qcd_mjj_unc_w_ewkcorr_overz_common_up_{year}'] = export1d(h_ph_unc_ewk_down)

        # EWK V
        h_z = acc['mjj'][re.compile(f'.*EWKZ.*{year}')].integrate('region', 'sr_vbf').integrate('dataset')
        f[f'z_ewk_mjj_nominal_{year}'] = export1d(h_z)

        h_w = acc['mjj'][re.compile(f'.*EWKW.*{year}')].integrate('region', 'sr_vbf').integrate('dataset')
        f[f'w_ewk_mjj_nominal_{year}'] = export1d(h_w)

        h_ph = acc['mjj'][re.compile(f'GJets_SM_5f_EWK.*{year}')].integrate('region', 'cr_g_vbf').integrate('dataset')
        f[f'gjets_ewk_mjj_nominal_{year}'] = export1d(h_ph)
        print(h_ph.values())

        # Scale + PDF variations for EWK Z
        h_z_unc = acc['mjj_unc'][re.compile(f'.*EWKZ.*{year}')].integrate('region', 'sr_vbf').integrate('dataset')
        for unc in map(str, h_z_unc.axis('uncertainty').identifiers()):
            if 'goverz' in unc or 'ewkcorr' in unc:
                continue
            h = h_z_unc.integrate(h_z_unc.axis('uncertainty'), unc)
            f[f'z_ewk_mjj_{unc}_{year}'] = export1d(h)

        # Scale + PDF variations for EWK photons
        h_ph_unc = acc['mjj_unc'][re.compile(f'GJets_SM.*{year}')].integrate('region', 'cr_g_vbf').integrate('dataset')
        for unc in map(str, h_ph_unc.axis('uncertainty').identifiers()):
            if 'zoverw' in unc or 'ewkcorr' in unc:
                continue
            h = h_ph_unc.integrate(h_ph_unc.axis('uncertainty'), unc)
            f[f'gjets_ewk_mjj_{unc}_{year}'] = export1d(h)

def make_ratios(infile):
    f = r.TFile(infile)
    of = r.TFile(infile.replace('.root','_ratio.root'),'RECREATE')
    of.cd()

    # Z / W ratios (scale + PDF variations)
    for source in ['ewk','qcd']:
        for year in [2017,2018]:
            denominator = f.Get(f'w_{source}_mjj_nominal_{year}')
            for name in map(lambda x:x.GetName(), f.GetListOfKeys()):
                if not name.startswith(f'z_{source}'):
                    continue
                if not f"{year}" in name or 'ewkcorr' in name:
                    continue
                ratio = f.Get(name).Clone(f'ratio_{name}')
                ratio.Divide(denominator)
                ratio.SetDirectory(of)
                ratio.Write()
    
    # Z / W ratios (up and down EWK variations)
    for year in [2017,2018]:
        for vartype in ['up', 'down']:
            varied_z_name = f'z_qcd_mjj_unc_w_ewkcorr_overz_common_{vartype}_{year}'
            varied_w = f.Get(f'w_qcd_mjj_unc_w_ewkcorr_overz_common_{vartype}_{year}')
            varied_ratio = f.Get(varied_z_name).Clone(f'ratio_{varied_z_name}')
            varied_ratio.Divide(varied_w)
            varied_ratio.SetDirectory(of)
            varied_ratio.Write()  

        nominal_z_name = f'z_qcd_mjj_nominal_{year}'
        nominal_w = f.Get(f'w_qcd_mjj_nominal_{year}')
        nominal_ratio = f.Get(nominal_z_name).Clone(f'ratio_{nominal_z_name}')
        nominal_ratio.Divide(nominal_w)
        nominal_ratio.SetDirectory(of)    
        nominal_ratio.Write()  

    # GJets / Z ratios (scale + PDF variations)
    for source in ['ewk','qcd']:
        for year in [2017,2018]:
            denominator = f.Get(f'z_{source}_mjj_nominal_{year}')
            for name in map(lambda x:x.GetName(), f.GetListOfKeys()):
                if not name.startswith(f'gjets_{source}'):
                    continue
                if not f"{year}" in name or 'ewkcorr' in name:
                    continue
                ratio = f.Get(name).Clone(f'ratio_{name}')
                ratio.Divide(denominator)
                ratio.SetDirectory(of)
                ratio.Write()

    # GJets / Z ratios (up and down EWK variations)
    for year in [2017,2018]:
        for vartype in ['up', 'down']:
            varied_g_name = f'gjets_qcd_mjj_unc_w_ewkcorr_overz_common_{vartype}_{year}'
            varied_z = f.Get(f'z_qcd_mjj_unc_w_ewkcorr_overz_common_{vartype}_{year}')
            varied_ratio = f.Get(varied_g_name).Clone(f'ratio_{varied_g_name}')
            varied_ratio.Divide(varied_z)
            varied_ratio.SetDirectory(of)
            varied_ratio.Write()

        nominal_g_name = f'gjets_qcd_mjj_nominal_{year}'
        nominal_z = f.Get(f'z_qcd_mjj_nominal_{year}')
        nominal_ratio = f.Get(nominal_g_name).Clone(f'ratio_{nominal_g_name}')
        nominal_ratio.Divide(nominal_z)
        nominal_ratio.SetDirectory(of)
        nominal_ratio.Write()

    of.Close()
    return str(of.GetName())

def make_uncertainties(infile):
    f = r.TFile(infile)
    of = r.TFile(infile.replace('_ratio','_ratio_unc'),'RECREATE')
    of.cd()

    # Uncertainty in Z / W ratios (scale + PDF variations)
    for source in ['ewk','qcd']:
        for year in [2017,2018]:
            nominal = f.Get(f'ratio_z_{source}_mjj_nominal_{year}')
            for name in map(lambda x:x.GetName(), f.GetListOfKeys()):
                m = bool(re.match(f'.*z_{source}_mjj_unc_(.*)_{year}', name)) and 'ewkcorr' not in name
                if not m:
                    continue
                ratio = f.Get(name)
                variation = ratio.Clone(f'uncertainty_{name}')

                # Content: Varied ratio / Nominal ratio
                variation.Divide(nominal)

                variation.SetDirectory(of)
                variation.Write()
                
                ratio.SetDirectory(of)
                ratio.Write()
    
    # Uncertainty in Z / W ratios (up and down EWK variations)
    for year in [2017,2018]:
        nominal = f.Get(f'ratio_z_qcd_mjj_nominal_{year}')
        for vartype in ['up', 'down']:
            varied_name = f'ratio_z_qcd_mjj_unc_w_ewkcorr_overz_common_{vartype}_{year}'
            varied = f.Get(varied_name)
            # Variation: (varied Z / W) / (nominal Z / W)
            variation = varied.Clone(f'uncertainty_{varied_name}')
            variation.Divide(nominal)
    
            variation.SetDirectory(of)
            variation.Write()

            varied.SetDirectory(of)
            varied.Write()
    
    # Copy the EWK uncertainty for QCD Z / W ratio
    # apply to the EWK Z / W ratio (for now)
    for year in [2017,2018]:
        for vartype in ['up', 'down']:
            qcd_unc_name = f'uncertainty_ratio_z_qcd_mjj_unc_w_ewkcorr_overz_common_{vartype}_{year}'
            qcd_unc = of.Get(qcd_unc_name)
            ewk_unc = qcd_unc.Clone(f'{qcd_unc_name.replace("qcd", "ewk")}')

            ewk_unc.SetDirectory(of)
            ewk_unc.Write()

    # Uncertainty in GJets / Z ratios (scale + PDF variations)
    for source in ['ewk','qcd']:
        for year in [2017,2018]:
            nominal = f.Get(f'ratio_gjets_{source}_mjj_nominal_{year}')
            for name in map(lambda x:x.GetName(), f.GetListOfKeys()):
                m = bool(re.match(f'.*gjets_{source}_mjj_unc_(.*)_{year}', name)) and 'ewkcorr' not in name
                if not m:
                    continue
                ratio = f.Get(name)
                variation = ratio.Clone(f'uncertainty_{name}')

                # Content: Varied ratio / Nominal ratio
                variation.Divide(nominal)

                variation.SetDirectory(of)
                variation.Write()
                
                ratio.SetDirectory(of)
                ratio.Write()

    # Uncertainty in GJets / Z ratios (up and down EWK variations)
    for year in [2017,2018]:
        nominal = f.Get(f'ratio_gjets_qcd_mjj_nominal_{year}')
        for vartype in ['up', 'down']:
            varied_name = f'ratio_gjets_qcd_mjj_unc_w_ewkcorr_overz_common_{vartype}_{year}'
            varied = f.Get(varied_name)
            # Variation: (varied Z / W) / (nominal Z / W)
            variation = varied.Clone(f'uncertainty_{varied_name}')
            variation.Divide(nominal)
    
            variation.SetDirectory(of)
            variation.Write()

            varied.SetDirectory(of)
            varied.Write()
    
    # Copy the EWK uncertainty for QCD GJets / Z ratio
    # apply to the EWK GJets / Z ratio (for now)
    for year in [2017,2018]:
        for vartype in ['up', 'down']:
            qcd_unc_name = f'uncertainty_ratio_gjets_qcd_mjj_unc_w_ewkcorr_overz_common_{vartype}_{year}'
            qcd_unc = of.Get(qcd_unc_name)
            ewk_unc = qcd_unc.Clone(f'{qcd_unc_name.replace("qcd", "ewk")}')

            ewk_unc.SetDirectory(of)
            ewk_unc.Write()

    of.Close()
    return str(of.GetName())

import os
pjoin = os.path.join

def main():
    inpath = sys.argv[1] 

    # Get the output tag for output directory naming
    if inpath.endswith('/'):
        outtag = inpath.split('/')[-2]
    else:
        outtag = inpath.split('/')[-1]
    
    outdir = f'./output/{outtag}'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outfile = pjoin(outdir, f'vbf_z_w_gjets_theory_unc.root')
    from_coffea(inpath, outfile)
    outfile = make_ratios(outfile)
    make_uncertainties(outfile)
if __name__ == "__main__":
    main()
