#!/usr/bin/env python

import os
import re
import sys
from pprint import pprint
from bucoffea.plot.util import merge_datasets, merge_extensions, scale_xs_lumi 
from bucoffea.plot.stack_plot import Style, make_plot
from bucoffea.plot.cr_ratio_plot import cr_ratio_plot
from bucoffea.plot.style import plot_settings

from collections import defaultdict
from klepto.archives import dir_archive

def plot(inpath):
        indir=os.path.abspath(inpath)

        # The processor output is stored in an
        # 'accumulator', which in our case is
        # just a dictionary holding all the histograms
        # Put all your *coffea files into 'indir' and
        # pass the directory as an argument here.
        # All input files in the directory will
        # automatically be found, merged and read.
        # The merging only happens the first time
        # you run over a specific set of inputs.
        acc = dir_archive(
                          inpath,
                          serialized=True,
                          compression=0,
                          memsize=1e3
                          )
        # Get a settings dictionary that details
        # which plots to make for each region,
        # what the axis limits are, etc
        # Can add plots by extending the dictionary
        # Or modify axes ranges, etc
        settings = plot_settings()

        merged = set()

        # Separate plots per year
        for year in [2017,2018]:
            # The data to be used for each region
            # Muon regions use MET,
            # electron+photon regions use EGamma
            # ( EGamma = SingleElectron+SinglePhoton for 2017)
            data = {
                'sr_vbf' : None,
                'cr_1m_vbf' : f'MET_{year}',
                'cr_2m_vbf' : f'MET_{year}',
                'cr_1e_vbf' : f'EGamma_{year}',
                'cr_2e_vbf' : f'EGamma_{year}',
                'cr_g_vbf' : f'EGamma_{year}',
            }

            # Same for MC selection
            # Match datasets by regular expressions
            # Here for LO V samples (HT binned)
            mc_lo = {
                'sr_vbf' : re.compile(f'(ZJetsToNuNu.*|EW.*|TTJets.*FXFX.*|Diboson.*|ST.*|QCD_HT.*|.*DYJetsToLL_M-50_HT_MLM.*|.*WJetsToLNu.*HT.*).*{year}'),
                'cr_1m_vbf' : re.compile(f'(EW.*|TTJets.*FXFX.*|Diboson.*|ST.*|QCD_HT.*|.*DYJetsToLL_M-50_HT_MLM.*|.*WJetsToLNu.*HT.*).*{year}'),
                'cr_1e_vbf' : re.compile(f'(EW.*|TTJets.*FXFX.*|Diboson.*|ST.*|QCD_HT.*|.*DYJetsToLL_M-50_HT_MLM.*|.*WJetsToLNu.*HT.*).*{year}'),
                'cr_2m_vbf' : re.compile(f'(EW.*|TTJets.*FXFX.*|Diboson.*|ST.*|QCD_HT.*|.*DYJetsToLL_M-50_HT_MLM.*).*{year}'),
                'cr_2e_vbf' : re.compile(f'(EW.*|TTJets.*FXFX.*|Diboson.*|ST.*|QCD_HT.*|.*DYJetsToLL_M-50_HT_MLM.*).*{year}'),
                'cr_g_vbf' : re.compile(f'(GJets_(DR-0p4|SM).*|QCD_HT.*|WJetsToLNu.*HT.*).*{year}'),
            }

            # Want to compare LO and NLO,
            # so do same thing for NLO V samples
            # All non-V samples remain the same
            mc_nlo = {
                    'sr_vbf' : re.compile(f'(ZJetsToNuNu.*|EW.*|TTJets.*FXFX.*|Diboson.*|ST.*|QCD_HT.*|.*DYJetsToLL_M-50_HT_MLM.*|.*WJetsToLNu.*FXFX.*).*{year}'),
                    'cr_1m_vbf' : re.compile(f'(EW.*|TTJets.*FXFX.*|Diboson.*|ST.*|QCD_HT.*|.*DY.*FXFX.*|.*WJetsToLNu.*FXFX.*).*{year}'),
                    'cr_1e_vbf' : re.compile(f'(EW.*|TTJets.*FXFX.*|Diboson.*|ST.*|QCD_HT.*|.*DY.*FXFX.*|.*WJetsToLNu.*FXFX.*).*{year}'),
                    'cr_2m_vbf' : re.compile(f'(EW.*|TTJets.*FXFX.*|Diboson.*|ST.*|QCD_HT.*|.*DY.*FXFX.*).*{year}'),
                    'cr_2e_vbf' : re.compile(f'(EW.*|TTJets.*FXFX.*|Diboson.*|ST.*|QCD_HT.*|.*DY.*FXFX.*).*{year}'),
                    'cr_g_vbf' : re.compile(f'(GJets_(DR-0p4|SM).*|QCD_HT.*|W.*FXFX.*).*{year}'),
            }

            regions = list(mc_lo.keys())
            # Remove signal region, no need in ratio plots
            regions.remove('sr_vbf')

            # Make control region ratio plots for both
            # LO and NLO. Can be skipped if you only
            # want data / MC agreement plots.
            outdir = f'./output/{os.path.basename(indir)}/ratios'

            # Load ingredients from cache
            acc.load('mjj')
            acc.load('sumw')
            acc.load('sumw_pileup')
            acc.load('nevents')
            cr_ratio_plot(acc, year=year,tag='losf',outdir=outdir, mc=mc_lo, regions=regions, distribution='mjj')
            cr_ratio_plot(acc, year=year,tag='nlo',outdir=outdir, mc=mc_nlo, regions=regions, distribution='mjj')

            # Data / MC plots are made here
            # Loop over all regions
            for region in mc_lo.keys():
                ratio = True if region != 'sr_vbf' else False 
                # Make separate output direcotry for each region
                outdir = f'./output/{os.path.basename(indir)}/{region}'
                # Settings for this region
                plotset = settings[region]

                # Loop over the distributions
                for distribution in plotset.keys():
                    # Load from cache
                    if not distribution in merged:
                        acc.load(distribution)

                        if not distribution in acc.keys():
                            print(f"WARNING: Distribution {distribution} not found in input files.")
                            continue
                        acc[distribution] = merge_extensions(acc[distribution], acc, reweight_pu=not ('nopu' in distribution))
                        scale_xs_lumi(acc[distribution]) 
                        acc[distribution] = merge_datasets(acc[distribution])
                        acc[distribution].axis('dataset').sorting = 'integral'
                        merged.add(distribution)
                    try:
                        # The heavy lifting of making a plot is hidden
                        # in make_plot. We call it once using the LO MC
                        make_plot(acc,
                                region=region,
                                distribution=distribution,
                                year=year,
                                data=data[region],
                                mc=mc_lo[region],
                                ylim=plotset[distribution].get('ylim',None),
                                xlim=plotset[distribution].get('xlim',None),
                                tag = 'losf',
                                outdir=f'./output/{os.path.basename(indir)}/{region}',
                                output_format='pdf',
                                ratio=ratio)

                        # And then we also call it for the NLO MC
                        # The output files will be named according to the 'tag'
                        # argument, so we  will be able to tell them apart.
                        make_plot(acc,
                                region=region,
                                distribution=distribution,
                                year=year,
                                data=data[region],
                                mc=mc_nlo[region],
                                ylim=plotset[distribution].get('ylim',None),
                                xlim=plotset[distribution].get('xlim',None),
                                tag = 'nlo',
                                outdir=f'./output/{os.path.basename(indir)}/{region}',
                                output_format='pdf',
                                ratio=ratio)
                   
                    except KeyError:
                        continue

def main():
    inpath = sys.argv[1]
    plot(inpath)


if __name__ == "__main__":
    main()
