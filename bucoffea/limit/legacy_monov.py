#!/usr/bin/env python
import copy
import os
import re
from collections import defaultdict

import uproot
from coffea import hist
from coffea.hist.export import export1d

import ROOT as r
from bucoffea.plot.util import merge_datasets, merge_extensions, scale_xs_lumi
from legacy_monojet import legacy_dataset_name, datasets, legacy_region_name, suppress_negative_bins

pjoin = os.path.join

def recoil_bins_2016():
    return [250,300,350,400,500,600,750,1000]


def legacy_limit_input_monov(acc, outdir='./output', unblind=False):
    """Writes ROOT TH1s to file as a limit input

    :param acc: Accumulator (processor output)
    :type acc: coffea.processor.accumulator
    :param outdir: Output directory
    :type outdir: string
    """
    distribution = 'recoil'

    regions = [
                'cr_2m_v',
                'cr_1m_v',
                'cr_2e_v',
                'cr_1e_v',
                'cr_g_v',
                'sr_v_no_veto_all'
                ]
    if unblind:
        regions.append("sr_v")

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    for wp in ['tau21','loosemd','tightmd','loose','tight']:
        for year in [2017,2018]:
            signal = re.compile(f'.*(Hinv|HToInvisible).*{year}')
            f = uproot.recreate(pjoin(outdir, f'legacy_limit_monov_{wp}_{year}.root'))
            data, mc = datasets(year, unblind)
            for region in regions:
                if wp == 'tau21':
                    monov_region_name = region
                else:
                    if region.endswith("_v"):
                        monov_region_name = region.replace('_v',f'_{wp}_v')
                    else:
                        monov_region_name = region.replace('_v_',f'_{wp}_v_')
                print(f'Region {region}')
                # Rebin
                h = copy.deepcopy(acc[distribution])
                
                newax = hist.Bin('recoil','Recoil (GeV)', recoil_bins_2016())

                h = h.rebin(h.axis(newax.name), newax)

                h = merge_extensions(h, acc)
                scale_xs_lumi(h)

                h = merge_datasets(h)

                h = h.integrate(h.axis('region'),monov_region_name)
                
                for dataset in map(str, h.axis('dataset').identifiers()):
                    if not (data[region].match(dataset) or mc[region].match(dataset) or signal.match(dataset)):
                        continue
                    print(f"   Dataset: {dataset}")

                    th1 = export1d(h.integrate('dataset', dataset))
                    try:
                        histo_name = f'{legacy_region_name(region)}_{legacy_dataset_name(dataset)}'
                    except:
                        print(f"Skipping {dataset}")
                        continue
                    f[histo_name] = th1
            if not unblind:
                f[f'{legacy_region_name("sr_v")}_data'] = f[f'{legacy_region_name("sr_v")}_zjets']
    merge_legacy_inputs(outdir)


def merge_legacy_inputs(outdir):
    '''
    Workaround for uproot's lack of subdirectory support.
    '''

    files = defaultdict(dict)
    for fname in os.listdir(outdir):
        m = re.match('legacy_limit_monov_([0-9a-z]*)_(\d+).root', fname)
        if not m:
            continue
        wp, year = m.groups()
        files[wp][year] = pjoin(outdir, fname)

    outfile_nominal = r.TFile(pjoin(outdir, f'merged_legacy_limit_nominal_monov.root'),'RECREATE')
    outfile_MD = r.TFile(pjoin(outdir, f'merged_legacy_limit_MD_monov.root'),'RECREATE')
    for wp, ifiles in files.items():
        outfile = r.TFile(pjoin(outdir, f'merged_legacy_limit_monov_{wp}.root'),'RECREATE')
        for year, file in ifiles.items():
            subdir = outfile.mkdir(f'category_monov_{year}')
            infile = r.TFile(file)
            for key in infile.GetListOfKeys():
                print(key)
                h = key.ReadObj().Clone()
                h.SetTitle(h.GetName())
                h.SetDirectory(subdir)
                h.GetXaxis().SetTitle('met')
                suppress_negative_bins(h)
                subdir.Write()

            if wp == 'tau21':
                continue
            elif 'md' in wp:
                outfile_MD.cd()
                subdir = outfile_MD.mkdir(f'category_monov{wp.replace("md","")}_{year}')
            else:
                outfile_nominal.cd()
                subdir = outfile_nominal.mkdir(f'category_monov{wp}_{year}')
            infile = r.TFile(file)
            for key in infile.GetListOfKeys():
                print(key)
                h = key.ReadObj().Clone()
                h.SetTitle(h.GetName())
                h.SetDirectory(subdir)
                h.GetXaxis().SetTitle('met')
                suppress_negative_bins(h)
                subdir.Write()
