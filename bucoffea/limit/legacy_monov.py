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
from legacy_monojet import legacy_dataset_name, datasets, legacy_region_name

pjoin = os.path.join

def recoil_bins_2016():
    return [250,300,350,400,500,600,750,1000]


def legacy_limit_input_monov(acc, outdir='./output'):
    """Writes ROOT TH1s to file as a limit input

    :param acc: Accumulator (processor output)
    :type acc: coffea.processor.accumulator
    :param outdir: Output directory
    :type outdir: string
    """
    distribution = 'recoil'

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    for wp in ['tau21','loosemd','tightmd','loose','tight']:
        for year in [2017,2018]:
            signal = re.compile(f'.*(Hinv|HToInvisible).*{year}')
            f = uproot.recreate(pjoin(outdir, f'legacy_limit_monov_{wp}_{year}.root'))
            data, mc = datasets(year)
            for region in ['cr_2m_v','cr_1m_v','cr_2e_v','cr_1e_v','cr_g_v','sr_v']:
                if wp == 'tau21':
                    monov_region_name = region
                else:
                    monov_region_name = region.replace('_v',f'_{wp}_v')
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
                        print(f"Skip dataset: {dataset}")
                        continue
                    print(f"   Dataset: {dataset}")

                    th1 = export1d(h.integrate('dataset', dataset))
                    try:
                        histo_name = f'{legacy_region_name(region)}_{legacy_dataset_name(dataset)}'
                    except:
                        print(f"Skipping {dataset}")
                        continue
                    f[histo_name] = th1
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
        files[year][wp] = pjoin(outdir, fname)

    for year, ifiles in files.items():
        for wp, file in ifiles.items():
            outfile = r.TFile(pjoin(outdir, f'merged_legacy_limit_monov_{wp}_{year}.root'),'RECREATE')
            subdir = outfile.mkdir(f'category_monov')
            infile = r.TFile(file)
            for key in infile.GetListOfKeys():
                print(key)
                h = key.ReadObj().Clone()
                h.SetTitle(h.GetName())
                h.SetDirectory(subdir)
                h.GetXaxis().SetTitle('met')
                # h.Write()
                subdir.Write()
        # produce a combined version
        outfile_nominal = r.TFile(pjoin(outdir, f'merged_legacy_limit_nominal_monov_{year}.root'),'RECREATE')
        outfile_MD = r.TFile(pjoin(outdir, f'merged_legacy_limit_MD_monov_{year}.root'),'RECREATE')
        for wp, file in ifiles.items():
            if wp == 'tau21':
                continue
            elif 'md' in wp:
                outfile_MD.cd()
                subdir = outfile_MD.mkdir('category_monov'+(wp.replace('md','')))
            else:
                outfile_nominal.cd()
                subdir = outfile_nominal.mkdir('category_monov'+wp)
            infile = r.TFile(file)
            for key in infile.GetListOfKeys():
                print(key)
                h = key.ReadObj().Clone()
                h.SetTitle(h.GetName())
                h.SetDirectory(subdir)
                h.GetXaxis().SetTitle('met')
                # h.Write()
                subdir.Write()
