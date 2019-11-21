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

pjoin = os.path.join

def datasets(year):
    data = {
                    'cr_1m_v' : f'MET_{year}',
                    'cr_2m_v' : f'MET_{year}',
                    'cr_1e_v' : f'EGamma_{year}',
                    'cr_2e_v' : f'EGamma_{year}',
                    'cr_g_v' : f'EGamma_{year}',
                    # 'sr_v' : f'MET_{year}',
                    'sr_v' : f'nomatch',
                }
    tmp = {}
    for k, v in data.items():
        tmp[k] = re.compile(v)
    data.update(tmp)


 
    mc = {
                'cr_1m_v' : re.compile(f'(TTJets.*FXFX.*|Diboson.*|ST.*|QCD_HT.*|.*DYJetsToLL_M-50_HT_MLM.*|.*WJet.*HT.*).*{year}'),
                'cr_1e_v' : re.compile(f'(TTJets.*FXFX.*|Diboson.*|ST.*|QCD_HT.*|.*DYJetsToLL_M-50_HT_MLM.*|.*WJet.*HT.*).*{year}'),
                'cr_2m_v' : re.compile(f'(TTJets.*FXFX.*|Diboson.*|ST.*|QCD_HT.*|.*DYJetsToLL_M-50_HT_MLM.*).*{year}'),
                'cr_2e_v' : re.compile(f'(TTJets.*FXFX.*|Diboson.*|ST.*|QCD_HT.*|.*DYJetsToLL_M-50_HT_MLM.*).*{year}'),
                'cr_g_v' : re.compile(f'(GJets.*HT.*|QCD_HT.*|W.*HT.*).*{year}'),
                'sr_v' : re.compile(f'(.*WJ.*HT.*|.*ZJetsToNuNu.*HT.*|W.*HT.*|TTJets.*FXFX.*|Diboson.*|QCD_HT.*).*{year}'),
            }
    return data, mc


def legacy_dataset_name(dataset):
    patterns = {
        '.*DY.*' : 'zll',
        'QCD.*' : 'qcd',
        'TT.*' : 'top',
        'Diboson.*' : 'diboson',
        '(MET|EGamma).*' : 'data',
        'WN?J.*' : 'wjets',
        'WJ.*' : 'wjets',
        'ZJ.*' : 'zjets',
        'GJets.*HT' : 'gjets',
        #'.*(Hinv|HToInvisible).*' : 'signal',
        'WH.*Hinv.*' : 'wh',
        'ZH.*HToInvisible.*' : 'zh',
        'VBF.*HToInvisible.*' : 'vbf',
        'GluGlu.*HToInvisible.*' : 'ggh',
        'ggZH*HToInvisible.*' : 'ggzh',
    }

    for pat, ret in patterns.items():
        if re.match(pat, dataset):
            return ret
    raise RuntimeError(f'Cannot find legacy region name for dataset :"{dataset}"')

def legacy_region_name(region):
    patterns = {
        'cr_2m_.*' : 'Zmm',
        'cr_2e_.*' : 'Zee',
        'cr_1m_.*' : 'Wmn',
        'cr_1e_.*' : 'Wen',
        'cr_g_.*' : 'gjets',
        'sr_.*' : 'signal',
    }

    for pat, ret in patterns.items():
        if re.match(pat, region):
            return ret
    raise RuntimeError(f'Cannot find legacy region name for region :"{region}"')
def recoil_bins_2016():
    return [250,300,350,400,500,600,750,1000]


def legacy_limit_input(acc, outdir='./output'):
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
        year = 2017
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
