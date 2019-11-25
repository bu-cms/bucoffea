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
                    'cr_1m_vbf' : f'MET_{year}',
                    'cr_2m_vbf' : f'MET_{year}',
                    'cr_1e_vbf' : f'EGamma_{year}',
                    'cr_2e_vbf' : f'EGamma_{year}',
                    'cr_g_vbf' : f'EGamma_{year}',
                    # 'sr_vbf' : f'MET_{year}',
                    'sr_vbf' : f'nomatch',
                }
    mc = {
            'sr_vbf' : re.compile(f'W(minus|plus)H_.*|((VBF|GluGlu)_HToInvisible.*|ggZH.*|ZJetsToNuNu.*|EW.*|TTJets.*FXFX.*|Diboson.*|ST.*|QCD_HT.*|.*DYJetsToLL_M-50_HT_MLM.*|.*W.*HT.*).*{year}'),
            'cr_1m_vbf' : re.compile(f'(EW.*|TTJets.*FXFX.*|Diboson.*|ST.*|QCD_HT.*|.*DYJetsToLL_M-50_HT_MLM.*|.*W.*HT.*).*{year}'),
            'cr_1e_vbf' : re.compile(f'(EW.*|TTJets.*FXFX.*|Diboson.*|ST.*|QCD_HT.*|.*DYJetsToLL_M-50_HT_MLM.*|.*W.*HT.*).*{year}'),
            'cr_2m_vbf' : re.compile(f'(EW.*|TTJets.*FXFX.*|Diboson.*|ST.*|QCD_HT.*|.*DYJetsToLL_M-50_HT_MLM.*).*{year}'),
            'cr_2e_vbf' : re.compile(f'(EW.*|TTJets.*FXFX.*|Diboson.*|ST.*|QCD_HT.*|.*DYJetsToLL_M-50_HT_MLM.*).*{year}'),
            'cr_g_vbf' : re.compile(f'(GJets_(?!Mjj).*|AGJets.*|QCD_HT.*|W.*HT.*).*{year}'),
          }
        
    tmp = {}
    
    for k, v in data.items():
        tmp[k] = re.compile(v)
    data.update(tmp)

    return data, mc

def legacy_dataset_name_vbf(dataset):
    patterns = {
        '.*DY.*' : 'qcdzll',
        'EWK_ZToLL.*' : 'ewkzll',
        'EWK_ZToNuNu.*' : 'ewkzjets',
        'EWK_W.*' : 'ewkwjets',
        'QCD.*' : 'qcd',
        'TT.*' : 'top',
        'Diboson.*' : 'diboson',
        '(MET|EGamma).*' : 'data',
        'WN?J.*' : 'qcdwjets',
        'WJ.*' : 'qcdwlnu',
        'ZJ.*' : 'qcdzjets',
        'GJets.*HT' : 'qcdgjets',
        'AGJets.*' : 'ewkgjets',
        'VBF_HToInv.*' : 'vbf',
        'GluGlu_HToInv.*' : 'ggh',
        'ggZH_.*' : 'zh',
        'W(minus|plus)H_.*' : 'wh'
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
    return [ 250.,  280.,  310.,  340.,  370.,  400.,  
             430.,  470.,  510., 550.,  590.,  640.,  
             690.,  740.,  790.,  840.,  900.,  960., 
             1020., 1090., 1160., 1250., 1400.]

def mjj_bins_2016():
    return [200., 400., 600., 900., 1200., 1500.,
            2000., 2750., 3500., 5000.]

def legacy_limit_input_vbf(acc, outdir='./output'):
    """Writes ROOT TH1s to file as a limit input

    :param acc: Accumulator (processor output)
    :type acc: coffea.processor.accumulator
    :param outdir: Output directory
    :type outdir: string
    """
    distribution = 'mjj'

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    for year in [2017,2018]:
        signal = re.compile(f'VBF_HToInvisible.*{year}')
        f = uproot.recreate(pjoin(outdir, f'legacy_limit_vbf_{year}.root'))
        data, mc = datasets(year)
        for region in ['cr_2m_vbf','cr_1m_vbf','cr_2e_vbf','cr_1e_vbf','cr_g_vbf','sr_vbf']:
            print(f'Region {region}')
            # Rebin
            h = copy.deepcopy(acc[distribution])
            
            newax = hist.Bin('mjj','$M_{jj}$ (GeV)', mjj_bins_2016()) 

            h = h.rebin(h.axis(newax.name), newax)

            h = merge_extensions(h, acc)
            scale_xs_lumi(h)

            h = merge_datasets(h)

            h = h.integrate(h.axis('region'),region)
            
            for dataset in map(str, h.axis('dataset').identifiers()):
                if not (data[region].match(dataset) or mc[region].match(dataset) or signal.match(dataset)):
                    print(f"Skip dataset: {dataset}")
                    continue
                print(f"   Dataset: {dataset}")

                th1 = export1d(h.integrate('dataset', dataset))
                try:
                    histo_name = f'{legacy_region_name(region)}_{legacy_dataset_name_vbf(dataset)}'
                    print(histo_name)
                except:
                    print(f"Skipping {dataset}")
                    continue
                f[histo_name] = th1
        #f[f'{legacy_region_name("sr_vbf")}_data'] = f[f'{legacy_region_name("sr_vbf")}_zjets']
    merge_legacy_inputs(outdir)

def merge_legacy_inputs(outdir):
    '''
    Workaround for uproot's lack of subdirectory support.
    '''

    files = defaultdict(dict)
    for fname in os.listdir(outdir):
        m = re.match('legacy_limit_([a-z]*)_(\d+).root', fname)
        if not m:
            continue
        category, year = m.groups()
        files[year][category] = pjoin(outdir, fname)

    for year, ifiles in files.items():
        outfile = r.TFile(pjoin(outdir, f'legacy_limit_{year}.root'),'RECREATE')
        for category, file in ifiles.items():
            subdir = outfile.mkdir(f'category_{category}')
            infile = r.TFile(file)
            for key in infile.GetListOfKeys():
                print(key)
                h = key.ReadObj().Clone()
                h.SetTitle(h.GetName())
                h.SetDirectory(subdir)
                h.GetXaxis().SetTitle('mjj')
                # h.Write()
                subdir.Write()
