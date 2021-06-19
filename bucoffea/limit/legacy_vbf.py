#!/usr/bin/env python
import copy
import os
import re
import numpy as np
from collections import defaultdict

import uproot
from coffea import hist

import ROOT as r
from bucoffea.plot.util import merge_datasets, merge_extensions, scale_xs_lumi, URTH1
from legacy_monojet import suppress_negative_bins
pjoin = os.path.join

def datasets(year, unblind=False):

    data = {
                    'cr_1m_vbf' : f'MET_{year}',
                    'cr_2m_vbf' : f'MET_{year}',
                    'cr_1e_vbf' : f'EGamma_{year}',
                    'cr_2e_vbf' : f'EGamma_{year}',
                    'cr_g_vbf' : f'EGamma_{year}',
                    'sr_vbf_no_veto_all' : f'nomatch',
                }
    if unblind:
        data['sr_vbf'] = f'MET_{year}'
    mc = {
            'sr_vbf_no_veto_all' : re.compile(f'(WH_WToQQ_Hinv_M125.*|ZH_ZToQQ_HToInv.*M125.*|(VBF|GluGlu)_HToInvisible.*M125.*|ggZH.*|ZJetsToNuNu.*|EW.*|Top_FXFX.*|Diboson.*|QCD_HT.*|DYJetsToLL.*|WJetsToLNu.*HT.*).*{year}'),
            'cr_1m_vbf' : re.compile(f'(EW.*|Top_FXFX.*|Diboson.*|QCD_HT.*|.*DYJetsToLL_M-50_HT_MLM.*|WJetsToLNu.*HT.*).*{year}'),
            'cr_1e_vbf' : re.compile(f'(EW.*|Top_FXFX.*|Diboson.*|QCD_HT.*|.*DYJetsToLL_M-50_HT_MLM.*|WJetsToLNu.*HT.*).*{year}'),
            'cr_2m_vbf' : re.compile(f'(EW.*|Top_FXFX.*|Diboson.*|QCD_HT.*|.*DYJetsToLL_M-50_HT_MLM.*).*{year}'),
            'cr_2e_vbf' : re.compile(f'(EW.*|Top_FXFX.*|Diboson.*|QCD_HT.*|.*DYJetsToLL_M-50_HT_MLM.*).*{year}'),
            'cr_g_vbf' : re.compile(f'(GJets_((HT|DR-0p4)|SM).*|QCD_data.*|WJetsToLNu.*HT.*).*{year}'),
            'sr_vbf' : re.compile('nomatch')
          }

    tmp = {}

    for k, v in data.items():
        tmp[k] = re.compile(v)
    data.update(tmp)

    return data, mc

def legacy_dataset_name_vbf(dataset):
    # For the VBF samples, we have two processes:
    # VBF with dipole recoil ON / OFF
    # The default one (as provided in the function argument) will be mapped to "vbf"
    m = re.match("VBF_HToInvisible_M(\d+)_withDipoleRecoil(_PSweights)?_pow_pythia8_201[0-9]", dataset)
    if m:
        mh = m.groups()[0]
        if mh=="125":
            return f"vbf"
        else:
            return f"vbf{mh}"

    m = re.match("VBF_HToInvisible_M(\d+)_withDipoleRecoil(_PSweights)?_pow_pythia8_201[0-9]", dataset)
    if m:
        mh = m.groups()[0]
        if mh=="125":
            return f"vbf"
        else:
            return f"vbf{mh}"

    m = re.match("ZH_ZToQQ_HToInvisible_M(\d+)(_PSweights)?_pow_pythia8_201[0-9]", dataset)
    if m:
        mh = m.groups()[0]
        if mh=="125":
            return "zh"
        else:
            return f"zh{mh}"

    m = re.match("WH_WToQQ_Hinv_M(\d+)_201[0-9]", dataset)
    if m:
        mh = m.groups()[0]
        if mh=="125":
            return "wh"
        else:
            return f"wh{mh}"

    m = re.match("GluGlu_HToInvisible_M(\d+)_HiggspTgt190_pow_pythia8_201[0-9]", dataset)
    if m:
        mh = m.groups()[0]
        if mh=="125":
            return "ggh"
        else:
            return f"ggh{mh}"

    m = re.match("ggZH_ZToQQ_HToInvisible_M(\d+)(_PSweights)?_pow_pythia8_201[0-9]", dataset)
    if m:
        mh = m.groups()[0]
        if mh=="125":
            return "ggzh"
        else:
            return f"ggzh{mh}"

    m = re.match("WH_HToInv_JHU_ptH150_201[0-9]", dataset)
    if m:
        return "wh_jhu"
    m = re.match("ZH_HToInv_JHU_ptH150_201[0-9]", dataset)
    if m:
        return "zh_jhu"

    patterns = {
        'EWKZ\d?Jets.*ZToLL.*' : 'ewkzll',
        'EWKZ\d?Jets.*ZToNuNu.*' : 'ewkzjets',
        'EWKW.*' : 'ewkwjets',
        'QCD.*' : 'qcd',
        'Top.*' : 'top',
        'Diboson.*' : 'diboson',
        '(MET|EGamma).*' : 'data',
        'WJetsToLNu.*' : 'qcdwjets',
        'ZJetsToNuNu.*' : 'qcdzjets',
        'DYJets.*' : 'qcdzll',
        'GJets_DR-0p4.*' : 'qcdgjets',
        'GJets.*HT.*' : 'qcdgjets',
        'GJets_SM_5f_EWK.*' : 'ewkgjets',
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

def export_coffea_histogram(h, overflow='over'):
    '''Helper function to: coffea histogram -> (sumw, xedges) with the desired overflow behavior.'''
    if h.dim() != 1:
        raise RuntimeError('The dimension of the histogram must be 1.')
    
    sumw, sumw2 = h.values(overflow=overflow, sumw2=True)[()]
    xedges = h.axis('mjj').edges()

    # Add the contents of the overflow to the last bin
    if overflow == 'over':
        sumw[-2] += sumw[-1]
        sumw2[-2] += sumw2[-1]
        sumw = np.r_[0, sumw[:-1], 0]
        sumw2 = np.r_[0, sumw2[:-1], 0]

    return URTH1(edges=xedges, sumw=sumw, sumw2=sumw2)

def legacy_limit_input_vbf(acc, outdir='./output', unblind=False):
    """Writes ROOT TH1s to file as a limit input

    :param acc: Accumulator (processor output)
    :type acc: coffea.processor.accumulator
    :param outdir: Output directory
    :type outdir: string
    """
    distribution = 'mjj'

    regions = [
                'cr_2m_vbf',
                'cr_1m_vbf',
                'cr_2e_vbf',
                'cr_1e_vbf',
                'cr_g_vbf',
                'sr_vbf_no_veto_all'
                ]
    if unblind:
        regions.append("sr_vbf")

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # Rebin
    h = copy.deepcopy(acc[distribution])
    newax = hist.Bin('mjj','$M_{jj}$ (GeV)', mjj_bins_2016())
    h = h.rebin(h.axis(newax.name), newax)
    h = merge_extensions(h, acc)
    scale_xs_lumi(h)
    h = merge_datasets(h)

    for year in [2017,2018]:
        signal = re.compile(f'VBF_HToInvisible.*{year}')
        f = uproot.recreate(pjoin(outdir, f'legacy_limit_vbf_{year}.root'))
        data, mc = datasets(year, unblind=unblind)
        for region in regions:
            print('='*20)
            print(f'Region {region}')
            print('='*20)
            tag = region.split('_')[0]

            ih = h.integrate(h.axis('region'),region)

            for dataset in map(str, h.axis('dataset').identifiers()):
                if not (data[region].match(dataset) or mc[region].match(dataset)):
                    # Insert dummy data for the signal region
                    if region == 'sr_vbf' and re.match('ZJetsToNuNu.*', dataset) and not unblind:
                        th1 = export_coffea_histogram(ih.integrate('dataset', dataset))
                        histo_name = 'signal_data'
                        f[histo_name] = th1
                        continue
                    else:
                        continue
                print(f"Dataset: {dataset}")

                th1 = export_coffea_histogram(ih.integrate('dataset', dataset))
                try:
                    # Patch for now: We only have the DR sample for 2018
                    # So we modify the function "DR" argument based on which year we're looking at
                    histo_name = f'{legacy_region_name(region)}_{legacy_dataset_name_vbf(dataset)}'
                    print(f'Saved under histogram: {histo_name}')
                except:
                    print(f"Skipping {dataset}")
                    continue

                print('-'*20)
                f[histo_name] = th1
        if not unblind:
            f[f'{legacy_region_name("sr_vbf")}_data'] = f[f'{legacy_region_name("sr_vbf")}_qcdzjets']
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

    outfile = r.TFile(pjoin(outdir, f'legacy_limit_vbf.root'),'RECREATE')
    for year, ifiles in files.items():
        for category, file in ifiles.items():
            subdir = outfile.mkdir(f'category_{category}_{year}')
            infile = r.TFile(file)
            for key in infile.GetListOfKeys():
                print(key)
                h = key.ReadObj().Clone()
                h.SetTitle(h.GetName())
                h.SetDirectory(subdir)
                h.GetXaxis().SetTitle('mjj')
                suppress_negative_bins(h)
                subdir.Write()
