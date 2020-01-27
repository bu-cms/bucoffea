import hashlib
import os
import shutil
import random
import re
import string
from collections import defaultdict
from pprint import pprint

import numpy as np
import yaml
from coffea import hist
from coffea.processor.accumulator import dict_accumulator
from coffea.util import load, save
from matplotlib import pyplot as plt
from tqdm import tqdm

from bucoffea.execute.dataset_definitions import short_name
from bucoffea.helpers.dataset import extract_year, is_data
from bucoffea.helpers.paths import bucoffea_path

from klepto.archives import dir_archive

pjoin = os.path.join

def sha256sum(filelist):
    h  = hashlib.sha256()
    b  = bytearray(128*1024)
    mv = memoryview(b)
    for filename in filelist:
        with open(filename, 'rb', buffering=0) as f:
            for n in iter(lambda : f.readinto(mv), 0):
                h.update(mv[:n])
    return h.hexdigest()

def klepto_load(inpath):
    acc = dir_archive(
                    inpath,
                    serialized=True,
                    compression=0,
                    memsize=1e3,
                    )
    return acc

def acc_from_dir(indir):
    """Load Coffea accumulator from directory with *.coffea files

    :param indir: Directory to search for coffea files
    :type indir: string
    :return: Sum of all found accumulators
    :rtype: dict
    """
    files = filter(lambda x: x.endswith(".coffea") and not ('cache' in x), os.listdir(indir))
    files = list(map(lambda x: os.path.abspath(pjoin(indir, x)), files))
    listhash = sha256sum(files)
    cache = pjoin(indir, f'merged_cache_{listhash}.coffea')
    if os.path.exists(cache):
        return load(cache)
    else:
        # Progress bar
        t = tqdm(total=len(files), desc='Merging input files')

        # Recursive merging
        to_merge = files

        # Use temporary files to store intermediate
        # merger results
        tmp_files = []
        def load_and_remove(path):
            data = load(path)
            os.remove(path)
            return data


        def next():
            '''Get next item to merge'''
            x = to_merge.pop(0)
            if isinstance(x, str):
                if x in tmp_files:
                    tmp_files.remove(x)
                    x = load_and_remove(x)
                else:
                    x = load(x)
            return x

        while len(to_merge) > 1:
            # Remove first two items from list,
            # merge them and insert in the back
            t.update()

            x = next()
            y = next()

            tmp = "/tmp/tmp_bucoffea_merge_" + "".join(random.sample(string.ascii_uppercase+string.digits,24))
            merged = x+y
            # clean up to save memory
            x = None
            y = None
            save(merged, tmp)
            merged = None
            to_merge.append(tmp)
            tmp_files.append(tmp)

        t.update()
        assert(len(to_merge)==1)

        shutil.copy(to_merge[0], cache)
        return load(cache)



def merge_extensions(histogram, acc, reweight_pu=True, noscale=False):
    """Merge extension datasets into one and scale to 1/sumw

    :param histogram: The histogram to modify
    :type histogram: Coffea histogram
    :param acc: The accumulator dictionary holding sumw etc
    :type acc: dict_accumulator
    :param reweight_pu: Whether to renormalize to account for the sum of PU weights
    :type reweight_pu: bool (default: True)
    :return: Modified histogram
    :rtype: Coffea histogram
    """
    all_datasets = map(str, histogram.identifiers('dataset'))
    mapping = defaultdict(list)
    sumw = defaultdict(float)
    sumw_pileup = defaultdict(float)
    nevents = defaultdict(float)

    for d in all_datasets:
        base = d

        to_replace =[
            '.*(_ext\d+).*',
            '.*(_new_+pmx).*',
            '.*(_PSweights).*'
        ]
        for regex in to_replace:
            m = re.match(regex, base)
            if m:
                base = base.replace(m.groups()[0],"")

        mapping[base].append(d)
        if not is_data(d):
            sumw[base] += acc['sumw'][d]
            if reweight_pu:
                sumw_pileup[base] += acc['sumw_pileup'][d]
                nevents[base] += acc['nevents'][d]
    histogram = histogram.group("dataset", hist.Cat("dataset", "Primary dataset"), mapping)

    if not noscale:
        histogram.scale({k:1/v for k, v in sumw.items()}, axis='dataset')

        if reweight_pu:
            pu_renorm = { k : nevents[k] / sumw_pileup[k] for k in sumw_pileup.keys()}
            histogram.scale(pu_renorm, axis='dataset')

    return histogram


def merge_datasets(histogram):
    """Merge datasets that belong same physics process

    :param histogram: The histogram to modify
    :type histogram: Coffea histogram
    :return: Modified histogram
    :rtype: Coffea histogram
    """
    all_datasets = list(map(str, histogram.identifiers('dataset')))
    # TODO:
    #   * Factor mapping out to configuration file?
    #   * Fill in more data sets
    #   * lots of duplicate code (re.match etc) -> simplify
    mapping = {
        'SingleMuon_2016' : [x for x in all_datasets if re.match('SingleMuon_2016[A-Z]+',x)],
        'EGamma_2016' : [x for x in all_datasets if re.match('SingleElectron_.*2016[A-Z]+',x) or re.match('SinglePhoton_2016[A-Z]+',x)],
        'MET_2016' : [x for x in all_datasets if re.match('MET_.*2016[A-Z]+',x)],
        'JetHT_2016' : [x for x in all_datasets if re.match('JetHT_.*2016[A-Z]+',x)],

        'SingleMuon_2017' : [x for x in all_datasets if re.match('SingleMuon_.*2017[A-Z]+',x)],
        'EGamma_2017' : [x for x in all_datasets if re.match('SingleElectron_.*2017[A-Z]+',x) or re.match('SinglePhoton_2017[A-Z]+',x)],
        'MET_2017' : [x for x in all_datasets if re.match('MET_.*2017[A-Z]+',x)],
        'JetHT_2017' : [x for x in all_datasets if re.match('JetHT_.*2017[A-Z]+',x)],

        'SingleMuon_2018' : [x for x in all_datasets if re.match('SingleMuon_.*2018[A-Z]+',x)],
        'EGamma_2018' : [x for x in all_datasets if re.match('EGamma_.*2018[A-Z]+',x)],
        'MET_2018' : [x for x in all_datasets if re.match('MET_.*2018[A-Z]+',x)],
        'JetHT_2018' : [x for x in all_datasets if re.match('JetHT_.*2018[A-Z]+',x)],

        'GJets_SM_5f_EWK-mg_2017' : ['GJets_SM_5f_EWK-mg_2017'],
        'GJets_SM_5f_EWK-mg_2018' : ['GJets_SM_5f_EWK-mg_2017'],
        
        'G1Jet_Pt-amcatnlo_2016' : [x for x in all_datasets if re.match('G1Jet_Pt-.*-amcatnlo_2016',x)],

        'WNJetsToLNu_LHEWpT-FXFX_2017' : [x for x in all_datasets if re.match('W(\d+)JetsToLNu_LHEWpT_(\d+)-.*-FXFX_2017',x)],
        'WNJetsToLNu-FXFX_2018' : [x for x in all_datasets if re.match('WJetsToLNu_(\d+)J-amcatnloFXFX_2018',x)],

        'DYNJetsToLL_M-50_LHEZpT-FXFX_2017' : [x for x in all_datasets if re.match('DY(\d+)JetsToLL_M-50_LHEZpT_(\d+)-.*-FXFX_2017',x)],
        'DYNJetsToLL_M-50_LHEZpT-FXFX_2018' : [x for x in all_datasets if re.match('DY(\d+)JetsToLL_M-50_LHEZpT_(\d+)-.*-FXFX_2018',x)],

        'DYNJetsToLL_M-50-MLM_2017' : [x for x in all_datasets if re.match('DY(\d+)JetsToLL_M-50-MLM_2017',x)],
        'DYNJetsToLL_M-50-MLM_2018' : [x for x in all_datasets if re.match('DY(\d+)JetsToLL_M-50-MLM_2018',x)],

        'ZJetsToNuNu_HT_2017' : [x for x in all_datasets if re.match('ZJetsToNuNu_HT-(\d+)To.*-mg_2017',x)],
        'ZJetsToNuNu_HT_2018' : [x for x in all_datasets if re.match('ZJetsToNuNu_HT-(\d+)To.*-mg_2018',x)],

        'WNJetsToLNu-MLM_2017' : [x for x in all_datasets if re.match('W(\d+)JetsToLNu_2017',x)],
        'WNJetsToLNu-MLM_2018' : [x for x in all_datasets if re.match('W(\d+)JetsToLNu_2018',x)],

        'WH_WToQQ_Hinv_M125_2017' : [x for x in all_datasets if re.match('W.*H_WToQQ_HToInvisible_M125.*2017',x)],
        'WH_WToQQ_Hinv_M125_2018' : [x for x in all_datasets if re.match('W.*H_WToQQ_HToInvisible_M125.*2018',x)]
    }

    # Some combinations are the same for all years
    yearly = {
        'GJets_HT_MLM_{year}' : 'GJets_HT-(\d+)To.*-MLM_{year}',
        'GJets_DR-0p4_HT_MLM_{year}' : 'GJets_DR-0p4_HT-(\d+)To.*-MLM_.*{year}',
        'WJetsToQQ_HT_MLM_{year}' : 'WJetsToQQ_HT-?(\d+)(T|t)o.*-MLM_{year}',
        'DYJetsToLL_M-50_HT_MLM_{year}' : 'DYJetsToLL_M-50_HT-(\d+)to.*-MLM_{year}',
        'WJetsToLNu_HT_MLM_{year}' : 'WJetsToLNu_HT-(\d+)To.*-MLM_{year}',

        'TTJets-FXFX_{year}' : 'TTJets-amcatnloFXFX_{year}',
        'TTJets-MLM_{year}' : 'TTJets-MLM_{year}',
        'TT_pow_{year}' : 'TTTo.*pow.*{year}',
        'ST_{year}' : 'ST.*{year}',

        'QCD_HT_{year}' : 'QCD_HT.*_{year}',

        'EWKW2Jets_WToLNu_M-50-mg_{year}' : 'EWKW(Plus|Minus)2Jets.*-mg_{year}',
        'Diboson_{year}' : '(WW|WZ|ZZ|WW).*_{year}',

    }
    for year in [2016,2017,2018]:
        for name, regex in yearly.items():
            mapping[name.format(year=year)] = [x for x in all_datasets if re.match(regex.format(year=year), x)]

    # Remove empty lists
    tmp = {}
    for k, v in mapping.items():
        if len(v):
            tmp[k] = v
    mapping = tmp

    # Add datasets we didn't catch yet
    mapped_datasets =  []
    for val in mapping.values():
        mapped_datasets.extend(val)

    for ds in all_datasets:
        if ds in mapped_datasets:
            continue
        else:
            mapping[ds] = [ds]
    # Apply the mapping
    histogram = histogram.group("dataset",hist.Cat("dataset", "Primary dataset"),  mapping)

    return histogram


def load_xs():
    """Function to read per-sample cross sections from file.

    :return: Mapping dataset -> cross-section
    :rtype: dict
    """
    xsfile = bucoffea_path('data/datasets/xs/xs.yml')
    with open(xsfile,'r') as f:
        xs_yml = yaml.load(f, Loader=yaml.FullLoader)

    # See the documentation in data/datasets/xs/README.md
    # for more information on how the XS is chosen.
    xs = {}
    loading_priority = ['nnnlo','nnlo','nlo','lo','gen']
    for dataset, xs_dict in xs_yml.items():
        if 'use' in xs_dict:
            key_to_use = xs_dict['use']
        else:
            for key in loading_priority:
                if key in xs_dict:
                    key_to_use = key
                    break
        xs[dataset] = xs_dict[key_to_use]

    # Data sets that only exist as extensions
    # cause problems later on, so we duplicate the XS
    # for the base process.
    tmp = {}
    for k in xs.keys():
        base = re.sub('_ext(\d+)','',k)
        if base not in xs.keys():
            tmp[base] = xs[k]

    xs.update(tmp)
    return xs

def lumi(year):
    """Golden JSON luminosity per for given year

    :param year: Year of data taking
    :type year: int
    :return: Golden JSON luminosity for that year in pb (!)
    :rtype: float
    """
    if year==2018:
        return 59.7
    if year==2017:
        return 41.5
    if year==2016:
        return 35.9

def scale_xs_lumi(histogram):
    """MC normalization so that it's ready to compare to data

    :param histogram: Histogram to normalize
    :type histogram: coffea Hist
    """
    # Get the list of datasets and filter MC data sets
    datasets = list(map(str, histogram.axis('dataset').identifiers()))

    mcs = [x for x in datasets if not is_data(x)]

    # Normalize to XS * lumi/ sumw
    known_xs = load_xs()

    xs_map = {}
    for mc in mcs:
        try:
            ixs = known_xs[re.sub('_new_*pmx','',mc)]
        except KeyError:
            print(f"WARNING: Cross section not found for dataset {mc}. Using 0.")
            ixs = 0
        xs_map[mc] = ixs

    norm_dict = {mc : 1e3 * xs_map[mc] * lumi(extract_year(mc)) for mc in mcs}
    histogram.scale(norm_dict, axis='dataset')

# def merge_and_norm(histogram, acc):
#     histogram = merge_extensions(histogram, acc)
#     scale_xs_lumi(histogram)
#     histogram = merge_datasets(histogram)
#     return histogram

def fig_ratio():
    """Shortcut to create figure with ratio and main panels

    :return: Figure and axes for main and ratio panels
    :rtype: tuple(Figure, axes, axes)
    """
    fig, (ax, rax) = plt.subplots(2, 1, figsize=(7,7), gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
    return fig, ax, rax
