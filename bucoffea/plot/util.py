import hashlib
import os
import re
from collections import defaultdict
from pprint import pprint

import numpy as np
from coffea import hist
from coffea.processor.accumulator import dict_accumulator
from coffea.util import load, save
from matplotlib import pyplot as plt
from tqdm import tqdm

from bucoffea.execute.dataset_definitions import short_name
from bucoffea.helpers.dataset import extract_year, is_data
from bucoffea.helpers.paths import bucoffea_path

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
        to_merge = list(map(load,files))

        # Remove first two items from list,
        # merge them and insert in the back
        while len(to_merge) > 1:
            t.update()
            x = to_merge.pop(0)
            y = to_merge.pop(0)
            to_merge.append(x+y)
        t.update()
        assert(len(to_merge)==1)

        acc = to_merge[0]
        save(acc, cache)
        return acc



def merge_extensions(histogram, acc, reweight_pu=True):
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
        m = re.match('.*(_ext\d+).*', d)
        base = d
        if m:
            base = d.replace(m.groups()[0],"")
        mapping[base].append(d)
        if not is_data(d):
            sumw[base] += acc['sumw'][d]
            sumw_pileup[base] += acc['sumw_pileup'][d]
            nevents[base] += acc['nevents'][d]

    histogram = histogram.group(hist.Cat("dataset", "Primary dataset"), "dataset", mapping)
    histogram.scale({k:1/v for k, v in sumw.items()}, axis='dataset')

    pu_renorm = { k : nevents[k] / sumw_pileup[k] for k in sumw_pileup.keys()}

    if reweight_pu:
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
        'SingleMuon_2017' : [x for x in all_datasets if re.match('SingleMuon_2017[A-Z]+',x)],
        'EGamma_2017' : [x for x in all_datasets if re.match('SingleElectron_2017[A-Z]+',x) or re.match('SinglePhoton_2017[A-Z]+',x)],
        'MET_2017' : [x for x in all_datasets if re.match('MET_2017[A-Z]+',x)],
        'JetHT_2017' : [x for x in all_datasets if re.match('JetHT_2017[A-Z]+',x)],

        'SingleMuon_2018' : [x for x in all_datasets if re.match('SingleMuon_2018[A-Z]+',x)],
        'EGamma_2018' : [x for x in all_datasets if re.match('EGamma_2018[A-Z]+',x)],
        'MET_2018' : [x for x in all_datasets if re.match('MET_2018[A-Z]+',x)],
        'JetHT_2018' : [x for x in all_datasets if re.match('JetHT_2018[A-Z]+',x)],

        'GJets_HT_MLM_2017' : [x for x in all_datasets if re.match('GJets_HT-(\d+)To(\d+)-MLM_2017',x)],
        'GJets_HT_MLM_2018' : [x for x in all_datasets if re.match('GJets_HT-(\d+)To(\d+)-MLM_2018',x)],

        'WNJetsToLNu_LHEWpT-FXFX_2017' : [x for x in all_datasets if re.match('W(\d+)JetsToLNu_LHEWpT_(\d+)-.*-FXFX_2017',x)],
        'WNJetsToLNu-FXFX_2018' : [x for x in all_datasets if re.match('WJetsToLNu_(\d+)J-amcatnloFXFX_2018',x)],

        'DYNJetsToLL_M-50_LHEZpT-FXFX_2017' : [x for x in all_datasets if re.match('DY(\d+)JetsToLL_M-50_LHEZpT_(\d+)-.*-FXFX_2017',x)],
        'DYNJetsToLL_M-50_LHEZpT-FXFX_2018' : [x for x in all_datasets if re.match('DY(\d+)JetsToLL_M-50_LHEZpT_(\d+)-.*-FXFX_2018',x)],

        'TTJets-FXFX_2017' : [x for x in all_datasets if re.match('TTJets-amcatnloFXFX_2017',x)],
        'TTJets-FXFX_2018' : [x for x in all_datasets if re.match('TTJets-amcatnloFXFX_2018',x)],

        'ST_2017' : [x for x in all_datasets if re.match('ST.*2017',x)],
        'ST_2018' : [x for x in all_datasets if re.match('ST.*2018',x)],

        'DYNJetsToLL_M-50-MLM_2017' : [x for x in all_datasets if re.match('DY(\d+)JetsToLL_M-50-MLM_2017',x)],
        'DYNJetsToLL_M-50-MLM_2018' : [x for x in all_datasets if re.match('DY(\d+)JetsToLL_M-50-MLM_2018',x)],
        'DYJetsToLL_M-50_HT_MLM_2017' : [x for x in all_datasets if re.match('DYJetsToLL_M-50_HT-(\d+)to.*-MLM_2017',x)],
        'DYJetsToLL_M-50_HT_MLM_2018' : [x for x in all_datasets if re.match('DYJetsToLL_M-50_HT-(\d+)to.*-MLM_2018',x)],
        'WJetsToLNu_HT_MLM_2017' : [x for x in all_datasets if re.match('WJetsToLNu_HT-(\d+)To.*-MLM_2017',x)],
        'WJetsToLNu_HT_MLM_2018' : [x for x in all_datasets if re.match('WJetsToLNu_HT-(\d+)To.*-MLM_2018',x)],
        # 'ZJetsToNuNu_HT_2017' : [x for x in all_datasets if re.match('ZJetsToNuNu_HT-(\d+)To(\d+)-mg_2017',x)],
        # 'ZJetsToNuNu_HT_2018' : [x for x in all_datasets if re.match('ZJetsToNuNu_HT-(\d+)To(\d+)-mg_2018',x)],
        'WNJetsToLNu-MLM_2017' : [x for x in all_datasets if re.match('W(\d+)JetsToLNu_2017',x)],
        'WNJetsToLNu-MLM_2018' : [x for x in all_datasets if re.match('W(\d+)JetsToLNu_2018',x)],
    }

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
    histogram = histogram.group(hist.Cat("dataset", "Primary dataset"), "dataset", mapping)

    return histogram


def load_xs():
    """Function to read per-sample cross sections from fileself.

    :return: Mapping dataset -> cross-section
    :rtype: dict
    """
    xsraw = np.loadtxt(bucoffea_path('data/datasets/xs/xs.txt'),dtype=str)
    xs = {}

    # Convert from CMS-style dataset names to short names
    for full, val, _, _ in xsraw:
        xs[short_name(full)] = float(val)

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
        return 41.3

def scale_xs_lumi(histogram):
    """MC normalization so that it's ready to compare to data

    :param histogram: Histogram to normalize
    :type histogram: coffea Hist
    """
    # Get the list of datasets and filter MC data sets
    datasets = list(map(str, histogram.axis('dataset').identifiers()))

    mcs = [x for x in datasets if not is_data(x)]

    # Normalize to XS * lumi/ sumw
    xs = load_xs()
    norm_dict = {mc : 1e3 * xs[mc] * lumi(extract_year(mc)) for mc in mcs}
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
