import numpy as np
import re
from bucoffea.helpers import dasgowrapper, bucoffea_path
import os
import yaml
import socket
import subprocess
from collections import defaultdict
pjoin = os.path.join

def short_name(dataset):
    _, name, conditions, _ = dataset.split("/")

    # Remove useless info
    name = name.replace("_TuneCP5","")
    name = name.replace("_TuneCUETP8M1","")
    name = name.replace("_13TeV","")
    name = name.replace("-pythia8","")
    name = name.replace("madgraphMLM","MLM")
    name = name.replace("madgraph","mg")
    name = name.replace("amcnloFXFX","FXFX")
    name = name.replace("powheg","pow")

    # Detect extension
    m=re.match(r".*(ext\d+).*",conditions);
    if m:
        name = name + "_" + m.groups()[0]
    m=re.match(r".*(ver\d+).*",conditions);
    if m:
        name = name + "_" + m.groups()[0]
    if 'new_pmx' in conditions:
        name = name + '_new_pmx'
    if 'RunIISummer16' in conditions:
        name = name + "_2016"
    elif "RunIIFall17" in conditions:
        name = name + "_2017"
    elif 'RunIIAutumn18' in conditions:
        name = name + "_2018"

    m = re.match(r"Run(\d+[A-Z]*)", conditions)
    if m:
        name = name + "_" + m.groups()[0]

    return name



def load_lists():
    files = [
        bucoffea_path(f"data/datasets/datasets_2016.txt"),
        bucoffea_path(f"data/datasets/datasets_2017.txt"),
        bucoffea_path(f"data/datasets/datasets_2018.txt")
    ]
    lines = []
    for fpath in files:
        with open(fpath,"r") as f:
            lines.extend(f.readlines())

    lines = filter(lambda x: "NANOAOD" in x and not x.startswith("#"), lines)
    return lines

def files_from_das(regex):
    """Generate file list per dataset from DAS

    :param regex: Regular expression to match datasets
    :type regex: string
    :return: Mapping of dataset : [files]
    :rtype: dict
    """
    datasets = {}

    for line in load_lists():
        # Skip empty lines
        dataset = line.strip()
        if not len(line) or line.startswith("#") or not "/" in line:
            continue
        # Dataset 'nicknames'
        name = short_name(dataset)
        if not re.match(regex, name):
            continue

        # Get files from DAS
        files = dasgowrapper.das_go_query(f"file dataset={dataset}")
        datasets[name] = files.split()

    for dataset, filelist in datasets.items():
            newlist = []
            for file in filelist:
                file = file.decode("utf-8")
                if file.startswith("/store/"):
                    newlist.append("root://cms-xrd-global.cern.ch//" + file)
                else:
                    newlist.append(file)
            datasets[dataset] = newlist

    return datasets

def files_from_ac(regex):
    """Generate file list per dataset from T2_DE_RWTH

    :param regex: Regular expression to match datasets
    :type regex: string
    :return: Mapping of dataset : [files]
    :rtype: dict
    """
    path = bucoffea_path('data/datasets/crabfiles.yml')

    with open(path, 'r') as stream:
        try:
            fileset = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)

    for dataset, files in fileset.items():
        if not re.match(regex, dataset):
            continue
        for ifile in reversed(files):
            if not len(ifile):
                files.remove(ifile)
        fileset[dataset] = files
    return fileset

def find_files(directory, regex):
    fileset = {}
    for path, _, files in os.walk(directory):

        files = list(filter(lambda x: x.endswith('.root'), files))
        if not len(files):
            continue
        dataset = path.split('/')[-3]
        if not re.match(regex, dataset):
            continue
        files = [pjoin(path,x) for x in files]
        fileset[dataset] = files
    return fileset

def eosls(path):
    return subprocess.check_output(['xrdfs', 'root://cmseos.fnal.gov','ls','-l',path]).decode('utf-8')

def eosfind(path):
    cmd = ['eos', 'root://cmseos.fnal.gov/', 'find',  '--size', path]
    return subprocess.check_output(cmd).decode('utf-8')


def find_files_eos(directory, regex):
    fileset = defaultdict(list)
    lines = eosfind(re.sub('.*/store/','/store/',directory)).splitlines()
    # For files, lines are formatted as
    # path=(File path starting with /eos/uscms) size=(Size in bits)
    # For folders, the 'size' part is left out, so they can easily be filtered
    for line in lines:
        parts = line.split()

        # Ignore lines representing directories
        if len(parts) < 2:
            continue
        # Ensure we are not missing a part
        if len(parts) > 2:
            raise RuntimeError(f'Encountered malformed line: {line}')

        # The info we care about
        path = parts[0].replace('path=','')
        if not path.endswith('.root'):
            continue

        dataset = path.split('/')[9]
        if not re.match(regex, dataset):
            continue
        fileset[dataset].append(re.sub('.*/store','root://cmsxrootd-site.fnal.gov//store', path))
    return fileset

def files_from_eos(regex):
    """Generate file list per dataset from EOS

    :param regex: Regular expression to match datasets
    :type regex: string
    :return: Mapping of dataset : [files]
    :rtype: dict
    """

    host = socket.gethostname()
    if 'lxplus' in host:
        topdir = '/eos/cms/store/group/phys_exotica/monojet/aalbert/nanopost/'
        tag = '16Jul19'

        fileset_16jul = find_files(pjoin(topdir, tag), regex)

        topdir = '/eos/user/a/aalbert/nanopost/'
        tag = '10Aug19'

        fileset_10aug = find_files(pjoin(topdir, tag), regex)

        fileset = {}
        keys = set(list(fileset_16jul.keys()) + list(fileset_10aug.keys()))
        for key in keys:
            if key in fileset_10aug:
                fileset[key] = fileset_10aug[key]
            else:
                fileset[key] = fileset_16jul[key]
    elif 'lpc' in host:
        topdir = '/eos/uscms/store/user/aandreas/nanopost/'
        tag = '05Jun20v5'
        fileset = find_files_eos(pjoin(topdir, tag), regex)

    return fileset
