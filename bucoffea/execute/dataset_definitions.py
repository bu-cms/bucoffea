import numpy as np
import re
from bucoffea.helpers import dasgowrapper, bucoffea_path

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
    m=re.match(".*(ext\d+).*",conditions);
    if m:
        name = name + "_" + m.groups()[0]

    if "RunIIFall17" in conditions:
        name = name + "_2017"
    elif 'RunIIAutumn18' in conditions:
        name = name + "_2018"
    return name

def load_lists():
    files = [
        bucoffea_path(f"data/datasets/datasets_2017.txt"),
        bucoffea_path(f"data/datasets/datasets_2018.txt")
    ]
    lines = []
    for fpath in files:
        with open(fpath,"r") as f:
            lines.extend(f.readlines())
    return lines

def get_datasets(regex):
    datasets = {}

    for line in load_lists():
        # Skip empty lines
        dataset = line.strip()
        if not len(line):
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