import numpy as np
from bucoffea.helpers import dasgowrapper

def get_datasets():
    datasets = {}
    for line in np.loadtxt("../datasets_2016.txt",dtype=str):
        name, dataset = line.split(" ")
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