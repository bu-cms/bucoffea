import numpy as np
from bucoffea.helpers import dasgowrapper

def get_datasets():
    datasets = {}
    for dataset in np.loadtxt("../datasets_2016.txt",dtype=str):
        files = dasgowrapper.das_go_query(f"file dataset={dataset}")
        datasets[dataset] = files.split()

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