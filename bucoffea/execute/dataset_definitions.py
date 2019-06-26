import numpy as np
from bucoffea.helpers import dasgowrapper

def get_datasets():
    datasets = {}
    with open("datasets_2016.txt","r") as fin:
        lines = fin.readlines()
    for line in lines:
        name, dataset = line.strip().split(" ")
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