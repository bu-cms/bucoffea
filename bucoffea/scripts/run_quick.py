#!/usr/bin/env python

from bucoffea.helpers.dataset import extract_year
from bucoffea.processor.executor import run_uproot_job_nanoaod
from bucoffea.helpers.cutflow import print_cutflow
from coffea.util import save
import coffea.processor as processor
import argparse

def parse_commandline():

    parser = argparse.ArgumentParser()
    parser.add_argument('processor', type=str, help='The processor to be run. (monojet or vbfhinv)')
    args = parser.parse_args()

    return args

def main():

    fileset = {
        "ZJetsToNuNu_HT-400To600-mg_2018" : [
            "/eos/cms/store/group/phys_exotica/monojet/aalbert/nanopost/sync_27Oct19/ZJetsToNuNu_HT-400To600_13TeV-madgraph/ZJetsToNuNu_HT-400To600-mg_2018/191108_200239/0000/tree_1.root",
            "/eos/cms/store/group/phys_exotica/monojet/aalbert/nanopost/sync_27Oct19/ZJetsToNuNu_HT-400To600_13TeV-madgraph/ZJetsToNuNu_HT-400To600-mg_2018/191108_200239/0000/tree_2.root",
            "/eos/cms/store/group/phys_exotica/monojet/aalbert/nanopost/sync_27Oct19/ZJetsToNuNu_HT-400To600_13TeV-madgraph/ZJetsToNuNu_HT-400To600-mg_2018/191108_200239/0000/tree_3.root",
            "/eos/cms/store/group/phys_exotica/monojet/aalbert/nanopost/sync_27Oct19/ZJetsToNuNu_HT-400To600_13TeV-madgraph/ZJetsToNuNu_HT-400To600-mg_2018/191108_200239/0000/tree_4.root",
            "/eos/cms/store/group/phys_exotica/monojet/aalbert/nanopost/sync_27Oct19/ZJetsToNuNu_HT-400To600_13TeV-madgraph/ZJetsToNuNu_HT-400To600-mg_2018/191108_200239/0000/tree_5.root",
            "/eos/cms/store/group/phys_exotica/monojet/aalbert/nanopost/sync_27Oct19/ZJetsToNuNu_HT-400To600_13TeV-madgraph/ZJetsToNuNu_HT-400To600-mg_2018/191108_200239/0000/tree_6.root",  
            "/eos/cms/store/group/phys_exotica/monojet/aalbert/nanopost/sync_27Oct19/ZJetsToNuNu_HT-400To600_13TeV-madgraph/ZJetsToNuNu_HT-400To600-mg_2018/191108_200239/0000/tree_7.root",
            "/eos/cms/store/group/phys_exotica/monojet/aalbert/nanopost/sync_27Oct19/ZJetsToNuNu_HT-400To600_13TeV-madgraph/ZJetsToNuNu_HT-400To600-mg_2018/191108_200239/0000/tree_8.root",
            "/eos/cms/store/group/phys_exotica/monojet/aalbert/nanopost/sync_27Oct19/ZJetsToNuNu_HT-400To600_13TeV-madgraph/ZJetsToNuNu_HT-400To600-mg_2018/191108_200239/0000/tree_9.root"
        ]
    }

    years = list(set(map(extract_year, fileset.keys())))
    assert(len(years)==1)

    args = parse_commandline()
    processor_class = args.processor

    if processor_class == 'monojet':
        from bucoffea.monojet import monojetProcessor
        processorInstance = monojetProcessor(years[0])
    elif processor_class == 'vbfhinv':
        from bucoffea.vbfhinv import vbfhinvProcessor
        processorInstance = vbfhinvProcessor(years[0])
    elif processor_class == 'lhe':
        from bucoffea.gen.lheVProcessor import lheVProcessor
        processorInstance = lheVProcessor()
    elif args.processor == 'purity':
        from bucoffea.photon_purity import photonPurityProcessor
        processorInstance = photonPurityProcessor()

    for dataset, filelist in fileset.items():
        newlist = []
        for file in filelist:
            if file.startswith("/store/"):
                newlist.append("root://cms-xrd-global.cern.ch//" + file)
            else: newlist.append(file)
        fileset[dataset] = newlist

    for dataset, filelist in fileset.items():
        tmp = {dataset:filelist}
        output = run_uproot_job_nanoaod(tmp,
                                    treename='Events',
                                    processor_instance=processorInstance,
                                    executor=processor.futures_executor,
                                    executor_args={'workers': 4, 'flatten': True},
                                    chunksize=500000,
                                    )
        save(output, f"{processor_class}_{dataset}.coffea")
        # Debugging / testing output
        # debug_plot_output(output)
        print_cutflow(output, outfile=f'{processor_class}_cutflow_{dataset}.txt')

if __name__ == "__main__":
    main()
