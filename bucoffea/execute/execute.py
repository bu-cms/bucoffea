#!/usr/bin/env python
import os
import argparse
from coffea import processor
from bucoffea.processors.monojet import monojetProcessor
import lz4.frame as lz4f
import cloudpickle
import htcondor


pjoin = os.path.join
def do_run(args):
    # Create output directory
    if not os.path.exists(args.outpath):
        os.makedirs(args.outpath)
    
    # Run over all files associated to dataset
    datasets = get_datasets()
    fileset = { args.dataset : datasets[args.dataset]}
    output = processor.run_uproot_job(fileset,
                                  treename='Events',
                                  processor_instance=monojetProcessor(),
                                  executor=processor.futures_executor,
                                  executor_args={'workers': 1, 'function_args': {'flatten': True}},
                                  chunksize=500000,
                                 )

    # Save output
    outname_base = "hists_{}.cpkl.lz4".format(args.dataset)
    with lz4f.open(pjoin(args.outpath, outname_base), mode="wb", compression_level=5) as fout:
        cloudpickle.dump(output, fout)


def do_submit(args):
    datasets = get_datasets()

    schedd = htcondor.Schedd()
    for dataset in datasets.keys():
        sub = htcondor.Submit({
            "executable": "execute.py", 
            "arguments": " --outpath {} run --dataset {}".format(dataset, args.outpath),
            "log" : "/dev/null"
            })
        with schedd.transaction() as txn:
            print(sub.queue(txn))

def main():
    # create the top-level parser
    parser = argparse.ArgumentParser(prog='PROG')
    parser.add_argument('--outpath', type=str, help='Path to save output under.')
    
    subparsers = parser.add_subparsers(help='sub-command help')
    
    parser_run = subparsers.add_parser('run', help='run help')
    parser_run.add_argument('--dataset', type=str, help='Dataset name to run over.')
    parser_run.set_defaults(func=do_run)

    # create the parser for the "b" command
    parser_submit = subparsers.add_parser('submit', help='b help')
    parser_submit.add_argument('--baz', choices='XYZ', help='baz help')
    parser_submit.set_defaults(func=do_submit)


    args = parser.parse_args()

    args.func(args)

from dataset_definitions import get_datasets


if __name__ == "__main__":
    main()