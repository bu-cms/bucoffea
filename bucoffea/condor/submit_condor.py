import argparse
from coffea import hist
from coffea.analysis_objects import JaggedCandidateArray
import coffea.processor as processor
from awkward import JaggedArray
import numpy as np

import lz4.frame as lz4f
import cloudpickle




def do_run(args):
    datasets = get_datasets()
    fileset = { args.dataset : datasets[args.dataset]}
    output = processor.run_uproot_job(fileset,
                                  treename='Events',
                                  processor_instance=monojetProcessor(),
                                  executor=processor.futures_executor,
                                  executor_args={'workers': 1, 'function_args': {'flatten': True}},
                                  chunksize=500000,
                                 )

    
    with lz4f.open("hists_{}.cpkl.lz4".format(args.dataset), mode="wb", compression_level=5) as fout:
        cloudpickle.dump(output, fout)
    pass

def do_submit(args):
    datasets = get_datasets()
   

    with lz4f.open("hists.cpkl.lz4", mode="wb", compression_level=5) as fout:
        cloudpickle.dump(output, fout)


def main():
    # create the top-level parser
    parser = argparse.ArgumentParser(prog='PROG')
    
    subparsers = parser.add_subparsers(help='sub-command help')
    
    parser_run = subparsers.add_parser('run', help='run help')
    parser_run.add_argument('--dataset', type=str, help='Dataset name to run over.')
    parser_run.add_argument('bar', type=int, help='bar help')
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