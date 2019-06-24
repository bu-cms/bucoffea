#!/usr/bin/env python
import os
import argparse
from coffea import processor
from bucoffea.monojet import monojetProcessor
import lz4.frame as lz4f
import cloudpickle
from pathlib import Path
from dataset_definitions import get_datasets




pjoin = os.path.join
def do_run(args):
    """Run the analysis locally."""
    # Run over all files associated to dataset
    datasets = get_datasets()
    fileset = { args.dataset : datasets[args.dataset]}
    output = processor.run_uproot_job(fileset,
                                  treename='Events',
                                  processor_instance=monojetProcessor(),
                                  executor=processor.futures_executor,
                                  executor_args={'workers': args.jobs, 'function_args': {'flatten': False}},
                                  chunksize=500000,
                                 )

    # Save output
    outname_base = "hists_{}.cpkl.lz4".format(args.dataset)
    with lz4f.open(pjoin(args.outpath, outname_base), mode="wb", compression_level=5) as fout:
        cloudpickle.dump(output, fout)


def do_submit(args):
    """Submit the analysis to HTCondor."""
    import htcondor

    datasets = get_datasets()

    schedd = htcondor.Schedd()
    for dataset in datasets.keys():
        sub = htcondor.Submit({
            "executable": pjoin(str(Path(__file__).absolute().parent), "htcondor_wrap.sh"),
            "getenv" : "true",
            "arguments": "{} --outpath {} run --dataset {} -j {}".format(str(Path(__file__).absolute()), args.outpath, dataset, args.jobs),
            "Output" : pjoin(args.outpath,"out_{}.txt".format(dataset)),
            "Error" : pjoin(args.outpath,"err_{}.txt".format(dataset)),
            "log" :pjoin(args.outpath,"log_{}.txt".format(dataset)),
            "request_cpus" : args.jobs
            })
        with schedd.transaction() as txn:
            print(sub.queue(txn))
        break

def main():
    parser = argparse.ArgumentParser(prog='Execution wrapper for coffea analysis')
    parser.add_argument('--outpath', type=str, help='Path to save output under.')
    parser.add_argument('--jobs','-j', type=int, help='Number of cores to use / request.')

    subparsers = parser.add_subparsers(help='sub-command help')

    # Arguments passed to the "run" operation
    parser_run = subparsers.add_parser('run', help='Running help')
    parser_run.add_argument('--dataset', type=str, help='Dataset name to run over.')
    parser_run.set_defaults(func=do_run)

    # Arguments passed to the "submit" operation
    parser_submit = subparsers.add_parser('submit', help='Submission help')
    # parser_submit.add_argument('--baz', choices='XYZ', help='baz help')
    parser_submit.set_defaults(func=do_submit)


    args = parser.parse_args()

    # Create output directory
    if not os.path.exists(args.outpath):
        os.makedirs(args.outpath)

    args.func(args)



if __name__ == "__main__":
    main()