#!/usr/bin/env python
import os
import argparse
from coffea import processor
from bucoffea.monojet import monojetProcessor
import lz4.frame as lz4f
import cloudpickle
from pathlib import Path
from dataset_definitions import get_datasets
from coffea.util import save
from bucoffea.helpers import bucoffea_path, vo_proxy_path
import shutil

pjoin = os.path.join
def do_run(args):
    """Run the analysis locally."""
    # Run over all files associated to dataset
    datasets = get_datasets()
    fileset = { args.dataset : datasets[args.dataset]}
    if "2016" in args.dataset: year=2016
    elif "2017" in args.dataset: year=2017
    elif "2018" in args.dataset: year=2018
    else: raise RuntimeError("Cannot deduce year from dataset name.")

    fileset = { args.dataset : datasets[args.dataset][:1]}
    output = processor.run_uproot_job(fileset,
                                  treename='Events',
                                  processor_instance=monojetProcessor(year=year),
                                  executor=processor.futures_executor,
                                  executor_args={'workers': args.jobs, 'function_args': {'flatten': True}},
                                  chunksize=500000,
                                 )

    # Save output
    outpath = pjoin(args.outpath, f"monojet_{args.dataset}.coffea")
    save(output, outpath)


def do_submit(args):
    """Submit the analysis to HTCondor."""
    import htcondor

    datasets = get_datasets()

    subdir = os.path.abspath(".")
    if not os.path.exists(subdir):
        os.makedirs(subdir)

    # Get proxy and copy to a safe location on AFS
    proxy = vo_proxy_path()
    proxydir = os.path.expanduser("~/.voms/")
    if not os.path.exists(proxydir):
        os.makedirs(proxydir)
    shutil.copy2(proxy, proxydir)
    input_files = [
        # pjoin(proxydir, os.path.basename(proxy)),
        bucoffea_path("config.yaml"),

    ]

    schedd = htcondor.Schedd()
    for dataset in datasets.keys():
        print(f"Submitting dataset: {dataset}.")
        sub = htcondor.Submit({
            "Initialdir" : subdir,
            "executable": bucoffea_path("execute/htcondor_wrap.sh"),
            # "input": pjoin(str(Path(__file__).absolute().parent), "htcondor_wrap.sh"),
            "should_transfer_files" : "YES",
            "when_to_transfer_output" : "ON_EXIT",
            "transfer_input_files" : ", ".join(input_files),
            "getenv" : "true",
            "arguments": "{} {} --outpath {} run --dataset {} -j {}".format(pjoin(proxydir, os.path.basename(proxy)), str(Path(__file__).absolute()), args.outpath, dataset, args.jobs),
            "Output" : f"out_{dataset}.txt",
            "Error" : f"err_{dataset}.txt",
            "log" :f"/dev/null",
            "request_cpus" : str(args.jobs),
            "+MaxRuntime" : f"{60*60*8}"
            })
        with open("job.jdl","w") as f:
            f.write(str(sub))
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