#!/usr/bin/env python
import os
import argparse
from coffea import processor
from bucoffea.monojet import monojetProcessor
import lz4.frame as lz4f
import cloudpickle
from pathlib import Path
from bucoffea.execute.dataset_definitions import get_datasets
from coffea.util import save
from bucoffea.helpers import bucoffea_path, vo_proxy_path, xrootd_format
from bucoffea.helpers.condor import  condor_submit
import shutil
from datetime import datetime
pjoin = os.path.join
def do_run(args):
    """Run the analysis locally."""
    # Run over all files associated to dataset
    fileset = get_datasets(regex=args.dataset)

    if "2016" in args.dataset: year=2016
    elif "2017" in args.dataset: year=2017
    elif "2018" in args.dataset: year=2018
    else: raise RuntimeError("Cannot deduce year from dataset name.")

    output = processor.run_uproot_job(fileset,
                                  treename='Events',
                                  processor_instance=monojetProcessor(year=year),
                                  executor=processor.futures_executor,
                                  executor_args={'workers': args.jobs, 'flatten': True},
                                  chunksize=500000,
                                 )

    # Save output
    outpath = pjoin(args.outpath, f"monojet_{args.dataset}.coffea")
    save(output, outpath)

def do_worker(args):
    """Run the analysis on a worker node."""
    # Run over all files associated to dataset

    with open(args.filelist, "r") as f:
        files = [xrootd_format(x.strip()) for x in f.readlines()]

    fileset = { args.dataset : files}
    if "2016" in args.dataset: year=2016
    elif "2017" in args.dataset: year=2017
    elif "2018" in args.dataset: year=2018
    else: raise RuntimeError("Cannot deduce year from dataset name.")

    output = processor.run_uproot_job(fileset,
                                  treename='Events',
                                  processor_instance=monojetProcessor(year=year),
                                  executor=processor.futures_executor,
                                  executor_args={'workers': args.jobs, 'flatten': True},
                                  chunksize=500000,
                                 )

    # Save output
    outpath = pjoin(args.outpath, f"monojet_{args.dataset}_{args.chunk}.coffea")
    save(output, outpath)



def chunkify(items, nchunk):
    '''Split list of items into nchunk ~equal sized chunks'''
    chunks = [[] for _ in range(nchunk)]
    for i in range(len(items)):
        chunks[i % nchunk].append(items[i])
    return chunks

def do_submit(args):
    """Submit the analysis to HTCondor."""
    import htcondor

    dataset_files = get_datasets(regex=args.dataset)

    subdir = os.path.abspath("./submission/")
    if not os.path.exists(subdir):
        os.makedirs(subdir)
    filedir = 'files'

    if not os.path.exists(pjoin(subdir, filedir)):
        os.makedirs(pjoin(subdir, filedir))

    # Get proxy and copy to a safe location on AFS
    proxy = vo_proxy_path()
    proxydir = os.path.expanduser("~/.voms/")
    if not os.path.exists(proxydir):
        os.makedirs(proxydir)
    shutil.copy2(proxy, proxydir)

    # schedd = htcondor.Schedd()
    for dataset, files in dataset_files.items():
        print(f"Submitting dataset: {dataset}.")

        chunks = chunkify(files,int(len(files) / args.filesperjob + 1) )
        for ichunk, chunk in enumerate(chunks):
            # Save input files to a txt file and send to job
            tmpfile = pjoin(subdir, filedir, f"input_{dataset}_{ichunk}.txt")
            with open(tmpfile, "w") as f:
                for file in chunk:
                    f.write(f"{file}\n")

            arguments = [
                # pjoin(proxydir, os.path.basename(proxy)),
                "$(Proxy_path)",
                str(Path(__file__).absolute()),
                f'--outpath {args.outpath}',
                f'--jobs {args.jobs}',
                'worker',
                f'--dataset {dataset}',
                f'--filelist {os.path.basename(tmpfile)}',
                f'--chunk {ichunk}'
            ]
            input_files = [
                # bucoffea_path("config.yaml"),
                # os.path.join(proxydir, os.path.basename(proxy)),
                os.path.abspath(tmpfile),
            ]
            sub = htcondor.Submit({
                "Proxy_path" : pjoin(proxydir,os.path.basename(proxy)),
                "Initialdir" : subdir,
                "executable": bucoffea_path("execute/htcondor_wrap.sh"),
                "should_transfer_files" : "YES",
                "when_to_transfer_output" : "ON_EXIT",
                "transfer_input_files" : ", ".join(input_files),
                "getenv" : "true",
                "arguments": " ".join(arguments),
                "Output" : f"{filedir}/out_{dataset}_{ichunk}.txt",
                "Error" : f"{filedir}/err_{dataset}_{ichunk}.txt",
                "log" :f"/dev/null",
                "request_cpus" : str(args.jobs),
                "+MaxRuntime" : f"{60*60*8}"
                })

            jdl = pjoin(subdir,f'job_{dataset}_{ichunk}.jdl')
            with open(jdl,"w") as f:
                f.write(str(sub))
                f.write("\nqueue 1")
            # with schedd.transaction() as txn:
                # print(sub.queue(txn))
            jobid = condor_submit(jdl)
            print(f"Submitted job {jobid}")
            with open("submission_history.txt","a") as f:
                f.write(f"{datetime.now()} {jobid}\n")
        break

def main():
    parser = argparse.ArgumentParser(prog='Execution wrapper for coffea analysis')
    parser.add_argument('--outpath', type=str, help='Path to save output under.')
    parser.add_argument('--jobs','-j', type=int, default=1, help='Number of cores to use / request.')

    subparsers = parser.add_subparsers(help='sub-command help')

    # Arguments passed to the "run" operation
    parser_run = subparsers.add_parser('run', help='Running help')
    parser_run.add_argument('--dataset', type=str, help='Dataset name to run over.')
    parser_run.set_defaults(func=do_run)

    # Arguments passed to the "worker" operation
    parser_run = subparsers.add_parser('worker', help='Running help')
    parser_run.add_argument('--dataset', type=str, help='Dataset name to run over.')
    parser_run.add_argument('--filelist', type=str, help='Text file with file names to run over.')
    parser_run.add_argument('--chunk', type=str, help='Number of this chunk for book keeping.')
    parser_run.set_defaults(func=do_worker)

    # Arguments passed to the "submit" operation
    parser_submit = subparsers.add_parser('submit', help='Submission help')
    parser_submit.add_argument('--dataset', type=str, help='Dataset regex to use.')
    parser_submit.add_argument('--filesperjob', type=int, default=10, help='Number of files to process per job')
    parser_submit.set_defaults(func=do_submit)


    args = parser.parse_args()

    # Create output directory
    args.outpath = os.path.abspath(args.outpath)
    if not os.path.exists(args.outpath):
        os.makedirs(args.outpath)

    args.func(args)



if __name__ == "__main__":
    main()