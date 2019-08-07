#!/usr/bin/env python
import argparse
import os
import shutil
from datetime import datetime
from multiprocessing.pool import Pool
from pathlib import Path

import cloudpickle
import lz4.frame as lz4f
from coffea import processor
from coffea.util import save

from bucoffea.execute.dataset_definitions import (files_from_ac,
                                                  files_from_das,
                                                  files_from_eos)
from bucoffea.helpers import bucoffea_path, vo_proxy_path, xrootd_format
from bucoffea.helpers.condor import condor_submit
from bucoffea.helpers.git import git_rev_parse, git_diff
from bucoffea.processor.executor import run_uproot_job_nanoaod

pjoin = os.path.join

def choose_processor(args):
    if args.processor == 'monojet':
        from bucoffea.monojet import monojetProcessor
        return monojetProcessor
    elif args.processor == 'vbfhinv':
        from bucoffea.vbfhinv import vbfhinvProcessor
        return vbfhinvProcessor
    elif args.processor == 'pdfweight':
        from bucoffea.gen import pdfWeightProcessor
        return pdfWeightProcessor

def do_run(args):
    """Run the analysis locally."""
    # Run over all files associated to dataset
    if args.datasrc == 'das':
        fileset = files_from_das(regex=args.dataset)
    else:
        fileset = files_from_eos(regex=args.dataset)

    for dataset, files in fileset.items():
        output = run_uproot_job_nanoaod({dataset:files},
                                    treename='Events',
                                    processor_instance=choose_processor(args)(),
                                    executor=processor.futures_executor,
                                    executor_args={'workers': args.jobs, 'flatten': True},
                                    chunksize=500000,
                                    )

        # Save output
        try:
            os.makedirs(args.outpath)
        except FileExistsError:
            pass
        outpath = pjoin(args.outpath, f"monojet_{dataset}.coffea")
        save(output, outpath)

def do_worker(args):
    """Run the analysis on a worker node."""
    # Run over all files associated to dataset
    with open(args.filelist, "r") as f:
        files = [xrootd_format(x.strip()) for x in f.readlines()]
    fileset = {args.dataset : files}

    ndatasets = len(fileset)
    nfiles = sum([len(x) for x in fileset.values()])
    print(f"Running over {ndatasets} datasets with a total of {nfiles} files.")

    output = run_uproot_job_nanoaod(fileset,
                                  treename='Events',
                                  processor_instance=choose_processor(args)(),
                                  executor=processor.futures_executor,
                                  executor_args={'workers': args.jobs, 'flatten': True},
                                  chunksize=500000,
                                 )

    # Save output
    try:
        os.makedirs(args.outpath)
    except FileExistsError:
        pass
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

    if args.datasrc == 'das':
        dataset_files = files_from_das(regex=args.dataset)
    elif args.datasrc == 'ac':
        dataset_files = files_from_ac(regex=args.dataset)
    else:
        dataset_files = files_from_eos(regex=args.dataset)

    # Test mode: One file per data set
    if args.test:
        tmp = {}
        for k, v in dataset_files.items():
            tmp[k] = v[:1]
        dataset_files = tmp

    # Submission directory:
    # Uses tag from commandline if specified
    # Or just time tag
    timetag = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    if args.name:
        subdir = os.path.abspath(pjoin("./submission/", args.name))
        if os.path.exists(subdir) and not args.force:
            raise RuntimeError(f"Will not overwrite existing task directory unless '--force' is specified: {subdir}")
    else:
        subdir = os.path.abspath(pjoin("./submission/", timetag))
    if not os.path.exists(subdir):
        os.makedirs(subdir)

    # Repo version information
    with open(pjoin(subdir, 'version.txt'),'w') as f:
        f.write(git_rev_parse()+'\n')
        f.write(git_diff()+'\n')

    # Sub-directory to store submission files
    filedir = 'files'
    if not os.path.exists(pjoin(subdir, filedir)):
        os.makedirs(pjoin(subdir, filedir))

    # Get proxy and copy to a safe location on AFS
    proxy = vo_proxy_path()
    proxydir = os.path.expanduser("~/.voms/")
    if not os.path.exists(proxydir):
        os.makedirs(proxydir)
    shutil.copy2(proxy, proxydir)

    if args.async:
        jdl_to_submit = []

    for dataset, files in dataset_files.items():
        print(f"Submitting dataset: {dataset}.")

        chunks = chunkify(files,int(len(files) / args.filesperjob + 1) )
        for ichunk, chunk in enumerate(chunks):
            # Save input files to a txt file and send to job
            tmpfile = pjoin(subdir, filedir, f"input_{dataset}_{ichunk:03d}of{len(chunks):03d}.txt")
            with open(tmpfile, "w") as f:
                for file in chunk:
                    f.write(f"{file}\n")

            # Job file creation
            arguments = [
                # pjoin(proxydir, os.path.basename(proxy)),
                "$(Proxy_path)",
                str(Path(__file__).absolute()),
                args.processor,
                f'--outpath {pjoin(subdir, "output")}',
                f'--jobs {args.jobs}',
                'worker',
                f'--dataset {dataset}',
                f'--filelist {os.path.basename(tmpfile)}',
                f'--chunk {ichunk}'
            ]
            input_files = [
                os.path.abspath(tmpfile),
            ]
            environment = {
                "NOPREFETCH" : str(args.no_prefetch).lower()
            }
            sub = htcondor.Submit({
                "Proxy_path" : pjoin(proxydir,os.path.basename(proxy)),
                "Initialdir" : subdir,
                "executable": bucoffea_path("execute/htcondor_wrap.sh"),
                "should_transfer_files" : "YES",
                "when_to_transfer_output" : "ON_EXIT",
                "transfer_input_files" : ", ".join(input_files),
                "getenv" : "true",
                "environment" : '"' + ' '.join([f"{k}={v}" for k, v in environment.items()]) + '"',
                "arguments": " ".join(arguments),
                "Output" : f"{filedir}/out_{dataset}_{ichunk:03d}of{len(chunks):03d}.txt",
                "Error" : f"{filedir}/err_{dataset}_{ichunk:03d}of{len(chunks):03d}.txt",
                "log" : f"{filedir}/log_{dataset}_{ichunk:03d}of{len(chunks):03d}.txt",
                # "log" :f"/dev/null",
                "request_cpus" : str(args.jobs),
                "+MaxRuntime" : f"{60*60*8}"
                })

            jdl = pjoin(subdir,filedir,f'job_{dataset}_{ichunk}.jdl')
            with open(jdl,"w") as f:
                f.write(str(sub))
                f.write("\nqueue 1\n")

            # Submission
            if args.dry:
                jobid = -1
                print(f"Submitted job {jobid}")
            else:
                if args.async:
                    jdl_to_submit.append(jdl)
                else:
                    jobid = condor_submit(jdl)
                    print(f"Submitted job {jobid}")
    if args.async:
        p = Pool(processes=8)
        res = p.map_async(condor_submit, jdl_to_submit)
        res.wait()
        if res.successful():
            print(f"Asynchronous submission successful for {len(jdl_to_submit)} jobs.")
        else:
            print("Asynchronous submission failed.")


def main():
    parser = argparse.ArgumentParser(prog='Execution wrapper for coffea analysis')
    parser.add_argument('processor', type=str, help='Processor to run.', choices=['monojet','vbfhinv','pdfweight'])
    parser.add_argument('--outpath', type=str, help='Path to save output under.')
    parser.add_argument('--jobs','-j', type=int, default=1, help='Number of cores to use / request.')
    parser.add_argument('--datasrc', type=str, default='eos', help='Source of data files.', choices=['eos','das','ac'])

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
    parser_submit.add_argument('--name', type=str, default=None, help='Name to identify this submission')
    parser_submit.add_argument('--no-prefetch', action="store_true", default=False, help='Do not prefetch input files on worker but run over xrootd.')
    parser_submit.add_argument('--dry', action="store_true", default=False, help='Do not trigger submission, just dry run.')
    parser_submit.add_argument('--test', action="store_true", default=False, help='Only run over one file per dataset for testing.')
    parser_submit.add_argument('--force', action="store_true", default=False, help='Overwrite existing submission folder with same tag.')
    parser_submit.add_argument('--async', action="store_true", default=False, help='Submit asynchronously.')
    parser_submit.set_defaults(func=do_submit)


    args = parser.parse_args()

    args.func(args)



if __name__ == "__main__":
    main()
