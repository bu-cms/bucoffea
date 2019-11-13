"""Customized version of the Coffea executor submodule"""

from __future__ import print_function, division
import concurrent.futures
from functools import partial
from itertools import repeat
import time
import uproot
import pickle
import sys
import math
import copy
import cloudpickle
from tqdm.auto import tqdm
from collections import defaultdict
from cachetools import LRUCache
import lz4.frame as lz4f
from coffea.processor import ProcessorABC
from coffea.processor.accumulator import (
    AccumulatorABC,
    value_accumulator,
    set_accumulator,
    dict_accumulator,
)
from coffea.processor.dataframe import (
    LazyDataFrame,
)
from coffea.processor.executor import _normalize_fileset, _get_metadata, dask_executor
try:
    from collections.abc import Mapping, Sequence
except ImportError:
    from collections import Mapping, Sequence


_PICKLE_PROTOCOL = pickle.HIGHEST_PROTOCOL


# instrument xrootd source
if not hasattr(uproot.source.xrootd.XRootDSource, '_read_real'):
    def _read(self, chunkindex):
        self.bytesread = getattr(self, 'bytesread', 0) + self._chunkbytes
        return self._read_real(chunkindex)

    uproot.source.xrootd.XRootDSource._read_real = uproot.source.xrootd.XRootDSource._read
    uproot.source.xrootd.XRootDSource._read = _read

def _work_function_nanoaod(item, processor_instance, flatten=False, savemetrics=False, mmap=False):
    if processor_instance == 'heavy':
        item, processor_instance = item
    if not isinstance(processor_instance, ProcessorABC):
        processor_instance = cloudpickle.loads(lz4f.decompress(processor_instance))
    if mmap:
        localsource = {}
    else:
        opts = dict(uproot.FileSource.defaults)
        opts.update({'parallel': None})

        def localsource(path):
            return uproot.FileSource(path, **opts)

    file = uproot.open(item.filename, localsource=localsource)
    tree = file[item.treename]
    df = LazyDataFrame(tree, item.chunksize, item.index, flatten=flatten)
    # For NanoAOD, we have to look at the "Runs" TTree for info such as weight sums
    # The different cases in the loop represent the different formats and accordingly
    # different ways of dealing with the provided values.
    for name in map(lambda x: x.decode('utf-8'), file['Runs'].keys()):
        arr = file['Runs'][name].array()
        if name.startswith('n'):
            # Check that all instances are the same, then save that value
            tmp = set([])
            for entry in arr:
                tmp.add(entry)
            assert(len(tmp)==1)
            df[name] = list(tmp)[0]
        elif name in ['genEventCount','genEventSumw','genEventSumw2']:
            # One entry per run -> just sum
            df[name] = int(item.index==0) * arr.sum()
        elif name in ['LHEScaleSumw','LHEPdfSumw']:
            # Sum per variation, conserve number of variations
            tmp = 0 * arr[0]
            for i in range(len(arr)):
                for j in range(len(arr[i])):
                    tmp[j] += arr[i][j]
            df[name] = int(item.index==0) * tmp
    ### END NANOAOD
    df['dataset'] = item.dataset
    tic = time.time()
    out = processor_instance.process(df)
    toc = time.time()
    metrics = dict_accumulator()
    if savemetrics:
        if isinstance(file.source, uproot.source.xrootd.XRootDSource):
            metrics['bytesread'] = value_accumulator(int, file.source.bytesread)
            metrics['dataservers'] = set_accumulator({file.source._source.get_property('DataServer')})
        metrics['columns'] = set_accumulator(df.materialized)
        metrics['entries'] = value_accumulator(int, df.size)
        metrics['processtime'] = value_accumulator(float, toc - tic)
    wrapped_out = dict_accumulator({'out': out, 'metrics': metrics})
    file.source.close()
    return wrapped_out

def run_uproot_job_nanoaod(fileset,
                   treename,
                   processor_instance,
                   executor,
                   executor_args={},
                   pre_executor=None,
                   pre_args=None,
                   chunksize=200000,
                   maxchunks=None,
                   metadata_cache=LRUCache(100000)
                   ):
    '''A tool to run a processor using uproot for data delivery
    A convenience wrapper to submit jobs for a file set, which is a
    dictionary of dataset: [file list] entries.  Supports only uproot
    reading, via the LazyDataFrame class.  For more customized processing,
    e.g. to read other objects from the files and pass them into data frames,
    one can write a similar function in their user code.
    Parameters
    ----------
        fileset : dict
            A dictionary ``{dataset: [file, file], }``
            Optionally, if some files' tree name differ, the dictionary can be specified:
            ``{dataset: {'treename': 'name', 'files': [file, file]}, }``
        treename : str
            name of tree inside each root file, can be ``None``;
            treename can also be defined in fileset, which will override the passed treename
        processor_instance : ProcessorABC
            An instance of a class deriving from ProcessorABC
        executor : callable
            A function that takes 3 arguments: items, function, accumulator
            and performs some action equivalent to:
            ``for item in items: accumulator += function(item)``
        executor_args : dict, optional
            Arguments to pass to executor.  See `iterative_executor`,
            `futures_executor`, `dask_executor`, or `parsl_executor` for available options.
            Some options that affect the behavior of this function:
            'savemetrics' saves some detailed metrics for xrootd processing (default False);
            'flatten' removes any jagged structure from the input files (default False);
            'processor_compression' sets the compression level used to send processor instance
            to workers (default 1).
        pre_executor : callable
            A function like executor, used to calculate fileset metadata
            Defaults to executor
        pre_args : dict, optional
            Similar to executor_args, defaults to executor_args
        chunksize : int, optional
            Maximum number of entries to process at a time in the data frame.
        maxchunks : int, optional
            Maximum number of chunks to process per dataset
            Defaults to processing the whole dataset
        metadata_cache : mapping, optional
            A dict-like object to use as a cache for (file, tree) metadata that is used to
            determine chunking.  Defaults to a in-memory LRU cache that holds 100k entries
            (about 1MB depending on the length of filenames, etc.)  If you edit an input file
            (please don't) during a session, the session can be restarted to clear the cache.
    '''
    if not isinstance(fileset, Mapping):
        raise ValueError("Expected fileset to be a mapping dataset: list(files)")
    if not isinstance(processor_instance, ProcessorABC):
        raise ValueError("Expected processor_instance to derive from ProcessorABC")

    if pre_executor is None:
        pre_executor = executor
    if pre_args is None:
        pre_args = executor_args

    fileset = list(_normalize_fileset(fileset, treename))
    for filemeta in fileset:
        filemeta.maybe_populate(metadata_cache)

    chunks = []
    if maxchunks is None:
        # this is a bit of an abuse of map-reduce but ok
        to_get = set(filemeta for filemeta in fileset if not filemeta.populated)
        if len(to_get) > 0:
            out = set_accumulator()
            real_pre_args = {
                'desc': 'Preprocessing',
                'unit': 'file',
            }
            real_pre_args.update(pre_args)
            executor(to_get, _get_metadata, out, **real_pre_args)
            while out:
                item = out.pop()
                metadata_cache[item] = item.numentries
            for filemeta in fileset:
                filemeta.maybe_populate(metadata_cache)
        while fileset:
            filemeta = fileset.pop()
            for chunk in filemeta.chunks(chunksize):
                chunks.append(chunk)
    else:
        # get just enough file info to compute chunking
        nchunks = defaultdict(int)
        while fileset:
            filemeta = fileset.pop()
            if nchunks[filemeta.dataset] >= maxchunks:
                continue
            if not filemeta.populated:
                filemeta.numentries = _get_metadata(filemeta).pop().numentries
                metadata_cache[filemeta] = filemeta.numentries
            for chunk in filemeta.chunks(chunksize):
                chunks.append(chunk)
                nchunks[filemeta.dataset] += 1
                if nchunks[filemeta.dataset] >= maxchunks:
                    break

    # pop all _work_function args here
    savemetrics = executor_args.pop('savemetrics', False)
    flatten = executor_args.pop('flatten', False)
    mmap = executor_args.pop('mmap', False)
    pi_compression = executor_args.pop('processor_compression', 1)
    if pi_compression is None:
        pi_to_send = processor_instance
    else:
        pi_to_send = lz4f.compress(cloudpickle.dumps(processor_instance), compression_level=pi_compression)
    # hack around dask/dask#5503 which is really a silly request but here we are
    if executor is dask_executor:
        executor_args['heavy_input'] = pi_to_send
        closure = partial(_work_function_nanoaod,
                          processor_instance='heavy',
                          flatten=flatten,
                          savemetrics=savemetrics,
                          mmap=mmap,
                          )
    else:
        closure = partial(_work_function_nanoaod,
                          processor_instance=pi_to_send,
                          flatten=flatten,
                          savemetrics=savemetrics,
                          mmap=mmap,
                          )

    out = processor_instance.accumulator.identity()
    wrapped_out = dict_accumulator({'out': out, 'metrics': dict_accumulator()})
    exe_args = {
        'unit': 'chunk',
    }
    exe_args.update(executor_args)
    executor(chunks, closure, wrapped_out, **exe_args)
    wrapped_out['metrics']['chunks'] = value_accumulator(int, len(chunks))
    processor_instance.postprocess(out)
    if savemetrics:
        return out, wrapped_out['metrics']
    return out
