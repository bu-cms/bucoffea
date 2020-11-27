"""Customized version of the Coffea executor submodule"""

from __future__ import print_function, division
import concurrent.futures
from functools import partial
from itertools import repeat
import time
import uproot
import uproot4
import pickle
import sys
import math
import numpy as np
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
DEFAULT_METADATA_CACHE = LRUCache(100000)

# instrument xrootd source
if not hasattr(uproot.source.xrootd.XRootDSource, '_read_real'):
    def _read(self, chunkindex):
        self.bytesread = getattr(self, 'bytesread', 0) + self._chunkbytes
        return self._read_real(chunkindex)

    uproot.source.xrootd.XRootDSource._read_real = uproot.source.xrootd.XRootDSource._read
    uproot.source.xrootd.XRootDSource._read = _read

def _work_function_nanoaod(item, processor_instance, flatten=False, savemetrics=False,
                   mmap=False, cachestrategy=None, schema=None, skipbadfiles=False,
                   retries=0, xrootdtimeout=None):
    if processor_instance == 'heavy':
        item, processor_instance = item
    if not isinstance(processor_instance, ProcessorABC):
        processor_instance = cloudpickle.loads(lz4f.decompress(processor_instance))

    import warnings
    out = processor_instance.accumulator.identity()
    retry_count = 0
    while retry_count <= retries:
        try:
            filecontext = uproot4.open(
                item.filename,
                timeout=xrootdtimeout,
                file_handler=uproot4.MemmapSource if mmap else uproot4.MultithreadedFileSource,
            )
            with filecontext as file:
                # To deprecate
                tree = file[item.treename]
                events = LazyDataFrame(tree, item.entrystart, item.entrystop, flatten=flatten)

                # For NanoAOD, we have to look at the "Runs" TTree for info such as weight sums
                # The different cases in the loop represent the different formats and accordingly
                # different ways of dealing with the provided values.
                for name in file['Runs'].keys():
                    if name.startswith('n'):
                        arr = file['Runs'][name].array()
                        # Check that all instances are the same, then save that value
                        tmp = set([])
                        for entry in arr:
                            tmp.add(entry)
                        assert(len(tmp)==1)
                        events[name] = list(tmp)[0]
                    elif any([x in name for x in ['genEventSumw','genEventSumw2']]):
                        arr = file['Runs'][name].array()
                        # One entry per run -> just sum
                        events[name] = int(item.entrystart==0) * np.sum(arr)
                    elif any([x in name for x in ['LHEScaleSumw','LHEPdfSumw']]):
                        # Sum per variation, conserve number of variations
                        arr = file['Runs'][name].array()
                        tmp = np.zeros(len(arr[0]))
                        for i in range(len(arr)):
                            for j in range(len(arr[i])):
                                tmp[j] += arr[i][j]
                        events[name] = int(item.entrystart==0) * tmp
                        pass

                    events['dataset'] = item.dataset
                    events['filename'] = item.filename
                tic = time.time()
                try:
                    out = processor_instance.process(events)
                except Exception as e:
                    raise Exception(f"Failed processing file: {item.filename} ({item.entrystart}-{item.entrystop})") from e
                toc = time.time()
                metrics = dict_accumulator()
                if savemetrics:
                    if isinstance(file, uproot4.ReadOnlyDirectory):
                        metrics['bytesread'] = value_accumulator(int, file.file.source.num_requested_bytes)
                    metrics['columns'] = set_accumulator(events.materialized)
                    metrics['entries'] = value_accumulator(int, events.size)
                    metrics['processtime'] = value_accumulator(float, toc - tic)
                return dict_accumulator({'out': out, 'metrics': metrics})
            break
        # catch xrootd errors and optionally skip
        # or retry to read the file
        except OSError as e:
            if not skipbadfiles:
                raise e
            else:
                w_str = 'Bad file source %s.' % item.filename
                if retries:
                    w_str += ' Attempt %d of %d.' % (retry_count + 1, retries + 1)
                    if retry_count + 1 < retries:
                        w_str += ' Will retry.'
                    else:
                        w_str += ' Skipping.'
                else:
                    w_str += ' Skipping.'
                warnings.warn(w_str)
            metrics = dict_accumulator()
            if savemetrics:
                metrics['bytesread'] = value_accumulator(int, 0)
                metrics['dataservers'] = set_accumulator({})
                metrics['columns'] = set_accumulator({})
                metrics['entries'] = value_accumulator(int, 0)
                metrics['processtime'] = value_accumulator(float, 0)
            wrapped_out = dict_accumulator({'out': out, 'metrics': metrics})
        except Exception as e:
            if retries == retry_count:
                raise e
            w_str = 'Attempt %d of %d. Will retry.' % (retry_count + 1, retries + 1)
            warnings.warn(w_str)
        retry_count += 1

    return wrapped_out

def run_uproot_job_nanoaod(fileset,
                   treename,
                   processor_instance,
                   executor,
                   executor_args={},
                   pre_executor=None,
                   pre_args=None,
                   chunksize=100000,
                   maxchunks=None,
                   metadata_cache=None,
                   ):
    '''A tool to run a processor using uproot for data delivery
    A convenience wrapper to submit jobs for a file set, which is a
    dictionary of dataset: [file list] entries.  Supports only uproot TTree
    reading, via NanoEvents or LazyDataFrame.  For more customized processing,
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
            Some options are not passed to executors but rather and affect the behavior of the
            work function itself:
            - ``savemetrics`` saves some detailed metrics for xrootd processing (default False)
            - ``schema`` builds the dataframe as a `NanoEvents` object rather than `LazyDataFrame`
              (default ``None``); schema options include `NanoEvents`, `NanoAODSchema` and `TreeMakerSchema`
            - ``processor_compression`` sets the compression level used to send processor instance to workers (default 1)
            - ``skipbadfiles`` instead of failing on a bad file, skip it (default False)
            - ``retries`` optionally retry processing of a chunk on failure (default 0)
            - ``xrootdtimeout`` timeout for xrootd read (seconds)
            - ``tailtimeout`` timeout requirement on job tails (seconds)
            - ``align_clusters`` aligns the chunks to natural boundaries in the ROOT files (default False)
        pre_executor : callable
            A function like executor, used to calculate fileset metadata
            Defaults to executor
        pre_args : dict, optional
            Similar to executor_args, defaults to executor_args
        chunksize : int, optional
            Maximum number of entries to process at a time in the data frame, default: 100k
        maxchunks : int, optional
            Maximum number of chunks to process per dataset
            Defaults to processing the whole dataset
        metadata_cache : mapping, optional
            A dict-like object to use as a cache for (file, tree) metadata that is used to
            determine chunking.  Defaults to a in-memory LRU cache that holds 100k entries
            (about 1MB depending on the length of filenames, etc.)  If you edit an input file
            (please don't) during a session, the session can be restarted to clear the cache.
    '''

    import warnings

    if not isinstance(fileset, (Mapping, str)):
        raise ValueError("Expected fileset to be a mapping dataset: list(files) or filename")
    if not isinstance(processor_instance, ProcessorABC):
        raise ValueError("Expected processor_instance to derive from ProcessorABC")

    # make a copy since we modify in-place
    executor_args = dict(executor_args)

    if pre_executor is None:
        pre_executor = executor
    if pre_args is None:
        pre_args = dict(executor_args)
    else:
        pre_args = dict(pre_args)
    if metadata_cache is None:
        metadata_cache = DEFAULT_METADATA_CACHE

    fileset = list(_normalize_fileset(fileset, treename))
    for filemeta in fileset:
        filemeta.maybe_populate(metadata_cache)

    # pop _get_metdata args here (also sent to _work_function)
    skipbadfiles = executor_args.pop('skipbadfiles', False)
    if executor is dask_executor:
        # this executor has a builtin retry mechanism
        retries = 0
    else:
        retries = executor_args.pop('retries', 0)
    xrootdtimeout = executor_args.pop('xrootdtimeout', None)
    align_clusters = executor_args.pop('align_clusters', False)
    metadata_fetcher = partial(_get_metadata,
                               skipbadfiles=skipbadfiles,
                               retries=retries,
                               xrootdtimeout=xrootdtimeout,
                               align_clusters=align_clusters,
                               )

    chunks = []
    if maxchunks is None:
        # this is a bit of an abuse of map-reduce but ok
        to_get = set(filemeta for filemeta in fileset if not filemeta.populated(clusters=align_clusters))
        if len(to_get) > 0:
            out = set_accumulator()
            pre_arg_override = {
                'desc': 'Preprocessing',
                'unit': 'file',
                'compression': None,
                'tailtimeout': None,
                'worker_affinity': False,
            }
            pre_args.update(pre_arg_override)
            pre_executor(to_get, metadata_fetcher, out, **pre_args)
            while out:
                item = out.pop()
                metadata_cache[item] = item.metadata
            for filemeta in fileset:
                filemeta.maybe_populate(metadata_cache)
        while fileset:
            filemeta = fileset.pop()
            if skipbadfiles and not filemeta.populated(clusters=align_clusters):
                continue
            for chunk in filemeta.chunks(chunksize, align_clusters):
                chunks.append(chunk)
    else:
        # get just enough file info to compute chunking
        nchunks = defaultdict(int)
        while fileset:
            filemeta = fileset.pop()
            if nchunks[filemeta.dataset] >= maxchunks:
                continue
            if skipbadfiles and not filemeta.populated(clusters=align_clusters):
                continue
            if not filemeta.populated(clusters=align_clusters):
                filemeta.metadata = metadata_fetcher(filemeta).pop().metadata
                metadata_cache[filemeta] = filemeta.metadata
            for chunk in filemeta.chunks(chunksize, align_clusters):
                chunks.append(chunk)
                nchunks[filemeta.dataset] += 1
                if nchunks[filemeta.dataset] >= maxchunks:
                    break

    # pop all _work_function args here
    savemetrics = executor_args.pop('savemetrics', False)
    if "flatten" in executor_args:
        warnings.warn("Executor argument 'flatten' is deprecated, please refactor your processor to accept awkward arrays", DeprecationWarning)
    flatten = executor_args.pop('flatten', False)
    mmap = executor_args.pop('mmap', False)
    schema = executor_args.pop('schema', None)
    nano = executor_args.pop('nano', False)
 
    cachestrategy = executor_args.pop('cachestrategy', None)
    pi_compression = executor_args.pop('processor_compression', 1)
    if pi_compression is None:
        pi_to_send = processor_instance
    else:
        pi_to_send = lz4f.compress(cloudpickle.dumps(processor_instance), compression_level=pi_compression)
    closure = partial(
        _work_function_nanoaod,
        flatten=flatten,
        savemetrics=savemetrics,
        mmap=mmap,
        schema=schema,
        cachestrategy=cachestrategy,
        skipbadfiles=skipbadfiles,
        retries=retries,
        xrootdtimeout=xrootdtimeout,
    )
    # hack around dask/dask#5503 which is really a silly request but here we are
    if executor is dask_executor:
        executor_args['heavy_input'] = pi_to_send
        closure = partial(closure, processor_instance='heavy')
    else:
        closure = partial(closure, processor_instance=pi_to_send)

    out = processor_instance.accumulator.identity()
    wrapped_out = dict_accumulator({'out': out, 'metrics': dict_accumulator()})
    exe_args = {
        'unit': 'chunk',
        'function_name': type(processor_instance).__name__,
    }
    exe_args.update(executor_args)
    executor(chunks, closure, wrapped_out, **exe_args)

    wrapped_out['metrics']['chunks'] = value_accumulator(int, len(chunks))
    processor_instance.postprocess(out)
    if savemetrics:
        return out, wrapped_out['metrics']
    return out