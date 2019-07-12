"""Customized version of the Coffea executor submodule"""


from tqdm import tqdm
import uproot
from coffea.processor.executor import _get_chunking
from coffea.processor import LazyDataFrame, ProcessorABC
from coffea.processor.accumulator import value_accumulator, set_accumulator, dict_accumulator

try:
    from collections.abc import Mapping, Sequence
    from functools import lru_cache
except ImportError:
    from collections import Mapping, Sequence

    def lru_cache(maxsize):
        def null_wrapper(f):
            return f
        return null_wrapper

def _get_chunking_lazy_noempty(filelist, treename, chunksize):
    for fn in filelist:
        nentries = uproot.numentries(fn, treename)
        if not nentries:
            continue
        for index in range(nentries // chunksize + 1):
            yield (fn, chunksize, index)

def run_uproot_job_nanoaod(fileset, treename, processor_instance, executor, executor_args={}, chunksize=500000, maxchunks=None):
    '''
    A convenience wrapper to submit jobs for a file set, which is a
    dictionary of dataset: [file list] entries.  Supports only uproot
    reading, via the LazyDataFrame class.  For more customized processing,
    e.g. to read other objects from the files and pass them into data frames,
    one can write a similar function in their user code.

    Parameters
    ----------
        fileset:
            dictionary {dataset: [file, file], }
        treename:
            name of tree inside each root file
        processor_instance:
            an instance of a class deriving from ProcessorABC
        executor:
            any of `iterative_executor`, `futures_executor`, etc.

            In general, a function that takes 3 arguments: items, function accumulator
            and performs some action equivalent to:
            for item in items: accumulator += function(item)
        executor_args:
            extra arguments to pass to executor
            currently supported:
                workers: number of parallel processes for futures
                pre_workers: number of parallel threads for calculating chunking
                savemetrics: save some detailed metrics for xrootd processing
                flatten: flatten all branches returned by the dataframe (no jagged structure)
        chunksize:
            number of entries to process at a time in the data frame
        maxchunks:
            maximum number of chunks to process per dataset
    '''
    if not isinstance(fileset, Mapping):
        raise ValueError("Expected fileset to be a mapping dataset: list(files)")
    if not isinstance(processor_instance, ProcessorABC):
        raise ValueError("Expected processor_instance to derive from ProcessorABC")

    executor_args.setdefault('workers', 1)
    executor_args.setdefault('pre_workers', 4 * executor_args['workers'])
    executor_args.setdefault('savemetrics', False)

    items = []
    for dataset, filelist in tqdm(fileset.items(), desc='Preprocessing'):
        if maxchunks is not None:
            chunks = _get_chunking_lazy_noempty(tuple(filelist), treename, chunksize)
        else:
            chunks = _get_chunking(tuple(filelist), treename, chunksize, executor_args['pre_workers'])
        for ichunk, chunk in enumerate(chunks):
            if (maxchunks is not None) and (ichunk > maxchunks):
                break
            items.append((dataset, chunk[0], treename, chunk[1], chunk[2], processor_instance))

    out = processor_instance.accumulator.identity()
    wrapped_out = dict_accumulator({'out': out, 'metrics': dict_accumulator()})
    executor(items, _work_function_nanoaod, wrapped_out, **executor_args)
    processor_instance.postprocess(out)
    if executor_args['savemetrics']:
        return out, wrapped_out['metrics']
    return out


import time
def _work_function_nanoaod(item, flatten=False, savemetrics=False, mmap=False, **_):
    dataset, fn, treename, chunksize, index, processor_instance = item
    if mmap:
        localsource = {}
    else:
        opts = dict(uproot.FileSource.defaults)
        opts.update({'parallel': None})

        def localsource(path):
            return uproot.FileSource(path, **opts)

    file = uproot.open(fn, localsource=localsource)

    tree = file[treename]
    df = LazyDataFrame(tree, chunksize, index, flatten=flatten)
    for name in file['Runs'].keys():
        name = name.decode('utf-8')
        if index==0:
            df[name] = file['Runs'][name].array()
        else:
            df[name] = 0 * file['Runs'][name].array()
    df['dataset'] = dataset
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