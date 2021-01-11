#!/usr/bin/env python
import itertools
import time
import os
from tqdm import tqdm
from coffea.util import load
from klepto.archives import dir_archive

pjoin = os.path.join
import cachetools.func
import multiprocessing
def _load_keys(fn):
    '''Returns the keys saved in a coffea file'''
    return list(map(str, load(fn).keys()))

def _load_acc(args):
    '''Returns the accumulator saved in a coffea file'''
    fn, key = args
    return load(fn)[key]

def _load_and_sum(args):
    """
    merge item from list of coffea files and dump it to file

    For each file, the saved item corresponding to the
    same key is read out. The sum of the individual
    items for the individual files is dumped.

    :param args: Tuple (key to use, file list, output name)
    :type args: tuple
    :return: 0
    :rtype: int
    """

    # Args is a tuple for easy multiprocessing
    key, files, outname = args

    # Load the individual items
    items = []
    for fn in files:
        try:
            items.append(load(fn)[key])
        except KeyError:
            continue
    
    # Recursive merging
    while len(items) > 1:
        x = items.pop(0)
        y = items.pop(0)
        s = x + y
        items.append(s)
    
    assert(len(items)==1)
    
    # dump the content using klepto
    arc = dir_archive(
                    outname,
                    serialized=True,
                    compression=0,
                    memsize=1e3,
                    )
    arc[key] = items[0]
    arc.dump(key)
    arc.clear()
    return 0

class CoffeaMerger(object):
    '''
    Handles the merging of large numbers of coffea files.
    
    The results are stored using the klepto library.
    '''
    def __init__(self, indir, jobs=1, save_trees=False):
        files = filter(lambda x: x.endswith(".coffea") and not ('cache' in x), os.listdir(indir))
        files = list(map(lambda x: os.path.abspath(pjoin(indir, x)), files))
        self._files = files
        self._keys = set()
        self._save_trees = save_trees

        # Open a multiproc pool for various operations
        self._pool = multiprocessing.Pool(processes=jobs)
        self._init_keys()

    def _init_keys(self):
        '''
        Construct set of all keys.
        '''
        pool_result = self._pool.map_async(
            _load_keys,
            self._files
        )
        for keys in pool_result.get():
            for k in keys:
                if not self._save_trees and k.startswith("tree"):
                    continue
                self._keys.add(str(k))

    def to_klepto_dir(self, outname):
        '''
        Run the merging and save to a klepto dir.
        '''

        # Queue asynchronous jobs for each key
        results = []
        for key in self._keys:
            # Each job loads and merges
            # the items for a specific key
            # from all files
            args = [(key, self._files, outname)]
            result = self._pool.map_async(
                    _load_and_sum,
                    args
                    )
            results.append(result)
        
        # Waiting for results and saving to file
        t = tqdm(total=len(self._keys), desc='Merging inputs')
        while len(results):
            for res in reversed(results):
                if res.ready():
                    results.remove(res)
                    t.update()
            time.sleep(1)
