import os
import re
import copy
import numpy as np
import yaml
from pprint import pprint
from bucoffea.execute.dataset_definitions import short_name

"""Script to convert the GenXSecAnalyzer list into YAML

The YAML format is a simple nested dictionary, where each top-level
key is the short name of a data set. For each data set, the cross
section as read from GenXSecAnalyzer will be available by the 'gen'
key. Additional keys can be added, e.g. for higher-order cross sections.
This script will *always* overwrite the 'gen' key, but *never* overwrite
other keys. Therefore, manually entered higher-order cross sections
are persistent.

Example:
    DY_2018:
        gen: 1 # from GenXSecAnalyzer
    DY_2017:
        gen: 0.9
"""

infile = 'xs_UL.txt'
outfile = 'xs_UL.yml'

data = np.loadtxt(infile,dtype=str)

if os.path.exists(outfile):
    with open(outfile,'r') as f:
        xs_dict = yaml.load(f, Loader=yaml.FullLoader)
else:
    xs_dict = {}

for line in data:
    dataset = short_name(str(line[0]))
    xs = float(line[1])
    if not (dataset in xs_dict):
        xs_dict[dataset] = {}
    if not 'gen' in xs_dict[dataset]:
        xs_dict[dataset]['gen'] = xs

# tmp = {}
# for k, v in xs_dict.items():
#     k = re.sub('(_ext\d+|_new_+pmx|_PSweights)','',k)
#     if k not in xs_dict:
#         tmp[k] = copy.deepcopy(v)
# xs_dict.update(tmp)

with open(outfile,'w') as f:
    yaml.dump(xs_dict, f)
