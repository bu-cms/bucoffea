import os

import numpy as np
import yaml

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

data = np.loadtxt('xs.txt',dtype=str)

if os.path.exists('xs.yml'):
    with open('xs.yml','r') as f:
        xs_dict = yaml.load(f, Loader=yaml.FullLoader)
else:
    xs_dict = {}

for line in data:
    dataset = short_name(str(line[0]))
    xs = float(line[1])
    if not (dataset in xs_dict):
        xs_dict[dataset] = {}
    else:
        xs_dict[dataset]['gen'] = xs


with open('xs.yml','w') as f:
    yaml.dump(xs_dict, f)
