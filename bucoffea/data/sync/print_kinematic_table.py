from coffea.util import load
from tabulate import tabulate
import numpy 
from pprint import pprint
def unpack(val):
    try:
        return unpack(val[0])
    except:
        return val
    
import sys
output = load(sys.argv[1])

data=output['kinematics']

keys = [
'event',
'met',
'met_phi',
'recoil',
'recoil_phi',
'ak4pt0',
'ak4eta0',
'leadbtag',
'nLooseMu',
'nTightMu',
'mupt0',
'mueta0',
'nLooseEl',
'nTightEl',
'elpt0',
'eleta0',
'nLooseGam',
'nTightGam',
'gpt0',
'geta0'
]

pprint(data)

for i in range(len(data['event'])):
    table = []
    for k in keys:
        line = [k]
        val = unpack(data[k][i])
        if isinstance(val,numpy.ndarray):
            val = "-"
        if isinstance(val,numpy.int64):
            val = int(val)
        if isinstance(val,numpy.float64):
            val = float(val)
        if isinstance(val,numpy.float32):
            val = float(val)
        line.append(val)
        table.append(line)
    print(tabulate(table,floatfmt='%.3f',tablefmt='plain',headers=['Quantity','BU']))
    print()
