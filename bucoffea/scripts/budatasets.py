#!/usr/bin/env python

from bucoffea.execute.dataset_definitions import load_lists, short_name
from tabulate import tabulate
import argparse
import re
def parse_args():
    parser = argparse.ArgumentParser(prog='Inspect datasets lists')
    parser.add_argument('-r','--re', type=str,  default='.*',help='Regular expression to filter datasets based on short name')
    parser.add_argument('--relong', type=str,  default='.*',help='Regular expression to filter datasets based on long name')
    parser.add_argument('-g','--grep',type=str, default="", help='Grep-like snippet to filter datasets based on short name')
    parser.add_argument('--greplong',type=str, default="", help='Grep-like snippet to filter datasets based on long name')

    args = parser.parse_args()

    return args

def main():
    args = parse_args()

    datasets = load_lists()
    table = []
    for long in datasets:
        short = short_name(long)
        if not re.match(args.re, short):
            continue
        if not re.match(args.relong, long):
            continue
        if not (args.grep in short):
            continue
        if not (args.greplong in long):
            continue
        table.append([short, long])
    print(tabulate(table))
if __name__ == "__main__":
    main()