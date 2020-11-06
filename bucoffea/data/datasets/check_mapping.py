#!/usr/bin/env python

import os
import sys
import re
from pprint import pprint
from bucoffea.plot.util import create_dataset_mapping
from bucoffea.execute.dataset_definitions import short_name

# Script to quickly test the validity of dataset mapping in its current form

def get_list_of_datasets(inputfile):
    '''Read the input file, return the list of all the datasets in there.'''
    with open(inputfile, 'r') as f:
        lines = f.readlines()
        # Cleanup the comments/empty lines 
        dataset_list = []
        for line in lines:
            if not line.startswith('/'):
                continue
            dataset_list.append(line)
    
    return dataset_list

def shortened_dataset_list(list_of_datasets):
    '''Given the original list of datasets, return the list of shortened dataset names.'''
    short_name_list = list( map(short_name, list_of_datasets) )
    return short_name_list

def do_printout(mapping,year):
    '''Print out the mapping in a readable format.'''
    for key, dataset_list in mapping.items():
        # If the merged dataset is not from the year we're currently look at, don't bother
        if str(year) not in key:
            continue
        print(f'Merged dataset name: {key}')
        print('*'*20)
        if len(dataset_list) == 0:
            print('Empty')
            pass
        else:
            for dataset in dataset_list:
                print(dataset)
        print('*'*20)

def main():
    for year in [2017, 2018]:
        print('-'*20)
        print(f'INFO: Looking at year: {year}')
        print('-'*20)
        # Input file containing all datasets
        inputfile = f'datasets_nanoaod_v7_{year}.txt'

        dataset_list = get_list_of_datasets(inputfile)
        dataset_list_short = shortened_dataset_list(dataset_list)

        # Now that we have the shortened dataset names, get the mapping and print it
        mapping = create_dataset_mapping(dataset_list_short)
        
        # Do the printout so we can see which datasets are merged to which
        do_printout(mapping,year)

if __name__ == '__main__':
    main()


