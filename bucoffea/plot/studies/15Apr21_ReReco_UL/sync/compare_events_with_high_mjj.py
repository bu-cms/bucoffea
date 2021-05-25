#!/usr/bin/env python

import os
import sys
import re
import uproot
import pandas as pd
import numpy as np

from matplotlib import pyplot as plt
from coffea import hist
from bucoffea.plot.util import merge_datasets, merge_extensions, scale_xs_lumi, fig_ratio
from klepto.archives import dir_archive
from pprint import pprint

pjoin = os.path.join

columnts_to_rename = {
    'luminosityBlock' : 'lumi',
    'diCleanJet_M' : 'mjj',
    'Leading_jet_pt' : 'leadak4_pt',
    'Leading_jet_eta' : 'leadak4_eta',
    'Subleading_jet_pt' : 'trailak4_pt',
    'Subleading_jet_eta' : 'trailak4_eta',
}

columns_to_take = list(columnts_to_rename.values()) + ['run', 'event']

def merge_dfs(df_bu, df_ic):
    '''Merge dataframes based on run,event,lumi.'''
    merged_df = df_bu.merge(df_ic, on=['run', 'event', 'lumi'], suffixes=('_bu', '_ic'))
    return merged_df

def do_renaming(df_ic):
    '''Rename IC branches so that they're compatible with BU branch naming.'''
    return df_ic.rename(columns=columnts_to_rename)[columns_to_take]

def get_high_mjj_events(merged_df, thresh=2000.):
    mask = merged_df['mjj_bu'] > thresh
    print(merged_df[mask])

def main():
    tag = '21May21_sync'
    f_ic = uproot.open(f'input/{tag}/Skim_ZJetsToNuNu_HT-800To1200.root')['Events']
    f_bu = uproot.open(f'input/{tag}/tree_ZJetsToNuNu_HT-800To1200_MLM_UL_2017.root')['sr_vbf_no_veto_all']

    df_bu = f_bu.pandas.df()[columns_to_take]
    df_ic = f_ic.pandas.df()

    df_ic = do_renaming(df_ic)

    merged_df = merge_dfs(df_bu, df_ic)

    get_high_mjj_events(merged_df)

if __name__ == '__main__':
    main()