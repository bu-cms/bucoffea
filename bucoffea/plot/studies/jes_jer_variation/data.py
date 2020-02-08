#!/usr/bin/env python

####################################
# Hold the dictionaries containing information
# about the datasets used in JES/JER variation
# plotting, in ./plot_jes_jer_var.py
####################################

# Dict mapping tags to dataset pairs
# and corresponding regexps
tag_to_dataset_pairs = {
    'znunu_over_wlnu17' : {
        'qcd' : {
            'dataset1' : {'regex' : 'ZJetsToNuNu.*2017', 'region' : 'sr_vbf'},
            'dataset2' : {'regex' : 'WJetsToLNu.*2017', 'region' : 'sr_vbf'}
        },
        'ewk' : {
            'dataset1' : {'regex' : 'EWKZ2Jets_ZToNuNu.*2017', 'region' : 'sr_vbf'},
            'dataset2' : {'regex' : 'EWKW2Jets_WToLNu.*2017', 'region' : 'sr_vbf'}
        },
    },
    'znunu_over_wlnu18' : {
        'qcd' : {
            'dataset1' : {'regex' : 'ZJetsToNuNu.*2018', 'region' : 'sr_vbf'},
            'dataset2' : {'regex' : 'WJetsToLNu.*2018', 'region' : 'sr_vbf'}
        },
        'ewk' : {
            'dataset1' : {'regex' : 'EWKZ2Jets_ZToNuNu.*2018', 'region' : 'sr_vbf'},
            'dataset2' : {'regex' : 'EWKW2Jets_WToLNu.*2018', 'region' : 'sr_vbf'}
        },
    },
    'znunu_over_zmumu17' : {
        'qcd': {
            'dataset1' : {'regex' : 'ZJetsToNuNu.*2017', 'region' : 'sr_vbf'},
            'dataset2' : {'regex' : 'DYJetsToLL.*2017', 'region' : 'cr_2m_vbf'},
        },
        'ewk': {
            'dataset1' : {'regex' : 'EWKZ2Jets_ZToNuNu.*2017', 'region' : 'sr_vbf'},
            'dataset2' : {'regex' : 'EWKZ2Jets_ZToLL.*2017', 'region' : 'cr_2m_vbf'},
        },
    },
    'znunu_over_zmumu18' : {
        'qcd': {
            'dataset1' : {'regex' : 'ZJetsToNuNu.*2018', 'region' : 'sr_vbf'},
            'dataset2' : {'regex' : 'DYJetsToLL.*2018', 'region' : 'cr_2m_vbf'},
        },
        'ewk': {
            'dataset1' : {'regex' : 'EWKZ2Jets_ZToNuNu.*2018', 'region' : 'sr_vbf'},
            'dataset2' : {'regex' : 'EWKZ2Jets_ZToLL.*2018', 'region' : 'cr_2m_vbf'},
        },
    },
    'znunu_over_zee17' : {
        'qcd': {
            'dataset1' : {'regex' : 'ZJetsToNuNu.*2017', 'region' : 'sr_vbf'},
            'dataset2' : {'regex' : 'DYJetsToLL.*2017', 'region' : 'cr_2e_vbf'},
        },
        'ewk': {
            'dataset1' : {'regex' : 'EWKZ2Jets_ZToNuNu.*2017', 'region' : 'sr_vbf'},
            'dataset2' : {'regex' : 'EWKZ2Jets_ZToLL.*2017', 'region' : 'cr_2e_vbf'},
        },
    },
    'znunu_over_zee18' : {
        'qcd': {
            'dataset1' : {'regex' : 'ZJetsToNuNu.*2018', 'region' : 'sr_vbf'},
            'dataset2' : {'regex' : 'DYJetsToLL.*2018', 'region' : 'cr_2e_vbf'},
        },
        'ewk': {
            'dataset1' : {'regex' : 'EWKZ2Jets_ZToNuNu.*2018', 'region' : 'sr_vbf'},
            'dataset2' : {'regex' : 'EWKZ2Jets_ZToLL.*2018', 'region' : 'cr_2e_vbf'},
        },
    },
    'wlnu_over_wenu17' : {
        'qcd': {
            'dataset1' : {'regex' : 'WJetsToLNu.*2017', 'region' : 'sr_vbf'},
            'dataset2' : {'regex' : 'WJetsToLNu.*2017', 'region' : 'cr_1e_vbf'},
        },
        'ewk': {
            'dataset1' : {'regex' : 'EWKW2Jets_WToLNu.*2017', 'region' : 'sr_vbf'},
            'dataset2' : {'regex' : 'EWKW2Jets_WToLNu.*2017', 'region' : 'cr_1e_vbf'},
        },
    },
    'wlnu_over_wenu18' : {
        'qcd': {
            'dataset1' : {'regex' : 'WJetsToLNu.*2018', 'region' : 'sr_vbf'},
            'dataset2' : {'regex' : 'WJetsToLNu.*2018', 'region' : 'cr_1e_vbf'},
        },
        'ewk': {
            'dataset1' : {'regex' : 'EWKW2Jets_WToLNu.*2018', 'region' : 'sr_vbf'},
            'dataset2' : {'regex' : 'EWKW2Jets_WToLNu.*2018', 'region' : 'cr_1e_vbf'},
        },
    },
    'wlnu_over_wmunu17' : {
        'qcd': {
            'dataset1' : {'regex' : 'WJetsToLNu.*2017', 'region' : 'sr_vbf'},
            'dataset2' : {'regex' : 'WJetsToLNu.*2017', 'region' : 'cr_1m_vbf'},
        },
        'ewk': {
            'dataset1' : {'regex' : 'EWKW2Jets_WToLNu.*2017', 'region' : 'sr_vbf'},
            'dataset2' : {'regex' : 'EWKW2Jets_WToLNu.*2017', 'region' : 'cr_1m_vbf'},
        },
    },
    'wlnu_over_wmunu18' : {
        'qcd': {
            'dataset1' : {'regex' : 'WJetsToLNu.*2018', 'region' : 'sr_vbf'},
            'dataset2' : {'regex' : 'WJetsToLNu.*2018', 'region' : 'cr_1m_vbf'},
        },
        'ewk': {
            'dataset1' : {'regex' : 'EWKW2Jets_WToLNu.*2018', 'region' : 'sr_vbf'},
            'dataset2' : {'regex' : 'EWKW2Jets_WToLNu.*2018', 'region' : 'cr_1m_vbf'},
        },
    },
}

# Dict mappping dataset names to 
# corresponding regexps
dataset_regex = {
    ### W processes 
    'wlnu17'  : {
        'qcd' : {'title': r'QCD $W\rightarrow \ell \nu$', 'regex': 'WJetsToLNu.*2017', 'region': 'sr_vbf'},
        'ewk' : {'title': r'EWK $W\rightarrow \ell \nu$', 'regex': 'EWKW2Jets_WToLNu.*2017', 'region': 'sr_vbf'}
    },
    'wlnu18'  : {
        'qcd' : {'title': r'QCD $W\rightarrow \ell \nu$', 'regex': 'WJetsToLNu.*2018', 'region': 'sr_vbf'},
        'ewk' : {'title': r'EWK $W\rightarrow \ell \nu$', 'regex': 'EWKW2Jets_WToLNu.*2018', 'region': 'sr_vbf'}
    },
    'wenu17'  : {
        'qcd' : {'title': r'QCD $W\rightarrow e \nu$', 'regex': 'WJetsToLNu.*2017', 'region': 'cr_1e_vbf'},
        'ewk' : {'title': r'EWK $W\rightarrow e \nu$', 'regex': 'EWKW2Jets_WToLNu.*2017', 'region': 'cr_1e_vbf'}
    },
    'wenu18'  : {
        'qcd' : {'title': r'QCD $W\rightarrow e \nu$', 'regex': 'WJetsToLNu.*2018', 'region': 'cr_1e_vbf'},
        'ewk' : {'title': r'EWK $W\rightarrow e \nu$', 'regex': 'EWKW2Jets_WToLNu.*2018', 'region': 'cr_1e_vbf'}
    },
    'wmunu17'  : {
        'qcd' : {'title': r'QCD $W\rightarrow \mu \nu$', 'regex': 'WJetsToLNu.*2017', 'region': 'cr_1m_vbf'},
        'ewk' : {'title': r'EWK $W\rightarrow \mu \nu$', 'regex': 'EWKW2Jets_WToLNu.*2017', 'region': 'cr_1m_vbf'}
    },
    'wmunu18'  : {
        'qcd' : {'title': r'QCD $W\rightarrow \mu \nu$', 'regex': 'WJetsToLNu.*2018', 'region': 'cr_1m_vbf'},
        'ewk' : {'title': r'EWK $W\rightarrow \mu \nu$', 'regex': 'EWKW2Jets_WToLNu.*2018', 'region': 'cr_1m_vbf'}
    },
    ### Z processes
    'zmumu17'  : {
        'qcd' : {'title': r'QCD $Z\rightarrow \mu \mu$', 'regex': 'DYJetsToLL.*2017', 'region': 'cr_2m_vbf'},
        'ewk' : {'title': r'EWK $Z\rightarrow \mu \mu$', 'regex': 'EWKZ2Jets_ZToLL.*2017', 'region': 'cr_2m_vbf'}
    },
    'zmumu18'  : {
        'qcd' : {'title': r'QCD $Z\rightarrow \mu \mu$', 'regex': 'DYJetsToLL.*2018', 'region': 'cr_2m_vbf'},
        'ewk' : {'title': r'EWK $Z\rightarrow \mu \mu$', 'regex': 'EWKZ2Jets_ZToLL.*2018', 'region': 'cr_2m_vbf'}
    },
    'zee17'  : {
        'qcd' : {'title': r'QCD $Z\rightarrow ee$', 'regex': 'DYJetsToLL.*2017', 'region': 'cr_2e_vbf'},
        'ewk' : {'title': r'EWK $Z\rightarrow ee$', 'regex': 'EWKZ2Jets_ZToLL.*2017', 'region': 'cr_2e_vbf'}
    },
    'zee18'  : {
        'qcd' : {'title': r'QCD $Z\rightarrow ee$', 'regex': 'DYJetsToLL.*2018', 'region': 'cr_2e_vbf'},
        'ewk' : {'title': r'EWK $Z\rightarrow ee$', 'regex': 'EWKZ2Jets_ZToLL.*2018', 'region': 'cr_2e_vbf'}
    },
    'znunu17'  : {
        'qcd' : {'title': r'QCD $Z\rightarrow \nu \nu$', 'regex': 'ZJetsToNuNu.*2017', 'region': 'sr_vbf'},
        'ewk' : {'title': r'EWK $Z\rightarrow \nu \nu$', 'regex': 'EWKZ2Jets_ZToNuNu.*2017', 'region': 'sr_vbf'}
    },
    'znunu18'  : {
        'qcd' : {'title': r'QCD $Z\rightarrow \nu \nu$', 'regex': 'ZJetsToNuNu.*2018', 'region': 'sr_vbf'},
        'ewk' : {'title': r'EWK $Z\rightarrow \nu \nu$', 'regex': 'EWKZ2Jets_ZToNuNu.*2018', 'region': 'sr_vbf'}
    },
    
    
}

