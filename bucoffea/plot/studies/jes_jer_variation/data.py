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
            'dataset1' : {'regex' : 'ZJetsToNuNu.*2017', 'region' : 'sr'},
            'dataset2' : {'regex' : 'WJetsToLNu.*2017', 'region' : 'sr'}
        },
        'ewk' : {
            'dataset1' : {'regex' : 'EWKZ2Jets_ZToNuNu.*2017', 'region' : 'sr'},
            'dataset2' : {'regex' : 'EWKW2Jets_WToLNu.*2017', 'region' : 'sr'}
        },
    },
    'znunu_over_wlnu18' : {
        'qcd' : {
            'dataset1' : {'regex' : 'ZJetsToNuNu.*2018', 'region' : 'sr'},
            'dataset2' : {'regex' : 'WJetsToLNu.*2018', 'region' : 'sr'}
        },
        'ewk' : {
            'dataset1' : {'regex' : 'EWKZ2Jets_ZToNuNu.*2018', 'region' : 'sr'},
            'dataset2' : {'regex' : 'EWKW2Jets_WToLNu.*2018', 'region' : 'sr'}
        },
    },
    'znunu_over_zmumu17' : {
        'qcd': {
            'dataset1' : {'regex' : 'ZJetsToNuNu.*2017', 'region' : 'sr'},
            'dataset2' : {'regex' : 'DYJetsToLL.*2017', 'region' : 'cr_2m'},
        },
        'ewk': {
            'dataset1' : {'regex' : 'EWKZ2Jets_ZToNuNu.*2017', 'region' : 'sr'},
            'dataset2' : {'regex' : 'EWKZ2Jets_ZToLL.*2017', 'region' : 'cr_2m'},
        },
    },
    'znunu_over_zmumu18' : {
        'qcd': {
            'dataset1' : {'regex' : 'ZJetsToNuNu.*2018', 'region' : 'sr'},
            'dataset2' : {'regex' : 'DYJetsToLL.*2018', 'region' : 'cr_2m'},
        },
        'ewk': {
            'dataset1' : {'regex' : 'EWKZ2Jets_ZToNuNu.*2018', 'region' : 'sr'},
            'dataset2' : {'regex' : 'EWKZ2Jets_ZToLL.*2018', 'region' : 'cr_2m'},
        },
    },
    'znunu_over_zee17' : {
        'qcd': {
            'dataset1' : {'regex' : 'ZJetsToNuNu.*2017', 'region' : 'sr'},
            'dataset2' : {'regex' : 'DYJetsToLL.*2017', 'region' : 'cr_2e'},
        },
        'ewk': {
            'dataset1' : {'regex' : 'EWKZ2Jets_ZToNuNu.*2017', 'region' : 'sr'},
            'dataset2' : {'regex' : 'EWKZ2Jets_ZToLL.*2017', 'region' : 'cr_2e'},
        },
    },
    'gjets_over_znunu17' : {
        'qcd': {
            'dataset1' : {'regex' : 'GJets_DR-0p4.*2017', 'region' : 'cr_g'},
            'dataset2' : {'regex' : 'ZJetsToNuNu.*2017', 'region' : 'sr'},
        },
        'ewk': {
            'dataset1' : {'regex' : 'GJets_SM_5f_EWK.*2017', 'region' : 'cr_g'},
            'dataset2' : {'regex' : 'EWKZ2Jets_ZToNuNu.*2017', 'region' : 'sr'},
        }
    },
    'gjets_over_znunu18' : {
        'qcd': {
            'dataset1' : {'regex' : 'GJets_DR-0p4.*2018', 'region' : 'cr_g'},
            'dataset2' : {'regex' : 'ZJetsToNuNu.*2018', 'region' : 'sr'},
        },
        'ewk': {
            'dataset1' : {'regex' : 'GJets_SM_5f_EWK.*2018', 'region' : 'cr_g'},
            'dataset2' : {'regex' : 'EWKZ2Jets_ZToNuNu.*2018', 'region' : 'sr'},
        }
    },
    'znunu_over_zee18' : {
        'qcd': {
            'dataset1' : {'regex' : 'ZJetsToNuNu.*2018', 'region' : 'sr'},
            'dataset2' : {'regex' : 'DYJetsToLL.*2018', 'region' : 'cr_2e'},
        },
        'ewk': {
            'dataset1' : {'regex' : 'EWKZ2Jets_ZToNuNu.*2018', 'region' : 'sr'},
            'dataset2' : {'regex' : 'EWKZ2Jets_ZToLL.*2018', 'region' : 'cr_2e'},
        },
    },
    'wlnu_over_wenu17' : {
        'qcd': {
            'dataset1' : {'regex' : 'WJetsToLNu.*2017', 'region' : 'sr'},
            'dataset2' : {'regex' : 'WJetsToLNu.*2017', 'region' : 'cr_1e'},
        },
        'ewk': {
            'dataset1' : {'regex' : 'EWKW2Jets_WToLNu.*2017', 'region' : 'sr'},
            'dataset2' : {'regex' : 'EWKW2Jets_WToLNu.*2017', 'region' : 'cr_1e'},
        },
    },
    'wlnu_over_wenu18' : {
        'qcd': {
            'dataset1' : {'regex' : 'WJetsToLNu.*2018', 'region' : 'sr'},
            'dataset2' : {'regex' : 'WJetsToLNu.*2018', 'region' : 'cr_1e'},
        },
        'ewk': {
            'dataset1' : {'regex' : 'EWKW2Jets_WToLNu.*2018', 'region' : 'sr'},
            'dataset2' : {'regex' : 'EWKW2Jets_WToLNu.*2018', 'region' : 'cr_1e'},
        },
    },
    'wlnu_over_wmunu17' : {
        'qcd': {
            'dataset1' : {'regex' : 'WJetsToLNu.*2017', 'region' : 'sr'},
            'dataset2' : {'regex' : 'WJetsToLNu.*2017', 'region' : 'cr_1m'},
        },
        'ewk': {
            'dataset1' : {'regex' : 'EWKW2Jets_WToLNu.*2017', 'region' : 'sr'},
            'dataset2' : {'regex' : 'EWKW2Jets_WToLNu.*2017', 'region' : 'cr_1m'},
        },
    },
    'wlnu_over_wmunu18' : {
        'qcd': {
            'dataset1' : {'regex' : 'WJetsToLNu.*2018', 'region' : 'sr'},
            'dataset2' : {'regex' : 'WJetsToLNu.*2018', 'region' : 'cr_1m'},
        },
        'ewk': {
            'dataset1' : {'regex' : 'EWKW2Jets_WToLNu.*2018', 'region' : 'sr'},
            'dataset2' : {'regex' : 'EWKW2Jets_WToLNu.*2018', 'region' : 'cr_1m'},
        },
    },
    'wlnu_over_gjets17' : {
        'qcd': {
            'dataset1' : {'regex' : 'WJetsToLNu.*2017', 'region' : 'sr'},
            'dataset2' : {'regex' : 'GJets_DR-0p4.*2017', 'region' : 'cr_g'},
        },
        'ewk': {
            'dataset1' : {'regex' : 'EWKW2Jets_WToLNu.*2017', 'region' : 'sr'},
            'dataset2' : {'regex' : 'GJets_SM_5f_EWK.*2017', 'region' : 'cr_g'},
        }
    },
    'wlnu_over_gjets18' : {
        'qcd': {
            'dataset1' : {'regex' : 'WJetsToLNu.*2018', 'region' : 'sr'},
            'dataset2' : {'regex' : 'GJets_DR-0p4.*2018', 'region' : 'cr_g'},
        },
        'ewk': {
            'dataset1' : {'regex' : 'EWKW2Jets_WToLNu.*2018', 'region' : 'sr'},
            'dataset2' : {'regex' : 'GJets_SM_5f_EWK.*2018', 'region' : 'cr_g'},
        }
    },
}

# Dict mappping dataset names to 
# corresponding regexps
dataset_regex = {
    ### W processes 
    'wlnu17'  : {
        'qcd' : {'title': r'QCD $W\rightarrow \ell \nu$', 'regex': 'WJetsToLNu.*2017', 'region': 'sr'},
        'ewk' : {'title': r'EWK $W\rightarrow \ell \nu$', 'regex': 'EWKW2Jets_WToLNu.*2017', 'region': 'sr'}
    },
    'wlnu18'  : {
        'qcd' : {'title': r'QCD $W\rightarrow \ell \nu$', 'regex': 'WJetsToLNu.*2018', 'region': 'sr'},
        'ewk' : {'title': r'EWK $W\rightarrow \ell \nu$', 'regex': 'EWKW2Jets_WToLNu.*2018', 'region': 'sr'}
    },
    'wenu17'  : {
        'qcd' : {'title': r'QCD $W\rightarrow e \nu$', 'regex': 'WJetsToLNu.*2017', 'region': 'cr_1e'},
        'ewk' : {'title': r'EWK $W\rightarrow e \nu$', 'regex': 'EWKW2Jets_WToLNu.*2017', 'region': 'cr_1e'}
    },
    'wenu18'  : {
        'qcd' : {'title': r'QCD $W\rightarrow e \nu$', 'regex': 'WJetsToLNu.*2018', 'region': 'cr_1e'},
        'ewk' : {'title': r'EWK $W\rightarrow e \nu$', 'regex': 'EWKW2Jets_WToLNu.*2018', 'region': 'cr_1e'}
    },
    'wmunu17'  : {
        'qcd' : {'title': r'QCD $W\rightarrow \mu \nu$', 'regex': 'WJetsToLNu.*2017', 'region': 'cr_1m'},
        'ewk' : {'title': r'EWK $W\rightarrow \mu \nu$', 'regex': 'EWKW2Jets_WToLNu.*2017', 'region': 'cr_1m'}
    },
    'wmunu18'  : {
        'qcd' : {'title': r'QCD $W\rightarrow \mu \nu$', 'regex': 'WJetsToLNu.*2018', 'region': 'cr_1m'},
        'ewk' : {'title': r'EWK $W\rightarrow \mu \nu$', 'regex': 'EWKW2Jets_WToLNu.*2018', 'region': 'cr_1m'}
    },
    ### Z processes
    'zmumu17'  : {
        'qcd' : {'title': r'QCD $Z\rightarrow \mu \mu$', 'regex': 'DYJetsToLL.*2017', 'region': 'cr_2m'},
        'ewk' : {'title': r'EWK $Z\rightarrow \mu \mu$', 'regex': 'EWKZ2Jets_ZToLL.*2017', 'region': 'cr_2m'}
    },
    'zmumu18'  : {
        'qcd' : {'title': r'QCD $Z\rightarrow \mu \mu$', 'regex': 'DYJetsToLL.*2018', 'region': 'cr_2m'},
        'ewk' : {'title': r'EWK $Z\rightarrow \mu \mu$', 'regex': 'EWKZ2Jets_ZToLL.*2018', 'region': 'cr_2m'}
    },
    'zee17'  : {
        'qcd' : {'title': r'QCD $Z\rightarrow ee$', 'regex': 'DYJetsToLL.*2017', 'region': 'cr_2e'},
        'ewk' : {'title': r'EWK $Z\rightarrow ee$', 'regex': 'EWKZ2Jets_ZToLL.*2017', 'region': 'cr_2e'}
    },
    'zee18'  : {
        'qcd' : {'title': r'QCD $Z\rightarrow ee$', 'regex': 'DYJetsToLL.*2018', 'region': 'cr_2e'},
        'ewk' : {'title': r'EWK $Z\rightarrow ee$', 'regex': 'EWKZ2Jets_ZToLL.*2018', 'region': 'cr_2e'}
    },
    'znunu17'  : {
        'qcd' : {'title': r'QCD $Z\rightarrow \nu \nu$', 'regex': 'ZJetsToNuNu.*2017', 'region': 'sr'},
        'ewk' : {'title': r'EWK $Z\rightarrow \nu \nu$', 'regex': 'EWKZ2Jets_ZToNuNu.*2017', 'region': 'sr'}
    },
    'znunu18'  : {
        'qcd' : {'title': r'QCD $Z\rightarrow \nu \nu$', 'regex': 'ZJetsToNuNu.*2018', 'region': 'sr'},
        'ewk' : {'title': r'EWK $Z\rightarrow \nu \nu$', 'regex': 'EWKZ2Jets_ZToNuNu.*2018', 'region': 'sr'}
    },
    ### Photon processes
    'gjets17' : {
        'qcd' : {'title' : r'QCD $\gamma$ + jets', 'regex' : 'GJets_DR-0p4.*2017', 'region' : 'cr_g'},
        'ewk' : {'title' : r'EWK $\gamma$ + jets', 'regex' : 'GJets_SM_5f_EWK.*2017', 'region' : 'cr_g'}
    },
    'gjets18' : {
        'qcd' : {'title' : r'QCD $\gamma$ + jets', 'regex' : 'GJets_DR-0p4.*2018', 'region' : 'cr_g'},
        'ewk' : {'title' : r'EWK $\gamma$ + jets', 'regex' : 'GJets_SM_5f_EWK.*2018', 'region' : 'cr_g'}
    }
    
}

indices_from_tags = {
        'znunu_over_wlnu17' : {
            'qcd' : r'QCD $Z(\nu\nu) / W(\ell\nu)$ 2017',
            'ewk' : r'EWK $Z(\nu\nu) / W(\ell\nu)$ 2017',
        },
        'znunu_over_wlnu18' : {
            'qcd' : r'QCD $Z(\nu\nu) / W(\ell\nu)$ 2018',
            'ewk' : r'EWK $Z(\nu\nu) / W(\ell\nu)$ 2018',
        },
        'znunu_over_zmumu17' : {
            'qcd' : r'QCD $Z(\nu\nu) / Z(\mu\mu)$ 2017',
            'ewk' : r'EWK $Z(\nu\nu) / Z(\mu\mu)$ 2017',
        },
        'znunu_over_zmumu18' : {
            'qcd' : r'QCD $Z(\nu\nu) / Z(\mu\mu)$ 2018',
            'ewk' : r'EWK $Z(\nu\nu) / Z(\mu\mu)$ 2018',
        },
        'znunu_over_zee17' : {
            'qcd' : r'QCD $Z(\nu\nu) / Z(ee)$ 2017',
            'ewk' : r'EWK $Z(\nu\nu) / Z(ee)$ 2017',
        },
        'znunu_over_zee18' : {
            'qcd' : r'QCD $Z(\nu\nu) / Z(ee)$ 2018',
            'ewk' : r'EWK $Z(\nu\nu) / Z(ee)$ 2018',
        },
        'gjets_over_znunu17' : {
            'qcd' : r'QCD $\gamma$ + jets / $Z(\nu\nu)$ 2017',
            'ewk' : r'EWK $\gamma$ + jets / $Z(\nu\nu)$ 2017',
        },
        'gjets_over_znunu18' : {
            'qcd' : r'QCD $\gamma$ + jets / $Z(\nu\nu)$ 2018',
            'ewk' : r'EWK $\gamma$ + jets / $Z(\nu\nu)$ 2018',
        },
        'wlnu_over_wenu17' : {
            'qcd' : r'QCD $W(\ell\nu) / W(e\nu)$ 2017',
            'ewk' : r'EWK $W(\ell\nu) / W(e\nu)$ 2017',
        },
        'wlnu_over_wenu18' : {
            'qcd' : r'QCD $W(\ell\nu) / W(e\nu)$ 2018',
            'ewk' : r'EWK $W(\ell\nu) / W(e\nu)$ 2018',
        },
        'wlnu_over_wmunu17' : {
            'qcd' : r'QCD $W(\ell\nu) / W(\mu\nu)$ 2017',
            'ewk' : r'EWK $W(\ell\nu) / W(\mu\nu)$ 2017',
        },
        'wlnu_over_wmunu18' : {
            'qcd' : r'QCD $W(\ell\nu) / W(\mu\nu)$ 2018',
            'ewk' : r'EWK $W(\ell\nu) / W(\mu\nu)$ 2018',
        },
        'wlnu_over_gjets17' : {
            'qcd' : r'QCD $W(\ell\nu) / \gamma$ + jets 2017',
            'ewk' : r'EWK $W(\ell\nu) / \gamma$ + jets 2017',
        },
        'wlnu_over_gjets18' : {
            'qcd' : r'QCD $W(\ell\nu) / \gamma$ + jets 2018',
            'ewk' : r'EWK $W(\ell\nu) / \gamma$ + jets 2018',
        }
    }