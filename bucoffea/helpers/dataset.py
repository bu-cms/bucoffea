import re

def is_lo_z(dataset):
    return bool(re.match('(DY|Z)(\d*)Jet.*(mg|MLM|madgraph).*', dataset))

def is_lo_z_ewk(dataset):
    return bool(re.match('EWKZ2Jets_ZTo.', dataset))

def is_lo_w(dataset):
    return bool(re.match('W(\d*)Jet.*(mg|MLM).*', dataset))

def is_lo_w_ewk(dataset):
    return bool(re.match('EWKW(Minus|Plus)2Jets_WToLNu.', dataset))

def is_lo_g(dataset):
    return bool(re.match('GJets.*', dataset))

def is_nlo_z(dataset):
    return bool(re.match('(DY|Z)(\d*)Jet.*FXFX.*', dataset))

def is_nlo_w(dataset):
    return bool(re.match('W(\d*)Jet.*FXFX.*', dataset))

def is_data(dataset):
    tags = ['EGamma','MET','SingleElectron','SingleMuon','SinglePhoton','JetHT']
    if any([dataset.startswith(itag) for itag in tags ]):
        return True
    return False


def extract_year(dataset):
    for x in [6,7,8]:
        if f"201{x}" in dataset:
            return 2010+x
    raise RuntimeError("Could not determine dataset year")
