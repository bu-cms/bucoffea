import re

def is_lo_z(dataset):
    return bool(re.match('(DY|Z)(\d+)Jet.*(mg|MLM).*', dataset))

def is_lo_w(dataset):
    return bool(re.match('W(\d+)Jet.*(mg|MLM).*', dataset))

def is_data(dataset):
    tags = ['EGamma','MET','SingleElectron','SingleMuon']
    if any([dataset.startswith(itag) for itag in tags ]):
        return True
    return False