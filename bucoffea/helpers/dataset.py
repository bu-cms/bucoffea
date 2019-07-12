import re

def is_lo_z(dataset):
    return bool(re.match('(DY|Z)(\d+)Jet.*(mg|MLM).*', dataset))

def is_lo_w(dataset):
    return bool(re.match('W(\d+)Jet.*(mg|MLM).*', dataset))

def is_data(dataset):
    tags = ['EGamma','MET','SingleElectron','SingleMuon','SinglePhoton']
    if any([dataset.startswith(itag) for itag in tags ]):
        return True
    return False


def extract_year(dataset):
    for x in [6,7,8]:
        if f"201{x}" in dataset:
            return 2010+x
    raise RuntimeError("Could not determine dataset year")
