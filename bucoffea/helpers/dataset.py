import re
def is_lo_znunu(dataset):
    return bool(re.match(r'Z(\d*)Jet.*(mg|MLM|madgraph).*', dataset))

def is_lo_z(dataset):
    return bool(re.match(r'(DY|Z)(\d*)Jet.*(mg|MLM|madgraph).*', dataset))

def is_lo_z_ewk(dataset):
    return bool(re.match(r'EWKZ2Jets_ZTo.', dataset))

def is_lo_w(dataset):
    return bool(re.match(r'W(\d*)Jet.*(mg|MLM).*', dataset))

def is_lo_w_ewk(dataset):
    return bool(re.match(r'EWKW(Minus|Plus)2Jets_WToLNu.', dataset))

def is_lo_g(dataset):
    return bool(re.match(r'GJets.*HT.*', dataset))

def is_lo_g_ewk(dataset):
    return bool(re.match(r'GJets.*EWK.*', dataset))

def is_nlo_g(dataset):
    return bool(re.match(r'G(\d)*Jet.*(amc|NLO).*', dataset))

def is_nlo_g_ewk(dataset):
    return bool(re.match(r'AJJ.*amc.*', dataset))

def is_nlo_z(dataset):
    return bool(re.match(r'(DY|Z)(\d*)Jet.*FXFX.*', dataset))

def is_nlo_w(dataset):
    return bool(re.match(r'W(\d*)Jet.*FXFX.*', dataset))

def has_v_jet(dataset):
    return bool(re.match(r'(WW|WZ|ZZ|TTJets|TTToHadronic|.*WToQQ|.*ZToQQ).*', dataset))

def is_data(dataset):
    tags = ['EGamma','MET','SingleElectron','SingleMuon','SinglePhoton','JetHT']
    if any([dataset.startswith(itag) for itag in tags ]):
        return True
    if re.match('QCD_data_(\d)+',dataset):
        return True
    return False


def extract_year(dataset):
    for x in [6,7,8]:
        if f"201{x}" in dataset:
            return 2010+x
    raise RuntimeError("Could not determine dataset year")
