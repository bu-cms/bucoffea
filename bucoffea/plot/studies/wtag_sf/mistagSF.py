from array import array
from bucoffea.plot.stack_plot import *
from coffea.hist.plot import clopper_pearson_interval
from klepto.archives import dir_archive
import ROOT

#inpath = "../../../input/merged_2020_11_18_03Sep20v7_use_1p0_for_wtag_medium_wp"
#inpath = "../../../input/merged_2020_10_14_update_mistag"
#inpath = "../../../input/merged_2020_11_20_03Sep20v7_hasjmrjms_Wmass70To120"
inpath = "../../../input/merged_2020_11_23_03Sep20v7_Wmass65To120"
#inpath = "../../../input/merged_2021_01_18_03Sep20v7_monojv_splitWkfac_splitmistag"

# construct an output dir name
if inpath[-1]=='/':
    inpath = inpath[:-1]
inpath_tag = inpath.split('/')[-1]
bin_scheme=5   # choose which binning scheme to use, see bin_schemes
massden=True   # True if we apply mass cut on denominator in the efficiency calculation
massnum=True   # True if we apply mass cut on numerator in the efficiency calculation
nlogjet=True   # True if we use NLO GJets sample for the mistag measurement in the photon CR
splitsys=True  # True if splited systemtics sources calculated separately, won't affect nominal values, only generate additional hitograms
realVSF=1.0    # This number can be used to scale the matched real V events up/down before doing realV extraction from data
#outdir = 'output_bin2_gvsothers_nojmr'
outdir = f'output_mistag/{inpath_tag}/output_bin{bin_scheme}'
if massden:
    outdir+="_massden"
if massnum:
    outdir+="_massnum"
if nlogjet:
    outdir+="_nlogjet"
if splitsys:
    outdir+="_splitsys"
outdir+="_realVSF"+( str(realVSF).replace(".","p"))
if not os.path.exists(outdir):
    os.makedirs(outdir)

#outfilename = "../../../data/sf/ak8/wtag_mistag_SF.root"
outfilename = f"{outdir}/wtag_mistag_SF.root"
bin_schems = {
        1 : hist.Bin('jetpt',r'AK8 jet $p_{T}$ (GeV)', [200,300,400,800]),
        2 : hist.Bin('jetpt',r'AK8 jet $p_{T}$ (GeV)', [240,1000]),
        3 : hist.Bin('jetpt',r'AK8 jet $p_{T}$ (GeV)', [240,360,500,1000]),
        5 : hist.Bin('jetpt',r'AK8 jet $p_{T}$ (GeV)', [240,280,340,400,500,600,1000]),
        }
newbin = bin_schems[bin_scheme]
#newbin = hist.Bin('jetpt',r'AK8 jet $p_{T}$ (GeV)', [240,360,500,1000])

if splitsys:
    all_sysvar = ["nominal","sysUp","sysDn","topNormUp","topNormDn","vvNormUp","vvNormDn","vgNormUp","vgNormDn","topVTagUp","topVTagDn","vvVTagUp","vvVTagDn","vgVTagUp","vgVTagDn","topBVetoUp","topBVetoDn","vvBVetoUp","vvBVetoDn","vgBVetoUp","vgBVetoDn"]
else:
    all_sysvar = ["nominal","sysUp","sysDn"]

def divide_sumw2(sumw_a, sumw2_a, sumw_b, sumw2_b): #return (sumw_c, sumw2_c) for c=a/b
    # check that the inputs have the same shape
    assert len(set((sumw_a.shape, sumw_b.shape, sumw2_a.shape, sumw2_b.shape))) == 1
    sumw_c  = np.zeros(sumw_a.shape)
    sumw2_c = np.zeros(sumw2_a.shape)
    for ientry in range(len(sumw2_c)):
        if sumw2_b[ientry]>0:
            sumw_c[ientry] = sumw_a[ientry] / sumw_b[ientry]
            sumw2_c[ientry] = (sumw_a[ientry]**2 * sumw2_b[ientry] + sumw_b[ientry]**2 * sumw2_a[ientry]) / (sumw_b[ientry]**2 * sumw_b[ientry]**2)
    return (sumw_c, sumw2_c)


# takes the coffea hists, calculate the efficiency using ROOT and return a TEfficiency containing the efficiencies calculated
def get_mistag_rate(hist, region_all, region_pass, flag='', isData=False): #flag is for histogram naming only
    sumw_all , sumw2_all  = hist.values(sumw2=True, overflow=True)[(region_all,)]
    sumw_pass, sumw2_pass = hist.values(sumw2=True, overflow=True)[(region_pass,)]
    # construct root th1f
    edges = hist.axis('jetpt').edges()
    th1_all = ROOT.TH1F(f'h_all_{flag}',f'h_all{flag}',len(edges)-1, array('d',edges))
    th1_pass= ROOT.TH1F(f'h_pass_{flag}',f'h_all{flag}',len(edges)-1, array('d',edges))
    for ibin in range(len(edges)-1):
        th1_all.SetBinContent(ibin+1, sumw_all[ibin])
        th1_pass.SetBinContent(ibin+1, sumw_pass[ibin])
        # data events are un-weighted
        if isData: continue
        th1_all.SetBinError(ibin+1, np.sqrt(sumw2_all[ibin]))
        th1_pass.SetBinError(ibin+1, np.sqrt(sumw2_pass[ibin]))
    # use TEfficiency
    teff_mistag_rate = ROOT.TEfficiency(th1_pass, th1_all)
    return teff_mistag_rate

# Convert a TEfficiency to TH1F
def efficiency_to_histogram(teff):
    heff = teff.GetCopyTotalHisto()
    heff.SetNameTitle('th1f_'+teff.GetName(), teff.GetTitle())
    Nbins = heff.GetNbinsX()
    for ibin in range(Nbins+1): #including overflow bins
        heff.SetBinContent(ibin, teff.GetEfficiency(ibin))
        heff.SetBinError(ibin, max(teff.GetEfficiencyErrorLow(ibin), teff.GetEfficiencyErrorUp(ibin)))
    return heff


# Divide two TEfficiencies and return a TH1F
# TH1F doesn't support assymetric error bars, using the larger one
def ratio_of_efficiencies(name, title, numerator, denominator):
    # first convert TEfficiencies into histograms
    h_numerator = efficiency_to_histogram(numerator)
    h_denominator = efficiency_to_histogram(denominator)
    h_ratio = h_numerator.Clone(name)
    h_ratio.SetTitle(title)
    h_ratio.Divide(h_denominator)
    return h_ratio


def main():
    # set to True if want to update the mistag root file
    # otherwise just make the plots
    if True:
        outfile = ROOT.TFile.Open(outfilename,'recreate')
    else:
        outfile = None

    # Prepare the acc
    acc = dir_archive(
        inpath,
        serialized=True,
        compression=0,
        memsize=1e3,
        )
    acc.load('sumw')
    acc.load('sumw_pileup')
    acc.load('nevents')
    distribution = 'ak8_pt0'
    distribution_Vmatched = 'ak8_Vmatched_pt0'
    distribution_mass = 'ak8_mass0'
    acc.load(distribution)
    acc.load(distribution_Vmatched)
    acc.load(distribution_mass)
    
    # merge datasets and scale with lumi xs
    htmp = acc[distribution]
    htmp_Vmatched = acc[distribution_Vmatched]
    htmp_mass = acc[distribution_mass]
    htmp = merge_extensions(htmp, acc, reweight_pu=True)
    scale_xs_lumi(htmp)
    htmp = merge_datasets(htmp)
    acc[distribution]=htmp
    htmp_Vmatched = merge_extensions(htmp_Vmatched, acc, reweight_pu=True)
    scale_xs_lumi(htmp_Vmatched)
    htmp_Vmatched = merge_datasets(htmp_Vmatched)
    acc[distribution_Vmatched]=htmp_Vmatched
    htmp_mass = merge_extensions(htmp_mass, acc, reweight_pu=True)
    scale_xs_lumi(htmp_mass)
    htmp_mass = merge_datasets(htmp_mass)
    acc[distribution_mass]=htmp_mass

    acc[distribution].axis('dataset').sorting = 'integral'
    acc[distribution_Vmatched].axis('dataset').sorting = 'integral'
    acc[distribution_mass].axis('dataset').sorting = 'integral'
    
    #binning stuff
    if newbin:
        htmp = htmp.rebin(htmp.axis('jetpt'),newbin)
        htmp_Vmatched = htmp_Vmatched.rebin(htmp_Vmatched.axis('jetpt'),newbin)
    edges = htmp.axis('jetpt').edges()
    centers = htmp.axis('jetpt').centers()
    halfwidth = [centers[i]-edges[i] for i in range(len(centers))]

    for lepton_flag in ['1m','2m','1e','2e','g']:
    #for lepton_flag in ['g']:
        for year in [2017,2018]:
            mc_map = {
                'cr_1m_v'      : re.compile(f'(Top.*FXFX|Diboson|QCD_HT|DYJetsToLL_M-50_HT_MLM|WJetsToLNu.*HT).*{year}'),
                'cr_1e_v'      : re.compile(f'(Top.*FXFX|Diboson|QCD_HT|DYJetsToLL_M-50_HT_MLM|WJetsToLNu.*HT|GJets_DR.*HT).*{year}'),
                'cr_2m_v'      : re.compile(f'(Top.*FXFX|Diboson|DYJetsToLL_M-50_HT_MLM).*{year}'),
                'cr_2e_v'      : re.compile(f'(Top.*FXFX|Diboson|DYJetsToLL_M-50_HT_MLM).*{year}'),
                'cr_g_v'       : re.compile(f'(Diboson|QCD_HT|GJets_DR.*HT|VQQGamma_FXFX|WJetsToLNu.*HT).*{year}'),
                'cr_nobveto_v' : re.compile(f'(Top.*FXFX|Diboson|QCD_HT|DYJetsToLL_M-50_HT_MLM|WJetsToLNu.*HT|GJets_DR.*HT|ZJetsToNuNu).*{year}'),
                'sr_v'         : re.compile(f'(Top.*FXFX|Diboson|QCD_HT|DYJetsToLL_M-50_HT_MLM|WJetsToLNu.*HT|GJets_DR.*HT|ZJetsToNuNu).*{year}'),
            }
            mc_map_noV = {
                'cr_1m_v'      : re.compile(f'(QCD_HT|DYJetsToLL_M-50_HT_MLM|WJetsToLNu.*HT).*{year}'),
                'cr_1e_v'      : re.compile(f'(QCD_HT|DYJetsToLL_M-50_HT_MLM|WJetsToLNu.*HT|GJets_DR.*HT).*{year}'),
                'cr_2m_v'      : re.compile(f'(QCD_HT|DYJetsToLL_M-50_HT_MLM).*{year}'),
                'cr_2e_v'      : re.compile(f'(QCD_HT|DYJetsToLL_M-50_HT_MLM).*{year}'),
                'cr_g_v'       : re.compile(f'(QCD_HT|GJets_DR.*HT|WJetsToLNu.*HT).*{year}'),
                'cr_nobveto_v' : re.compile(f'(QCD_HT|DYJetsToLL_M-50_HT_MLM|WJetsToLNu.*HT|GJets_DR.*HT|ZJetsToNuNu).*{year}'),
                'sr_v'         : re.compile(f'(QCD_HT|DYJetsToLL_M-50_HT_MLM|WJetsToLNu.*HT|GJets_DR.*HT|ZJetsToNuNu).*{year}'),
            }
            mc_map_realV = {
                'cr_1m_v'      : re.compile(f'(Top.*FXFX|Diboson).*{year}'),
                'cr_1e_v'      : re.compile(f'(Top.*FXFX|Diboson).*{year}'),
                'cr_2m_v'      : re.compile(f'(Top.*FXFX|Diboson).*{year}'),
                'cr_2e_v'      : re.compile(f'(Top.*FXFX|Diboson).*{year}'),
                'cr_g_v'       : re.compile(f'(Diboson|VQQGamma_FXFX).*{year}'),
                'cr_nobveto_v' : re.compile(f'(Top.*FXFX|Diboson).*{year}'),
                'sr_v'         : re.compile(f'(Top.*FXFX|Diboson).*{year}'),
            }
            # use NLO GJets for the measurement if need to
            if nlogjet:
                mc_map['cr_g_v']     = re.compile(f'(Diboson|QCD_HT|GJets_1j|VQQGamma_FXFX|WJetsToLNu.*HT).*{year}')
                mc_map_noV['cr_g_v'] = re.compile(f'(QCD_HT|GJets_1j|WJetsToLNu.*HT).*{year}')
            for wp in ['loose','tight','medium']:
                region_all = f'cr_{lepton_flag}_hasmass_inclusive_v'
                region_all_nomass = f'cr_{lepton_flag}_inclusive_v'
                region_pass= f'cr_{lepton_flag}_nomistag_{wp}_v'
                region_pass_nomass= f'cr_{lepton_flag}_nomistag_nomass_{wp}_v'
                mc_All = mc_map[f'cr_{lepton_flag}_v']
                mc_False = mc_map_noV[f'cr_{lepton_flag}_v']
                mc_Real = mc_map_realV[f'cr_{lepton_flag}_v']
                if lepton_flag in ['1e','2e','g']:
                    data = re.compile(f'EGamma_{year}')
                else:
                    data = re.compile(f'MET_{year}')
                    
                ### DEBUG ###
                # print(acc[distribution][mc_All].integrate("region",region_all).values())
                # print(acc[distribution][mc_All].integrate("region",region_pass).values())
                # print(acc[distribution_Vmatched][mc_All].integrate("region",region_pass).values())
                #############
                #make stack_plot for all and pass
                try:
                    acc["alskjxkjo"]
                    make_plot(acc, region=region_all, distribution=distribution, year=year, data=data, mc=mc_All, outdir=f'{outdir}/stack_plots', output_format='png', ylim=(10e-4,5e3))
                    make_plot(acc, region=region_all_nomass, distribution=distribution, year=year, data=data, mc=mc_All, outdir=f'{outdir}/stack_plots', output_format='png', ylim=(10e-4,5e3))
                    make_plot(acc, region=region_all_nomass, distribution=distribution_mass, year=year, data=data, mc=mc_All, outdir=f'{outdir}/stack_plots', output_format='png', ylim=(10e-4,5e3))
                    make_plot(acc, region=region_all, distribution=distribution, year=year, data=None, mc=mc_Real, outdir=f'{outdir}/stack_plots', output_format='png', ylim=(10e-4,5e3), ratio=False, tag="MCHasV")
                    make_plot(acc, region=region_all, distribution=distribution_Vmatched, year=year, data=None, mc=mc_Real, outdir=f'{outdir}/stack_plots', output_format='png', ylim=(10e-4,5e3), ratio=False, tag="MCHasV")
                    make_plot(acc, region=region_all, distribution=distribution, year=year, data=None, mc=mc_False, outdir=f'{outdir}/stack_plots', output_format='png', ylim=(10e-4,5e3), ratio=False, tag="MCNoV")
                    make_plot(acc, region=region_all, distribution=distribution, year=year, data=data, mc=None, outdir=f'{outdir}/stack_plots', output_format='png', ylim=(10e-4,5e3), ratio=False, tag="data")
                    make_plot(acc, region=region_pass, distribution=distribution, year=year, data=data, mc=mc_All, outdir=f'{outdir}/stack_plots', output_format='png', ylim=(10e-4,5e3))
                    make_plot(acc, region=region_pass_nomass, distribution=distribution, year=year, data=data, mc=mc_All, outdir=f'{outdir}/stack_plots', output_format='png', ylim=(10e-4,5e3))
                    make_plot(acc, region=region_pass_nomass, distribution=distribution_mass, year=year, data=data, mc=mc_All, outdir=f'{outdir}/stack_plots', output_format='png', ylim=(10e-4,5e3))
                    make_plot(acc, region=region_pass, distribution=distribution, year=year, data=None, mc=mc_Real, outdir=f'{outdir}/stack_plots', output_format='png', ylim=(10e-4,5e3), ratio=False, tag="MCHasV")
                    make_plot(acc, region=region_pass, distribution=distribution_Vmatched, year=year, data=None, mc=mc_Real, outdir=f'{outdir}/stack_plots', output_format='png', ylim=(10e-4,5e3), ratio=False, tag="MCHasV")
                    make_plot(acc, region=region_pass, distribution=distribution, year=year, data=None, mc=mc_False, outdir=f'{outdir}/stack_plots', output_format='png', ylim=(10e-4,5e3), ratio=False, tag="MCNoV")
                    make_plot(acc, region=region_pass, distribution=distribution, year=year, data=data, mc=None, outdir=f'{outdir}/stack_plots', output_format='png', ylim=(10e-4,5e3), ratio=False, tag="data")
                except ValueError:
                    print(f"Warning(ValueError): skipping plots for lepton_flag={lepton_flag} year={year} wp={wp} due to negative or zero bins")
                except AssertionError:
                    print(f"Warning(AssertionError): skipping plots for lepton_flag={lepton_flag} year={year} wp={wp} due to negative or zero bins")
                except KeyError:
                    print(f"Warning(KeyError): skipping plots for lepton_flag={lepton_flag} year={year} wp={wp} due to negative or zero bins")
                try:
                    make_plot(acc, region=region_pass_nomass, distribution=distribution, year=year, data=data, mc=mc_All, outdir=f'{outdir}/stack_plots', output_format='png', ylim=(10e-4,5e3))
                    make_plot(acc, region=region_pass_nomass, distribution=distribution_mass, year=year, data=data, mc=mc_All, outdir=f'{outdir}/stack_plots', output_format='png', ylim=(10e-4,5e3))
                except:
                    pass
    
    
        
                # extract mistag rate for data and mc
                selector_region_all,selector_region_pass = region_all, region_pass
                if not massden:
                    selector_region_all = region_all_nomass
                if not massnum:
                    selector_region_pass = region_pass_nomass

                for sysvar in all_sysvar:
                    if sysvar=="nominal":
                        sysvar_tag = ""
                    else:
                        sysvar_tag = "_"+sysvar

                    # background substraction from data: remove real Vs
                    h_data = htmp[data].integrate('dataset')
                    #h_mc_Real  = htmp[mc_Real].integrate('dataset')
                    h_mc_False = htmp[mc_False].integrate('dataset')

                    h_mc_Real  = htmp_Vmatched[mc_Real]
                    # vary within systematics
                    # norm for both 10%, bveto unc for top 6% for diboson 2%, vtag unc for both 10% (approx)
                    if sysvar=="sysUp": h_mc_Real.scale(1.15)
                    if sysvar=="sysDn": h_mc_Real.scale(0.85)
                    if sysvar=="topNormUp": h_mc_Real.scale  ( { "Top_FXFX_2017"      : 1.10, "Top_FXFX_2018"      : 1.10} , axis="dataset" ) 
                    if sysvar=="topNormDn": h_mc_Real.scale  ( { "Top_FXFX_2017"      : 0.90, "Top_FXFX_2018"      : 0.90} , axis="dataset" ) 
                    if sysvar=="vvNormUp":  h_mc_Real.scale  ( { "Diboson_2017"       : 1.10, "Diboson_2018"       : 1.10} , axis="dataset" ) 
                    if sysvar=="vvNormDn":  h_mc_Real.scale  ( { "Diboson_2017"       : 0.90, "Diboson_2018"       : 0.90} , axis="dataset" ) 
                    if sysvar=="vgNormUp":  h_mc_Real.scale  ( { "VQQGamma_FXFX_2017" : 1.10, "VQQGamma_FXFX_2018" : 1.10} , axis="dataset" ) 
                    if sysvar=="vgNormDn":  h_mc_Real.scale  ( { "VQQGamma_FXFX_2017" : 0.90, "VQQGamma_FXFX_2018" : 0.90} , axis="dataset" ) 
                    if sysvar=="topVTagUp": h_mc_Real.scale  ( { "Top_FXFX_2017"      : 1.10, "Top_FXFX_2018"      : 1.10} , axis="dataset" ) 
                    if sysvar=="topVTagDn": h_mc_Real.scale  ( { "Top_FXFX_2017"      : 0.90, "Top_FXFX_2018"      : 0.90} , axis="dataset" ) 
                    if sysvar=="vvVTagUp":  h_mc_Real.scale  ( { "Diboson_2017"       : 1.10, "Diboson_2018"       : 1.10} , axis="dataset" ) 
                    if sysvar=="vvVTagDn":  h_mc_Real.scale  ( { "Diboson_2017"       : 0.90, "Diboson_2018"       : 0.90} , axis="dataset" ) 
                    if sysvar=="vgVTagUp":  h_mc_Real.scale  ( { "VQQGamma_FXFX_2017" : 1.10, "VQQGamma_FXFX_2018" : 1.10} , axis="dataset" ) 
                    if sysvar=="vgVTagDn":  h_mc_Real.scale  ( { "VQQGamma_FXFX_2017" : 0.90, "VQQGamma_FXFX_2018" : 0.90} , axis="dataset" ) 
                    if sysvar=="topBVetoUp": h_mc_Real.scale ( { "Top_FXFX_2017"      : 1.06, "Top_FXFX_2018"      : 1.06} , axis="dataset" ) 
                    if sysvar=="topBVetoDn": h_mc_Real.scale ( { "Top_FXFX_2017"      : 0.94, "Top_FXFX_2018"      : 0.94} , axis="dataset" ) 
                    if sysvar=="vvBVetoUp":  h_mc_Real.scale ( { "Diboson_2017"       : 1.02, "Diboson_2018"       : 1.02} , axis="dataset" ) 
                    if sysvar=="vvBVetoDn":  h_mc_Real.scale ( { "Diboson_2017"       : 0.98, "Diboson_2018"       : 0.98} , axis="dataset" ) 
                    if sysvar=="vgBVetoUp":  h_mc_Real.scale ( { "VQQGamma_FXFX_2017" : 1.02, "VQQGamma_FXFX_2018" : 1.02} , axis="dataset" ) 
                    if sysvar=="vgBVetoDn":  h_mc_Real.scale ( { "VQQGamma_FXFX_2017" : 0.98, "VQQGamma_FXFX_2018" : 0.98} , axis="dataset" ) 
                    h_mc_Real  = h_mc_Real.integrate('dataset')
                    h_mc_Real.scale(-1*realVSF) # just for background substraction
                    h_data.add(h_mc_Real)
    
                    teff_mistag_rate_data = get_mistag_rate(h_data, selector_region_all, selector_region_pass, flag=f'data_{lepton_flag}_{wp}_{year}{sysvar_tag}', isData=True)
                    teff_mistag_rate_data.SetNameTitle(f'mistag_rate_data_{lepton_flag}_{wp}_{year}{sysvar_tag}','mistagging rate')
                    teff_mistag_rate_mc = get_mistag_rate(h_mc_False, selector_region_all, selector_region_pass, flag=f'mc_{lepton_flag}_{wp}_{year}{sysvar_tag}', isData=False)
                    teff_mistag_rate_mc.SetNameTitle(f'mistag_rate_mc_{lepton_flag}_{wp}_{year}{sysvar_tag}','mistagging rate')
    
                    # get the scale factors
                    # note that it's impossible to divide two TEfficiency in ROOT, have to do that manually
                    th1_mistag_SF = ratio_of_efficiencies(f'mistag_SF_{lepton_flag}_{wp}_{year}{sysvar_tag}', 'mistag scale factor', teff_mistag_rate_data, teff_mistag_rate_mc)
                    
                    # save the mistag rate and SF histograms into root file
                    if outfile:
                        teff_mistag_rate_data.Write()
                        teff_mistag_rate_mc.Write()
                        th1_mistag_SF.Write()
    
    # soup togather all CR using a weighted_average between the regions:
    for year in [2017,2018]:
        for wp in ['loose','tight','medium']:
            for sysvar in all_sysvar:
                if sysvar=="nominal":
                    sysvar_tag = ""
                else:
                    sysvar_tag = "_"+sysvar
                teff_mistag_rate_data_1e = outfile.Get(f'mistag_rate_data_1e_{wp}_{year}{sysvar_tag}')
                teff_mistag_rate_data_2e = outfile.Get(f'mistag_rate_data_2e_{wp}_{year}{sysvar_tag}')
                teff_mistag_rate_data_1m = outfile.Get(f'mistag_rate_data_1m_{wp}_{year}{sysvar_tag}')
                teff_mistag_rate_data_2m = outfile.Get(f'mistag_rate_data_2m_{wp}_{year}{sysvar_tag}')
                teff_mistag_rate_data_g = outfile.Get(f'mistag_rate_data_g_{wp}_{year}{sysvar_tag}')
                teff_mistag_rate_mc_1e = outfile.Get(f'mistag_rate_mc_1e_{wp}_{year}{sysvar_tag}')
                teff_mistag_rate_mc_2e = outfile.Get(f'mistag_rate_mc_2e_{wp}_{year}{sysvar_tag}')
                teff_mistag_rate_mc_1m = outfile.Get(f'mistag_rate_mc_1m_{wp}_{year}{sysvar_tag}')
                teff_mistag_rate_mc_2m = outfile.Get(f'mistag_rate_mc_2m_{wp}_{year}{sysvar_tag}')
                teff_mistag_rate_mc_g = outfile.Get(f'mistag_rate_mc_g_{wp}_{year}{sysvar_tag}')
                # souped SF for all W/Z regions
                teff_mistag_rate_data_wz = teff_mistag_rate_data_1e + teff_mistag_rate_data_2e\
                        + teff_mistag_rate_data_1m + teff_mistag_rate_data_2m
                teff_mistag_rate_data_wz.SetNameTitle(f'mistag_rate_data_wz_{wp}_{year}{sysvar_tag}', 'souped mistagging rate for W and Z')
                teff_mistag_rate_mc_wz = teff_mistag_rate_mc_1e + teff_mistag_rate_mc_2e\
                        + teff_mistag_rate_mc_1m + teff_mistag_rate_mc_2m
                teff_mistag_rate_mc_wz.SetNameTitle(f'mistag_rate_mc_wz_{wp}_{year}{sysvar_tag}', 'souped mistagging rate for W and Z')
                th1_mistag_SF_wz = ratio_of_efficiencies(f'mistag_SF_wz_{wp}_{year}{sysvar_tag}', 'souped mistag scale factor for W and Z', teff_mistag_rate_data_wz, teff_mistag_rate_mc_wz)
                # souped SF for all W regions
                teff_mistag_rate_data_w = teff_mistag_rate_data_1e + teff_mistag_rate_data_1m
                teff_mistag_rate_data_w.SetNameTitle(f'mistag_rate_data_w_{wp}_{year}{sysvar_tag}', 'souped mistagging rate for W')
                teff_mistag_rate_mc_w = teff_mistag_rate_mc_1e + teff_mistag_rate_mc_1m
                teff_mistag_rate_mc_w.SetNameTitle(f'mistag_rate_mc_w_{wp}_{year}{sysvar_tag}', 'souped mistagging rate for W')
                th1_mistag_SF_w = ratio_of_efficiencies(f'mistag_SF_w_{wp}_{year}{sysvar_tag}', 'souped mistag scale factor for W', teff_mistag_rate_data_w, teff_mistag_rate_mc_w)
                # souped SF for all Z regions
                teff_mistag_rate_data_z = teff_mistag_rate_data_2e + teff_mistag_rate_data_2m
                teff_mistag_rate_data_z.SetNameTitle(f'mistag_rate_data_z_{wp}_{year}{sysvar_tag}', 'souped mistagging rate for Z')
                teff_mistag_rate_mc_z = teff_mistag_rate_mc_2e + teff_mistag_rate_mc_2m
                teff_mistag_rate_mc_z.SetNameTitle(f'mistag_rate_mc_z_{wp}_{year}{sysvar_tag}', 'souped mistagging rate for Z')
                th1_mistag_SF_z = ratio_of_efficiencies(f'mistag_SF_z_{wp}_{year}{sysvar_tag}', 'souped mistag scale factor for Z', teff_mistag_rate_data_z, teff_mistag_rate_mc_z)
                # souped SF for all regions including photon
                teff_mistag_rate_data_all = teff_mistag_rate_data_1e + teff_mistag_rate_data_2e\
                        + teff_mistag_rate_data_1m + teff_mistag_rate_data_2m + teff_mistag_rate_data_g
                teff_mistag_rate_data_all.SetNameTitle(f'mistag_rate_data_all_{wp}_{year}{sysvar_tag}', 'souped mistagging rate for all')
                teff_mistag_rate_mc_all = teff_mistag_rate_mc_1e + teff_mistag_rate_mc_2e\
                        + teff_mistag_rate_mc_1m + teff_mistag_rate_mc_2m + teff_mistag_rate_mc_g
                teff_mistag_rate_mc_all.SetNameTitle(f'mistag_rate_mc_all_{wp}_{year}{sysvar_tag}', 'souped mistagging rate for all')
                th1_mistag_SF_all = ratio_of_efficiencies(f'mistag_SF_all_{wp}_{year}{sysvar_tag}', 'souped mistag scale factor for all', teff_mistag_rate_data_all, teff_mistag_rate_mc_all)
                if outfile:
                    teff_mistag_rate_data_wz.Write()
                    teff_mistag_rate_mc_wz.Write()
                    th1_mistag_SF_wz.Write()
                    teff_mistag_rate_data_w.Write()
                    teff_mistag_rate_mc_w.Write()
                    th1_mistag_SF_w.Write()
                    teff_mistag_rate_data_z.Write()
                    teff_mistag_rate_mc_z.Write()
                    th1_mistag_SF_z.Write()
                    teff_mistag_rate_data_all.Write()
                    teff_mistag_rate_mc_all.Write()
                    th1_mistag_SF_all.Write()
    
    if outfile:
        outfile.Close()
        

if __name__ == "__main__":
    main()
