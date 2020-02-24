from array import array
from bucoffea.plot.stack_plot import *
from coffea.hist.plot import clopper_pearson_interval
from klepto.archives import dir_archive
import ROOT

#outfilename = "../../../data/sf/ak8/wtag_mistag_SF.root"
outfilename = "wtag_mistag_SF.root"
# set to True if want to update the mistag root file
# otherwise just make the plots
if True:
    outfile = ROOT.TFile.Open(outfilename,'recreate')
else:
    outfile = None
newbin = hist.Bin('jetpt',r'AK8 jet $p_{T}$ (GeV)', [200,300,400,800])
inpath = "../../input/merged"

plot_dir = 'output/'
if not os.path.exists(plot_dir):
    os.makedirs(plot_dir)

acc = dir_archive(
    inpath,
    serialized=True,
    compression=0,
    memsize=1e3,
    )
acc.load('sumw')
acc.load('sumw_pileup')
acc.load('nevents')

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
    sumw_all , sumw2_all  = hist.values(sumw2=True)[(region_all,)]
    sumw_pass, sumw2_pass = hist.values(sumw2=True)[(region_pass,)]
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


for lepton_flag in ['1m','2m','1e','2e']:
    for year in [2017,2018]:
        for wp in ['loose','loosemd','tightmd','tight']:
            region_all = f'cr_{lepton_flag}_hasmass_inclusive_v'
            region_pass= f'cr_{lepton_flag}_nomistag_{wp}_v'
            mc_Real  = re.compile(f'(ST_|TTJets-MLM_|Diboson_){year}')
            mc_False = re.compile(f'(DY.*HT.*|QCD.*|W.*HT.*|GJets_HT.*|ZJetsToNuNu.*){year}')
            mc_All = re.compile(f'(DY.*HT.*|QCD.*|W.*HT.*|ST_|TTJets-MLM_|Diboson_|GJets.*HT.*|ZJetsToNuNu.*){year}')
            if lepton_flag=='1e' or lepton_flag=='2e':
                data = re.compile(f'EGamma_{year}')
            else:
                data = re.compile(f'MET_{year}')
            distribution = f'ak8_pt0'
            acc.load(distribution)
            acc.load('ak8_Vmatched_pt0')
            htmp = acc[distribution]
            htmp_Vmatched = acc['ak8_Vmatched_pt0']
    
            htmp = merge_extensions(htmp, acc, reweight_pu=True)
            scale_xs_lumi(htmp)
            htmp = merge_datasets(htmp)
            acc[distribution]=htmp

            htmp_Vmatched = merge_extensions(htmp_Vmatched, acc, reweight_pu=True)
            scale_xs_lumi(htmp_Vmatched)
            htmp_Vmatched = merge_datasets(htmp_Vmatched)
            acc['ak8_Vmatched_pt0']=htmp_Vmatched
    
            #make stack_plot for all and pass
            make_plot(acc, region=region_all, distribution=distribution, year=year, data=data, mc=mc_All, outdir='./output/stack_plots', output_format='png')
            make_plot(acc, region=region_pass, distribution=distribution, year=year, data=data, mc=mc_All, outdir='./output/stack_plots', output_format='png', ylim=(10e-3,10e3))
            make_plot(acc, region=region_pass, distribution='ak8_Vmatched_pt0', year=year, data=None, mc=mc_All, outdir='./output/stack_plots', output_format='png', ylim=(10e-4,5e3))
    
            #binning stuff
            if newbin:
                htmp = htmp.rebin(htmp.axis('jetpt'),newbin)
                htmp_Vmatched = htmp_Vmatched.rebin(htmp_Vmatched.axis('jetpt'),newbin)
            edges = htmp.axis('jetpt').edges()
            centers = htmp.axis('jetpt').centers()
            halfwidth = [centers[i]-edges[i] for i in range(len(centers))]
    
            # background substraction from data: remove real Vs
            h_data = htmp[data].integrate('dataset')
            h_mc_Real  = htmp_Vmatched[mc_All].integrate('dataset')
            h_mc_False = htmp[mc_False].integrate('dataset')
            h_mc_Real.scale(-1.) # just for background substraction
            h_data.add(h_mc_Real)

            # extract mistag rate for data and mc
            teff_mistag_rate_data = get_mistag_rate(h_data, region_all, region_pass, flag=f'data_{lepton_flag}_{wp}_{year}', isData=True)
            teff_mistag_rate_data.SetNameTitle(f'mistag_rate_data_{lepton_flag}_{wp}_{year}','mistagging rate')
            teff_mistag_rate_mc = get_mistag_rate(h_mc_False, region_all, region_pass, flag=f'mc_{lepton_flag}_{wp}_{year}', isData=False)
            teff_mistag_rate_mc.SetNameTitle(f'mistag_rate_mc_{lepton_flag}_{wp}_{year}','mistagging rate')

            # get the scale factors
            # note that it's impossible to divide two TEfficiency in ROOT, have to do that manually
            th1_mistag_SF = ratio_of_efficiencies(f'mistag_SF_{lepton_flag}_{wp}_{year}', 'mistag scale factor', teff_mistag_rate_data, teff_mistag_rate_mc)
            
            # save the mistag rate and SF histograms into root file
            if outfile:
                teff_mistag_rate_data.Write()
                teff_mistag_rate_mc.Write()
                th1_mistag_SF.Write()

# soup togather all CRs:
for year in [2017,2018]:
    for wp in ['loose','loosemd','tightmd','tight']:
        teff_mistag_rate_data_1e = outfile.Get(f'mistag_rate_data_1e_{wp}_{year}')
        teff_mistag_rate_data_2e = outfile.Get(f'mistag_rate_data_2e_{wp}_{year}')
        teff_mistag_rate_data_1m = outfile.Get(f'mistag_rate_data_1m_{wp}_{year}')
        teff_mistag_rate_data_2m = outfile.Get(f'mistag_rate_data_2m_{wp}_{year}')
        teff_mistag_rate_mc_1e = outfile.Get(f'mistag_rate_mc_1e_{wp}_{year}')
        teff_mistag_rate_mc_2e = outfile.Get(f'mistag_rate_mc_2e_{wp}_{year}')
        teff_mistag_rate_mc_1m = outfile.Get(f'mistag_rate_mc_1m_{wp}_{year}')
        teff_mistag_rate_mc_2m = outfile.Get(f'mistag_rate_mc_2m_{wp}_{year}')
        teff_mistag_rate_data = teff_mistag_rate_data_1e + teff_mistag_rate_data_2e\
                + teff_mistag_rate_data_1m + teff_mistag_rate_data_2m
        teff_mistag_rate_data.SetNameTitle(f'mistag_rate_data_{wp}_{year}', 'souped mistagging rate')
        teff_mistag_rate_mc = teff_mistag_rate_mc_1e + teff_mistag_rate_mc_2e\
                + teff_mistag_rate_mc_1m + teff_mistag_rate_mc_2m
        teff_mistag_rate_mc.SetNameTitle(f'mistag_rate_mc_{wp}_{year}', 'souped mistagging rate')
        th1_mistag_SF = ratio_of_efficiencies(f'Wmistag_{year}_{wp}_ak8_pt', 'souped mistag scale factor', teff_mistag_rate_data, teff_mistag_rate_mc)
        if outfile:
            teff_mistag_rate_data.Write()
            teff_mistag_rate_mc.Write()
            th1_mistag_SF.Write()

    

if outfile:
    outfile.Close()
        


