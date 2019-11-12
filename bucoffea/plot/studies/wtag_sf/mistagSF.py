from array import array
from bucoffea.plot.stack_plot import *
from klepto.archives import dir_archive
import ROOT

outfilename = "../../../data/sf/ak8/wtag_mistag_SF.root"
outfile = ROOT.TFile.Open(outfilename,'recreate')
newbin = hist.Bin('jetpt',r'AK8 jet $p_{T}$ (GeV)', [200,300,400,800])

inpath = "../../input/merged"
region = 'cr_2m_1ak8_inclusive_v'

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


def get_mistag_rate(hist, dataset):
    h_all = htmp[dataset].integrate('dataset').integrate('wppass')
    h_pass= htmp[dataset].integrate('dataset').integrate('wppass',slice(0.5,1.5))
    sumw_all , sumw2_all  = h_all .values(sumw2=True)[('cr_2m_1ak8_inclusive_v',)]
    sumw_pass, sumw2_pass = h_pass.values(sumw2=True)[('cr_2m_1ak8_inclusive_v',)]
    sumw_ratio, sumw2_ratio = divide_sumw2(sumw_pass, sumw2_pass, sumw_all, sumw2_all)
    return (sumw_ratio, sumw2_ratio)

for year in [2017,2018]:
    for wp in ['loose','loosemd','tight','tightmd']:
        mc = re.compile(f'(VDY.*HT.*|QCD.*|W.*HT.*|ST_|TTJets-MLM_|Diboson_|GJets.*HT.*|ZJetsToNuNu.*){year}')
        data = re.compile(f'MET_{year}')
        distribution = f'ak8_pass{wp}_pt0'
        acc.load(distribution)
        htmp = acc[distribution]

        htmp = merge_extensions(htmp, acc, reweight_pu=True)
        scale_xs_lumi(htmp)
        htmp = merge_datasets(htmp)

        #make stack_plot for all and pass
        acc[distribution+'_all'] = htmp.integrate('wppass')
        acc[distribution+'_pass'] = htmp.integrate('wppass', slice(0.5,1.5))
        make_plot(acc, region=region, distribution=distribution+'_all', year=year, data=data, mc=mc, outdir='./output/stack_plots', output_format='png')
        make_plot(acc, region=region, distribution=distribution+'_pass', year=year, data=data, mc=mc, outdir='./output/stack_plots', output_format='png', ylim=(10e-3,10e3))

        if newbin:
            htmp = htmp.rebin(htmp.axis('jetpt'),newbin)

        edges = htmp.axis('jetpt').edges()
        centers = htmp.axis('jetpt').centers()
        halfwidth = [centers[i]-edges[i] for i in range(len(centers))]

        # extract mistag rate for data and mc
        mistagrate_data,sumw2_data = get_mistag_rate(htmp, data)
        mistagerr_data = np.sqrt(sumw2_data)
        mistagrate_mc,sumw2_mc = get_mistag_rate(htmp, mc)
        mistagerr_mc = np.sqrt(sumw2_mc)
        sf, sumw2_sf = divide_sumw2(mistagrate_data, sumw2_data, mistagrate_mc, sumw2_mc)
        error_sf = np.sqrt(sumw2_sf)

        # make plots for the mistag rates
        fig = plt.figure()
        plt.errorbar(centers, mistagrate_data, xerr=halfwidth, yerr=mistagerr_data)
        plt.xlabel('ak8 jet pt'); plt.ylabel('mistag rate')
        fig.savefig(plot_dir+f'mistagrate_data_{wp}_{year}.png')
        plt.close()
        fig = plt.figure()
        plt.errorbar(centers, mistagrate_mc, xerr=halfwidth, yerr=mistagerr_mc)
        plt.xlabel('ak8 jet pt'); plt.ylabel('mistag rate')
        fig.savefig(plot_dir+f'mistagrate_mc_{wp}_{year}.png')
        plt.close()
        fig = plt.figure()
        plt.errorbar(centers, sf, xerr=halfwidth, yerr=error_sf)
        plt.xlabel('ak8 jet pt'); plt.ylabel('mistag rate data/MC SF'); plt.ylim(0,1.5)
        fig.savefig(plot_dir+f'mistagSF_{wp}_{year}.png')
        plt.close()
        
        # save the SFs as root histograms
        histname = f'Wmistag_{year}_{wp}_ak8_pt'
        h_sf = ROOT.TH1F(histname, histname, len(edges)-1, array('d',edges))
        for ibin in range(len(sf)):
            h_sf.SetBinContent(ibin+1, sf[ibin])
            h_sf.SetBinError(ibin+1, error_sf[ibin])
        h_sf.Write()

outfile.Close()
        


