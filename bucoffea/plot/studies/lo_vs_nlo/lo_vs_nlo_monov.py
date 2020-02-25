#!/usr/bin/env python
#import sys
from bucoffea.plot.stack_plot import *
from klepto.archives import dir_archive
np.seterr(divide='ignore', invalid='ignore')
import argparse
import multiprocessing
import time
import copy
import sys,os


# set acc as global variable
acc = None

def commandline():
    parser = argparse.ArgumentParser(prog='Plotter.')
    parser.add_argument('inpath', type=str, help='Input folder to use.')
    parser.add_argument('--region', type=str, default='.*', help='Region to plot.')
    parser.add_argument('--distribution', type=str, default='.*', help='Distribution to plot.')
    parser.add_argument('--njobs','-j', type=int, default=6, help='number of jobs running plotter in parallell')
    args = parser.parse_args()
    return args

def _load_and_plot(args):
    '''
    args: a dict containing keys like region, distribution, data, mc, etc. for stack_plot function, also with an inpath argument
    calls the stack_plot(...)
    returns nothing
    '''
    global acc
    # make_plot(acc, region=region,distribution=distribution, year=year, data=data, mc=mc, signal=signal, outdir="./output/LO", output_format="png")
    try:
        inpath = args['inpath']
        region = args['region']
        distribution = args['distribution']
        year = args['year']
        data = args['data']
        mc = args['mc']
        signal = args['signal']
        outdir = args['outdir']
        output_format = args['output_format']
    except KeyError:
        print("invalid input: "+str(args))
        return -1
    try:
        make_plot(acc, region=region,distribution=distribution, year=year, data=data, mc=mc, signal=signal, outdir=outdir, output_format=output_format)
        make_plot(acc, region=region,distribution=distribution, year=year, data=data, mc=mc, signal=signal, outdir=outdir+'/'+region, output_format="pdf")
    except AssertionError:
        print("invalid input: "+str(args))
        raise AssertionError
    return 0

    

def main(sysargs):
    global acc
    inpath = sysargs.inpath
    ##### load the archive as global variable
    acc = dir_archive(
        inpath,
        serialized=True,
        compression=0,
        memsize=1e3,
        )
    acc.load('sumw')
    acc.load('sumw_pileup')
    acc.load('nevents')

    distributions_0=['recoil','met','ak8_pt0','ak8_eta0','ak8_tau210','ak8_mass0','ak8_wvsqcd0','ak8_wvsqcdmd0','dimuon_pt','muon_pt0','muon_eta0','muon_phi0','electron_pt0','electron_eta0','electron_phi0','ak8_wvstqcd0','ak8_wvstqcdmd0','ak8_tvsqcd0','ak8_tvsqcdmd0']
    print("preparing the acc for plotting ... ")
    for distribution in distributions_0:
        if not re.match(sysargs.distribution, distribution):
            continue
        try:
            acc.load(distribution)
            acc[distribution] = merge_extensions(acc[distribution], acc, reweight_pu=not ('nopu' in distribution))
            scale_xs_lumi(acc[distribution])
            acc[distribution] = merge_datasets(acc[distribution])
            acc[distribution].axis('dataset').sorting = 'integral'
        except KeyError:
            print("key error with distribution: "+str(distribution))
            return -2

    
    args_list = []
    for year in [2017,2018]:
        mc_map = {
            'cr_1m_v'      : re.compile(f'(Top.*FXFX|Diboson|QCD_HT|DYJetsToLL_M-50_HT_MLM|WJetsToLNu.*HT|GJets_DR.*HT).*{year}'),
            'cr_1e_v'      : re.compile(f'(Top.*FXFX|Diboson|QCD_HT|DYJetsToLL_M-50_HT_MLM|WJetsToLNu.*HT|GJets_DR.*HT).*{year}'),
            'cr_2m_v'      : re.compile(f'(Top.*FXFX|Diboson|QCD_HT|DYJetsToLL_M-50_HT_MLM).*{year}'),
            'cr_2e_v'      : re.compile(f'(Top.*FXFX|Diboson|QCD_HT|DYJetsToLL_M-50_HT_MLM|GJets_DR.*HT).*{year}'),
            'cr_g_v'       : re.compile(f'(Diboson|QCD_HT|GJets_DR.*HT|WJets.*HT).*{year}'),
            #'cr_g_v'       : re.compile(f'(TTJets.*FXFX|Diboson|ST|QCD_HT|GJets_DR.*HT|WJets.*HT).*{year}'),
            'cr_nobveto_v' : re.compile(f'(Top.*FXFX|Diboson|QCD_HT|DYJetsToLL_M-50_HT_MLM|WJetsToLNu.*HT|GJets_DR.*HT|ZJetsToNuNu).*{year}'),
            'sr_v'         : re.compile(f'(Top.*FXFX|Diboson|QCD_HT|DYJetsToLL_M-50_HT_MLM|WJetsToLNu.*HT|GJets_DR.*HT|ZJetsToNuNu).*{year}'),
        }
        #for raw_region in ['cr_1m_v','cr_2m_v','cr_1e_v','cr_2e_v','cr_g_v']:
        for raw_region in ['cr_1m_v','cr_2m_v','cr_1e_v','cr_2e_v','cr_g_v','sr_v','cr_nobveto_v']:
        #for raw_region in ['sr_v']:
            for wp in ['','_inclusive','_loose','_tight','_loosemd','_tightmd']:
            #for wp in ['']:
                region = raw_region.replace('_v',wp+'_v')
                if not re.match(sysargs.region, region):
                    continue
                signal=None

                if ("_1e_" in region) or ("_2e_" in region) or ("_g_" in region): data = re.compile(f'EGamma_{year}')
                elif 'sr_' in region: data=None; signal=re.compile(f'WH.*{year}')
                else: data = re.compile(f'MET_{year}')

                mc = mc_map[raw_region]

                distributions=['recoil','met','ak8_pt0','ak8_eta0','ak8_tau210','ak8_mass0','ak8_wvsqcd0','ak8_wvsqcdmd0','ak8_wvstqcd0','ak8_wvstqcdmd0','ak8_tvsqcd0','ak8_tvsqcdmd0']
                if '_2m_' in region: distributions.append('dimuon_pt')
                if '_1m_' in region: distributions+=['muon_pt0','muon_eta0','muon_phi0']
                if '_1e_' in region: distributions+=['electron_pt0','electron_eta0','electron_phi0']
                for distribution in distributions:
                    if not re.match(sysargs.distribution, distribution):
                        continue
                    args={}
                    args["inpath"]=inpath
                    args["region"]=region
                    args["distribution"]=distribution
                    args["year"]=year
                    args["data"]=data
                    args["mc"]=mc
                    if distribution=="recoil" and '_g_' in region:
                        pattern = mc.pattern.replace("QCD_HT","QCD_data")
                        args["mc"] = re.compile(pattern)
                    args["signal"]=signal
                    #args["outdir"]="./output/" + inpath.split
                    args["outdir"] = pjoin('./output/',list(filter(lambda x:x,inpath.split('/')))[-1])
                    args["output_format"]="png"
                    args_list.append(args)


    pool = multiprocessing.Pool(processes=sysargs.njobs)
    for i in pool.imap_unordered(_load_and_plot, args_list, chunksize=1):
        if i==len(args_list)-1:
            print('terminate')
            pool.terminate()
            break
    print("done")
    


if __name__ == "__main__":
    args = commandline()
    main(args)
