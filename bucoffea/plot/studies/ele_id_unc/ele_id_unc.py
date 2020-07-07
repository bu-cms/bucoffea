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
import uproot


# set acc as global variable
acc = None

def commandline():
    parser = argparse.ArgumentParser(prog='Plotter.')
    parser.add_argument('inpath', type=str, help='Input folder to use.')
    parser.add_argument('--region', type=str, default='.*', help='Region to plot.')
    args = parser.parse_args()
    return args


def main(sysargs):
    global acc
    inpath = sysargs.inpath
    outdir = f'./output/{os.path.basename(inpath)}'
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

    print("preparing the acc for plotting ... ")
    recoil_binning = [ 250., 280., 310., 340., 370., 400., 430., 470., 510., 550., 590., 640., 690., 740., 790., 840., 900., 960., 1020., 1090., 1160., 1250., 1400.]
    recoil_bin = hist.Bin("recoil", "Recoil (GeV)", recoil_binning)
    for distribution in ['recoil','recoil_ele_id_nm','recoil_ele_id_up','recoil_ele_id_dn']:
        try:
            acc.load(distribution)
            acc[distribution] = merge_extensions(acc[distribution], acc, reweight_pu=not ('nopu' in distribution))
            scale_xs_lumi(acc[distribution])
            acc[distribution] = merge_datasets(acc[distribution])
            acc[distribution].axis('dataset').sorting = 'integral'
            # rebin
            acc[distribution] = acc[distribution].rebin("recoil",recoil_bin)
        except KeyError:
            print("key error with distribution: "+str(distribution))
            return -2

    ############Test block###############
    Debug=False
    if Debug:
        region = "cr_1e_j"
        mc = re.compile(f'(Top.*FXFX|Diboson|QCD_HT|DYJetsToLL_M-50_HT_MLM|WJetsToLNu.*HT|GJets_DR.*HT).*2018')
        h_df =           acc["recoil"].integrate("region",region).integrate("dataset",mc)
        h_nm = acc["recoil_ele_id_nm"].integrate("region",region).integrate("dataset",mc)
        h_up = acc["recoil_ele_id_up"].integrate("region",region).integrate("dataset",mc)
        h_dn = acc["recoil_ele_id_dn"].integrate("region",region).integrate("dataset",mc)
        
        fig,axs = plt.subplots(2,1,sharex=True)
        x = h_nm.axis("recoil").centers()
        y_df = h_df.values()[()]
        
        y_nm = h_nm.values()[()]
        line=axs[0].plot(x,y_nm,label="nominal")
        axs[1].plot(x, y_nm/y_df, color=line[0].get_color())
        
        y_up = h_up.values()[()] 
        line=axs[0].plot(x,y_up, label="up")
        axs[1].plot(x, y_up/y_df, color=line[0].get_color())
        
        y_dn = h_dn.values()[()] 
        line=axs[0].plot(x,y_dn, label="dn")
        axs[1].plot(x, y_dn/y_df, color=line[0].get_color())
        
        axs[1].set_xlabel("recoil")
        axs[0].set_xlim(250,1400)
        axs[1].set_ylim(0.9,1.05)
        axs[0].set_yscale('log')
        axs[0].legend()
        
        fig.savefig("test1.png")
        exit()
    #####################################

    outfile = uproot.recreate(outdir+"/ele_id_unc.root")
    
    for year in [2017,2018]:
        mc_map = {
            'cr_1e_v'      : re.compile(f'(Top.*FXFX|Diboson|QCD_HT|DYJetsToLL_M-50_HT_MLM|WJetsToLNu.*HT|GJets_DR.*HT).*{year}'),
            'cr_2e_v'      : re.compile(f'(Top.*FXFX|Diboson|DYJetsToLL_M-50_HT_MLM).*{year}'),
            'cr_1e_j'      : re.compile(f'(Top.*FXFX|Diboson|QCD_HT|DYJetsToLL_M-50_HT_MLM|WJetsToLNu.*HT|GJets_DR.*HT).*{year}'),
            'cr_2e_j'      : re.compile(f'(Top.*FXFX|Diboson|DYJetsToLL_M-50_HT_MLM).*{year}'),
        }
        #for raw_region in ['cr_1m_v','cr_2m_v','cr_1e_v','cr_2e_v','cr_g_v']:
        for raw_region in ['cr_1e_v','cr_2e_v','cr_1e_j','cr_2e_j']:
        #for raw_region in ['sr_v']:
            for wp in ['','_inclusive','_loose','_tight']:
            #for wp in ['']:
                region = raw_region.replace('_v',wp+'_v')
                if len(wp)>0 and '_j' in region:
                    continue
                if not re.match(sysargs.region, region):
                    continue

                ### name for histogram
                if "_v" in region:
                    basename = "monov"+wp.replace("_","")
                else:
                    basename = "monoj"
                if "1e" in region:
                    basename = basename+"_1e_"
                elif "2e" in region:
                    basename = basename+"_2e_"

                mc = mc_map[raw_region]
                hists = {
                        "nm"    :       acc["recoil_ele_id_nm"],
                        "up"    :       acc["recoil_ele_id_up"],
                        "dn"    :       acc["recoil_ele_id_dn"],
                        }
                x = hists["nm"].axis("recoil").centers()
                y_nm = hists["nm"].integrate("dataset",mc).integrate("region",region).values()[()]
                fig,ax = plt.subplots(1,1)
                for tag in ["up","dn"]:
                    y = hists[tag].integrate("dataset",mc).integrate("region",region).values()[()]
                    y = y/y_nm
                    ax.plot(x,y, label="electron ID weight "+tag)
                    ### store as histogram
                    outfile[basename+tag] = np.array(y), np.array(recoil_binning,dtype=float)
                ax.set_xlim(250,2000)
                ax.set_ylim(0.9,1.05)
                ax.legend()
                ax.set_xlabel("recoil (GeV)")
                ax.set_ylabel("weight scale")
                fig.savefig(outdir+"/weight_scale_"+region+".png")


    


if __name__ == "__main__":
    args = commandline()
    main(args)
