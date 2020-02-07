#!/usr/bin/env python

"""
Script to make EGamma photon SFs usable at high pt.

With the out-of-the-box binning provided by EGamma, bins
at pt > 200 GeV have horrible statistics. Since the values
@ > 200 GeV are statistically compatible with the 100-200 GeV
bin, this script takes the weighted average of the two bins
and assigns the averaged value as a constant to both bins.

The weighted average is calculated using the inverse of the
square of the bin error as its weight.

Usage:

./cap_sfs.py Photon_SF_file.root

--> produces Photon_SF_file_capped.root
"""


import numpy as np
import sys
import ROOT as r
infile = sys.argv[1]

# The output file
outfile = r.TFile(infile.replace(".root","_capped.root"),"RECREATE")

for hname in ["EGamma_SF2D"]:
    # Retrieve the original histogram
    f = r.TFile(infile)
    h = f.Get(hname).Clone()
    h.SetDirectory(0)
    f.Close()

    # Loop over X bins
    for ix in range(1,h.GetNbinsX()+1):

        # Loop over Y bins and collect values to be averaged
        vals = []
        errs = []
        for iy in range(1,h.GetNbinsY()+1):
            if h.GetYaxis().GetBinLowEdge(iy) >=100:
                vals.append(h.GetBinContent(ix, iy))
                errs.append(h.GetBinError(ix, iy))


        # Perform the averaging
        vals = np.array(vals)
        errs = np.array(errs)
        w = 1/errs**2
        new_val = np.sum(w * vals)/np.sum(w)
        new_err = np.sqrt(np.sum((w * errs)**2))/np.sum(w)

        # Loop again, write average values
        for iy in range(1,h.GetNbinsY()+1):
            counter = 0
            if h.GetYaxis().GetBinLowEdge(iy) >=100:
                h.SetBinError(ix, iy, new_err)
                h.SetBinContent(ix, iy, new_val)
                counter +=1

    # Save output histogram
    outfile.cd()
    h.Write(hname)
outfile.Close()