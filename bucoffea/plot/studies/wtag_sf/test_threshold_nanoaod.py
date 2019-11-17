import ROOT

inputfiles = ["root://cmsxrootd.fnal.gov//store/user/aandreas/nanopost/27Oct19/WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8/WJetsToLNu_HT-800To1200-MLM_2018/191028_001851/0000/tree_1.root"]

inputfile = ROOT.TFile.Open(inputfiles[0])
t = inputfile.Get('Events')

tagger="FatJet_deepTagMD_WvsQCD"
base_selections = "FatJet_pt_nom>250"
base_selections+="&& abs(FatJet_eta)<2.4"
base_selections+="&& (FatJet_jetId&2) ==2"
base_selections+="&& FatJet_msoftdrop_nom>65 && FatJet_msoftdrop_nom<105"
ntotal = t.GetEntries(base_selections)
print(f'total entries of ak8 in this file: {ntotal}')
print('cut value\t mistag rate')

for icut in range(50):
    cut_value = 2.0*icut/100
    selections = base_selections + f'&& FatJet_deepTagMD_WvsQCD>{cut_value}'
    ntagged = t.GetEntries(selections)
    percentage = 100.0*ntagged/ntotal
    print(f'{cut_value}\t {percentage}%')

