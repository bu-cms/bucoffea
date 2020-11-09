#!/bin/bash


expand(){
    # Loops over a file with one data set name stub per line
    # And finds all datasets from DAS matching the name 
    # with the given condition string
    INFILE=${1}
    COND=${2}

    # Ensure that the list ends with a newline
    # To avoid throwing off the 'read' statement below
    sed -i  -e '$a\' ${INFILE}

    # for STUB in $(cat ${INFILE}); do
    while read STUB; do
        # Skip commented or empty lines
        if [[ $STUB = \#* ]]; then
            echo $STUB
        elif [ -z "${STUB}" ]; then
            echo ""
        else
            # Expand the stub into a full dataset name in DAS
            dasgoclient --query="dataset=/${STUB}/${COND}/NANOAOD*"
        fi
    done < ${INFILE}
}

INFILE="dataset_names_mc.txt"
expand "${INFILE}" "RunIISummer16*02Apr2020*" | tee datasets_nanoaod_v7_2016.txt
expand "${INFILE}" "RunIIFall17*02Apr2020*"   | tee datasets_nanoaod_v7_2017.txt
expand "${INFILE}" "RunIIAutumn18*02Apr2020*" | tee datasets_nanoaod_v7_2018.txt

INFILE="dataset_names_data.txt"
expand "${INFILE}" "Run2016*02Apr2020*" | tee -a datasets_nanoaod_v7_2016.txt
expand "${INFILE}" "Run2017*02Apr2020*" | tee -a datasets_nanoaod_v7_2017.txt
expand "${INFILE}" "Run2018*02Apr2020*" | tee -a datasets_nanoaod_v7_2018.txt

# The wildcarded strings in our input lists sometimes match stuff we do not want
# So we do some postprocessing to throw out unwanted datasets. Can be adapted at will.
sed -i '/BGen/d' datasets_nanoaod_v7_201*.txt
sed -i '/DYBBJet/d' datasets_nanoaod_v7_201*.txt
sed -i '/DoubleEMEnriched/d' datasets_nanoaod_v7_201*.txt

sed -i '/WJetsToLNu_.*J_.*2017.*/d' datasets_nanoaod_v7_2017.txt
sed -i '/.*LHEWpT_0-50.*/d' datasets_nanoaod_v7_201*.txt
sed -i '/.*LHEWpT_50-150.*/d' datasets_nanoaod_v7_201*.txt
sed -i '/.*CP5up.*/d' datasets_nanoaod_v7_201*.txt
sed -i '/.*CP5down.*/d' datasets_nanoaod_v7_201*.txt
sed -i '/.*CH3.*/d' datasets_nanoaod_v7_201*.txt


sed -i '/\/WZ_TuneCP5_13TeV-pythia8\/RunIIFall17NanoAODv7-PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1\/NANOAODSIM/d' datasets_nanoaod_v7_2017.txt

# sed '/\(Run2016\|\/G.*Jet\)/!d' -i datasets_nanoaod_v7_2016.txt
sed -i '/.*JetHT.*/d' datasets_nanoaod_v7_2016.txt
sed -i '/.*SingleMuon.*/d' datasets_nanoaod_v7_2016.txt
