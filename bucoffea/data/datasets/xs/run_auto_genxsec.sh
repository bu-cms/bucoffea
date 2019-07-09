#!/bin/bash

get_xs(){
    INPUT=${1}
    echo "Input dataset: " $INPUT

    if [[ "$DATASET"=*"NANOAODSIM" ]]; then 
        DATASET=$(dasgoclient --query="parent dataset=${INPUT}")
        echo "Input is NANO, so I will run on parent dataset ${DATASET}"
    else
        DATASET=${INPUT}
    fi

    # Source CMSSW
    pushd /cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/cmssw/CMSSW_10_6_1/
    # popd /cvmfs/cms.cern.ch/slc7_amd64_gcc493/cms/cmssw/CMSSW_7_6_0/src;
    eval `scramv1 runtime -sh`
    popd;

    # Download GenXSecAnalyzer config
    if [ ! -e "./ana.py" ]; then
    curl "https://raw.githubusercontent.com/syuvivida/generator/master/cross_section/runJob/ana.py" -o "ana.py"
    fi

    # Get 10 files for dataset, sort by nevents, take file with largest nevents
    FILE=$(das_client --query="file dataset=${DATASET} | grep file.name, file.nevents" | sort -k3 | head -1 |  grep -P -oe "/store.*?root")
    OUTNAME=$(echo ${DATASET} | sed "s|/|_|g;s|RunII.*||g")


    # Output directory for logs
    mkdir -p "./output"

    # Go!
    cmsRun ana.py maxEvents=50000 inputFiles="root://cms-xrd-global.cern.ch//$FILE" 2>&1 | tee "output/${OUTNAME}.txt"

    echo $INPUT $(grep  ^Total "output/${OUTNAME}.txt"  | awk '{print $11, $13}') >> xs.txt
}
run_from_list() {
    LIST=${1}
    touch xs.txt
    while read DATASET; do
        if [[ $DATASET = \#* ]]; then
                echo $DATASET
                continue
        if [[ $DATASET = *Run201* ]]; then
                continue
        fi
        elif [ -z "${DATASET}" ]; then
                echo ""
                continue
        fi
        if [[ $(grep -c "$DATASET" xs.txt) -gt 0 ]]; then 
                continue
        fi
        
        get_xs $DATASET
    done < ${LIST}
}

run_from_list ../datasets_2017.txt
run_from_list ../datasets_2018.txt


