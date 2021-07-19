#!/bin/bash

get_xs(){
    INPUT=${1}
    echo "Input dataset: " $INPUT

    # if [[ "$INPUT" == *"JetsToLNu"*"FXFX"* ]]; then
    #     echo "Skipping broken XS dataset: ${INPUT}"
    #     return
    if [[ "$INPUT" == *"NANOAODSIM" ]]; then
        DATASET=$(dasgoclient --query="parent dataset=${INPUT}")
        echo "Input is NANO, so I will run on parent dataset ${DATASET}"
    elif [[ "$INPUT" == *"SIM" ]]; then
        DATASET=${INPUT}
    else
        echo "Skipping non-MC dataset: ${INPUT}"
        return
    fi

    # Source CMSSW
    pushd /cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/cmssw/CMSSW_10_6_21/
    eval `scramv1 runtime -sh`
    popd;
    # if [[ "$DATASET"=*"Autumn18"* ]]; then
        # pushd /cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/cmssw/CMSSW_10_6_21/src
        # eval `scramv1 runtime -sh`
        # popd;
    # elif [[ "$DATASET"=*"Fall17"* ]]; then
        # popd /cvmfs/cms.cern.ch/slc7_amd64_gcc630/cms/cmssw/CMSSW_10_6_21/src;
        # eval `scramv1 runtime -sh`
        # popd;
    # fi
    # Download GenXSecAnalyzer config
    if [ ! -e "./ana.py" ]; then
    curl "https://raw.githubusercontent.com/syuvivida/generator/master/cross_section/runJob/ana.py" -o "ana.py"
    fi

    # Get 10 files for dataset, sort by nevents, take file with largest nevents
    FILES=$(das_client --query="file dataset=${DATASET} | grep file.name, file.nevents" | sort -rnk2 | head -5 |  grep -P -oe "/store.*?root" |  xargs -I {} echo -n root://cms-xrd-global.cern.ch//{}, | sed "s|,$||")
    OUTNAME=$(echo ${DATASET} | sed "s|/|_|g;s|RunII.*||g")


    # Output directory for logs
    mkdir -p "./output"

    # Go!
    cmsRun ana.py maxEvents=100000 inputFiles="$FILES" 2>&1 | tee "output/${OUTNAME}.txt"
    echo 'DATASET=${INPUT}' >> "output/${OUTNAME}.txt"
    echo $INPUT $(grep 'final cross section' "output/${OUTNAME}.txt"  | sed 's|.*= ||;s|+-||;s|  | |') >> xs_UL.txt
    echo $DATASET $(grep 'final cross section' "output/${OUTNAME}.txt"  | sed 's|.*= ||;s|+-||;s|  | |') >> xs_mini.txt
}
run_from_list() {
    LIST=${1}
    touch xs_UL.txt
    while read DATASET; do
        if [[ $DATASET == \#* ]]; then
                echo $DATASET
                continue
        elif [ -z "${DATASET}" ]; then
                echo ""
                continue
        fi
        if [[ $(grep -c "$DATASET" xs_UL.txt) -gt 0 ]]; then
                continue
        fi

        get_xs $DATASET
    done < ${LIST}
}

run_from_list ../datasets_UL_mini_2017.txt
run_from_list ../datasets_UL_mini_2018.txt