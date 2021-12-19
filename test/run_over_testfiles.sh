#/bin/bash
set -e

test_processor(){
    PROCESSOR="${1}"
    FILELIST="${2}"

    while IFS=" " read -r DATASET REMOTE_PATH REMAINDER
    do
        FNAME=${DATASET}.root
        if [[ ! -f $FNAME ]]; then
            wget ${REMOTE_PATH}/${FNAME}
        fi
        echo ${FNAME} > files.txt

        buexec ${PROCESSOR} --outpath ./output/ --jobs 1 worker --dataset ${DATASET} --filelist files.txt --chunk 0
    done < "testfiles.txt"
}


test_processor monojet testfiles_eoy_v7.txt
test_processor vbfhinv testfiles_ul_v8.txt