#/bin/bash
while IFS=" " read -r DATASET REMOTE_PATH REMAINDER
do
    FNAME=$(basename $REMOTE_PATH)
    if [[ ! -f $FNAME ]]; then
        wget $REMOTE_PATH
    fi
    echo $FNAME > files.txt
    for PROCESSOR in monojet vbfhinv; do
        buexec $PROCESSOR --outpath ./output/ --jobs 1 worker --dataset WJetsToLNu_HT-400To600-MLM_2017 --filelist files.txt --chunk 0
    done
done < "testfiles.txt"