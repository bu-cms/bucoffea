#!/bin/bash
set -x
echo "Starting: $(date)"
source /cvmfs/sft.cern.ch/lcg/views/LCG_95apython3/x86_64-centos7-gcc8-opt/setup.sh

ARGS=("$@")
echo "Arguments: " ${ARGS[@]}
echo "Initiating VOMS proxy."
export X509_USER_PROXY=${1}
voms-proxy-info -all -file ${1}

if [ ! -z "${VIRTUAL_ENV}" ]; then
    echo "Found VIRTUAL_ENV variable."
    source ${VIRTUAL_ENV}/bin/activate
fi
echo "Using python at: $(which python)"

# Copy files to local disk before running
PREFETCH=true
if $PREFETCH; then
    echo "Prefetching."
    FLIST=$(readlink -e input*.txt)
    touch tmp.txt
    while read file; do
        LOCAL=$(echo "${file}" | md5sum | awk '{print $1}').root
        xrdcp $file ./$LOCAL
        echo $LOCAL >> tmp.txt
    done < $FLIST
    mv tmp.txt $FLIST
fi

executable=${2}
cp -v ${executable} .

echo "Directory content---"
ls -lah .
echo "===================="


echo "Setup done: $(date)"
time python $(basename ${executable}) ${ARGS[@]:2}
echo "Run done: $(date)"

echo "Cleaning up."
rm -v *.root

echo "End: $(date)"