#!/bin/bash
set -e
if [ "$BUCOFFEADEBUG" = true ]; then
    set -x
fi
echo "Starting: $(date)"
source /cvmfs/sft.cern.ch/lcg/views/LCG_95apython3/x86_64-centos7-gcc8-opt/setup.sh

ARGS=("$@")

if [ ! -z "${X509_USER_PROXY}" ]; then
    echo "Initiating VOMS proxy."
    voms-proxy-info -all -file ${X509_USER_PROXY}
fi

if [ ! -z "${VIRTUAL_ENV}" ]; then
    echo "Found VIRTUAL_ENV variable."
    source ${VIRTUAL_ENV}/bin/activate
fi
echo "Using python at: $(which python)"

# Copy files to local disk before running
if [ ! "$NOPREFETCH" = true ]; then
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

executable=${1}
cp -v ${executable} .

echo "Directory content---"
ls -lah .
echo "===================="


echo "Setup done: $(date)"
time python $(basename ${executable}) ${ARGS[@]:1}
echo "Run done: $(date)"

echo "Cleaning up."
rm -vf *.root
rm -vf ${FLIST}
rm -vf $(basename $executable)
echo "End: $(date)"

