#!/bin/bash
set -e
if [ "$BUCOFFEADEBUG" = true ]; then
    set -x
fi
echo "Starting: $(date)"
echo "Running on: $(hostname)"
echo "uname -a: $(uname -a)"

source /cvmfs/sft.cern.ch/lcg/views/LCG_95apython3/x86_64-centos7-gcc8-opt/setup.sh
echo "Using python at: $(which python)"
python --version

ARGS=("$@")

if [ ! -z "${X509_USER_PROXY}" ]; then
    echo "Initiating VOMS proxy."
    voms-proxy-info -all -file ${X509_USER_PROXY}
fi

if [ ! -z "${VIRTUAL_ENV}" ]; then
    echo "Found VIRTUAL_ENV variable."
    source ${VIRTUAL_ENV}/bin/activate
else
    tar xf *tgz
    rm -rvf *tgz
    ENVNAME="bucoffeaenv"
    python -m venv ${ENVNAME}
    source ${ENVNAME}/bin/activate
    python -m pip install -e bucoffea --no-cache-dir
    export PYTHONPATH="${PWD}/${ENVNAME}/lib/python3.6/site-packages":${PYTHONPATH}
fi

# Copy files to local disk before running
if [ "$BUCOFFEAPREFETCH" = true ]; then
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

echo "Directory content---"
ls -lah .
echo "===================="


echo "Setup done: $(date)"
time buexec ${ARGS[@]}
echo "Run done: $(date)"

echo "Cleaning up."
rm -vf *.root
rm -vf ${FLIST}
echo "End: $(date)"

