#!/bin/bash
echo "Starting: $(date)"
source /cvmfs/sft.cern.ch/lcg/views/LCG_95apython3/x86_64-centos7-gcc8-opt/setup.sh

ARGS=("$@")
echo "Arguments: " ${ARGS[@]}
echo "Initiating VOMS proxy."
export X509_USER_PROXY=${1}
voms-proxy-info -all
voms-proxy-info -all -file ${1}

if [ ! -z "${VIRTUAL_ENV}" ]; then
    echo "Found VIRTUAL_ENV variable."
    source ${VIRTUAL_ENV}/bin/activate
fi
echo "Using python at: $(which python)"

time ${ARGS[@]:1}

echo "Ending: $(date)"