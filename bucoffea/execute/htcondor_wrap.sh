#!/bin/bash
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

${ARGS[@]:1}