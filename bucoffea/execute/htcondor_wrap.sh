#!/bin/bash
ARGS=("$@")
echo "Arguments: " ${ARGS[@]}

echo "Initiating VOMS proxy."
export X509_USER_PROXY=$(readlink -e ./x509*)
voms-proxy-info

if [ ! -z "${VIRTUAL_ENV}" ]; then
    echo "Found VIRTUAL_ENV variable."
    source ${VIRTUAL_ENV}/bin/activate
fi

${ARGS[@]:1}