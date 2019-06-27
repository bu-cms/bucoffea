#!/bin/bash
echo "Initiating VOMS proxy."
echo "Arguments: " $ARGS

export X509_USER_PROXY=$1
ARGS=("$@")

if [ ! -z "${VIRTUAL_ENV}" ]; then
    echo "Found VIRTUAL_ENV variable."
    source ${VIRTUAL_ENV}/bin/activate
fi

${ARGS[@]:1}