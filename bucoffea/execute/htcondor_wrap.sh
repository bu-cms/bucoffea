#!/bin/bash
echo "Initiating VOMS proxy."
export X509_USER_PROXY=$(readlink -e ./x509up*)

echo "Arguments: " $@
if [ ! -z "${VIRTUAL_ENV}" ]; then
    echo "Found VIRTUAL_ENV variable."
    source ${VIRTUAL_ENV}/bin/activate
fi
$@
