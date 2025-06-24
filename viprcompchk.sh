#!/bin/bash

# usage: specify the build directory of vipr as first argument and the certificate as second.
# example: bash viprchk_comp.sh build examples/paper_eg3_weak.vipr
VIPR="$1"
CERTIFICATE="$2"

COMPLETE_CERTIFICATE="${CERTIFICATE/.vipr/_complete.vipr}"

if [ ! -f ${CERTIFICATE} ]; then
    echo "Certificate ${CERTIFICATE} does not exist."
    exit 
fi


if [ ! -e ${VIPR} ]; then
    echo "Directory $VIPR does not exist."
    exit
fi


echo ">>> Executing: ${VIPR}/viprcomp ${CERTIFICATE}"
eval "${VIPR}/viprcomp ${CERTIFICATE}"

echo ">>> Executing: ${VIPR}/viprchk ${COMPLETE_CERTIFICATE}"
eval "${VIPR}/viprchk ${COMPLETE_CERTIFICATE}"