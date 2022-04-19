#!/bin/bash
set -e

MASSFIX=$1

#./create_whs94_sai.sh SE C48 72 1 0 0 1 1 0   # spinup!  ---rr


if [ "${MASSFIX}" == '0' ]; then
    ./create_whs94_sai.sh SE C48 72 1 0 0 1 0 0   # tau, fix off, normal run (2mo) ---
    #./create_whs94_sai.sh SE C48 72 1 0 0 1 0 1   # tau, fix off, no diff ---
    #./create_whs94_sai.sh SE C48 72 1 1 0 1 0 1   # tau phys, fix off, no diff ---
    #./create_whs94_sai.sh SE C48 72 1 1 0 6 0 1   # tau phys, fix off, nsplit=6, no diff ---
fi

if [ "${MASSFIX}" == '1' ]; then
    ./create_whs94_sai.sh SE C48 72 1 0 1 1 0 0   # tau, fix on, normal run (2mo) ---
    #./create_whs94_sai.sh SE C48 72 1 0 1 1 0 1   # tau, fix on, no diff ---
    #./create_whs94_sai.sh SE C48 72 1 1 1 1 0 1   # tau phys, fix on, no diff ---
    #./create_whs94_sai.sh SE C48 72 1 1 1 6 0 1   # tau phys, fix on, nsplit=6, no diff ---
fi
