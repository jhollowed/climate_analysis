#!/bin/bash
set -e
RUN=$1

# ----spinup
if [ "${RUN}" == '0' ]; then
./create_hsw94_sai.sh 1 0 1 360 0 1 0
fi
# ---- runs 
if [ "${RUN}" == '1' ]; then
./create_hsw94_sai.sh 1 0 0 60 0 1 0
fi
if [ "${RUN}" == '2' ]; then
./create_hsw94_sai.sh 1 1 0 60 0 1 0
fi
