#!/bin/bash

# for create
#BUILD_FLAG=$1  # either 0,1,or 2;
               # 0 throws error if case exists,
               # 1 will rebuild (without clean-all) and submit,e.g. for namelist changes
               # 2 will nuke the case and create from scratch
#FIX=$2      # turn on mass fix?
#SPINUP=$3   # is spinup?
#STOP_N=$4

#TAUPHYS=$5  # injection timescale to match phys timestep?
#NSPLIT=$6   # nsplit
#NODIFF=$7   # turn off diffusion?

set -e
RUN=$1
STOPN=$2
NODIFF=$3

# ----spinup
if [ "${RUN}" == '0' ]; then
    ./create_hsw94_sai.sh 1 0 1 360 0 1 0
fi

# ---- runs 
if [ "${RUN}" == '1' ]; then
    ./create_hsw94_sai.sh 1 0 0 60 0 1 0
    #./create_hsw94_sai.sh 1 0 0 $2 0 1 $3
fi
if [ "${RUN}" == '2' ]; then
    ./create_hsw94_sai.sh 1 1 0 $2 0 1 $3
fi
