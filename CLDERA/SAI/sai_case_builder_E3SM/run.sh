#!/bin/bash

#COMPSET=$1
#SPINUP=$2
#BUILD_FLAG=$3
#RES=$4


# --- climatology spinups
#./create_AQP_AMIP_case_ne30_SAI.sh amip 1 1 ne30
#./create_AQP_AMIP_case_ne30_SAI.sh aqp 1 1 ne16

# --- plume runs
./create_AQP_AMIP_case_ne30_SAI.sh amip 0 1 ne30
#./create_AQP_AMIP_case_ne30_SAI.sh aqp 0 0 ne16
