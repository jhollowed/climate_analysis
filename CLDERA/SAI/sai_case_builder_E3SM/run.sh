#!/bin/bash

#COMPSET=$1
#SPINUP=$2
#BUILD_FLAG=$3
#RES=$4


# --- climatology spinups
#./create_AQP_AMIP_case_ne30_SAI.sh amip 1 1 ne30
#./create_AQP_AMIP_case_ne30_SAI.sh aqp 1 1 ne16

# --- plume runs
#./create_AQP_AMIP_case_ne30_SAI.sh amip 0 1 ne30
#./create_AQP_AMIP_case_ne30_SAI.sh aqp 0 0 ne16

#BUILD_FLAG=$1
#TOT_RUN_LENGTH=$2
#SUFFIX=$3
#CUSTOMNL=$4
./create_FIDEAL_ne16L72_SAI.sh 1 30
