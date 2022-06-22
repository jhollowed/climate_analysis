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


# --- hsw plume runs
#BUILD_FLAG=$1
#TOT_RUN_LENGTH=$2
#DO_HEATING=$3
#SUFFIX=$4
#CUSTOMNL=$5

#./create_FIDEAL_ne16L72_SAI.sh 1 90 0
#./create_FIDEAL_ne16L72_SAI.sh 1 90 1 '_90days_stratHeating'

#./create_FIDEAL_ne16L72_SAI.sh 1 1 0 '_90days_checkTimesteps' ~/repos/climate_analysis/CLDERA/SAI/sai_case_builder_E3SM/configs/user_nl_cam_aoa_SEne16_HSW_checkTimeSteps


# --- hsw plume runs with builtin SAI
#BUILD_FLAG=$1
#TOT_RUN_LENGTH=$2
#SUFFIX=$3
#CUSTOMNL=$4

./create_FIDEAL_ne16L72_builtin_SAI.sh 1 30
