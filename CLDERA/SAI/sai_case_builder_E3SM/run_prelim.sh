#!/bin/bash


# --- hsw plume runs with builtin SAI
#BUILD_FLAG=$1
#TOT_RUN_LENGTH=$2
#SUFFIX=$3
#CUSTOMNL=$4

./create_FIDEAL_ne16L72_prelim_SAI.sh 1 180 '_180day_newTeq_pthwy1' ~/repos/climate_analysis/CLDERA/SAI/sai_case_builder_E3SM/configs/user_nl_cam_SEne16_HSW_builtin_SAI_newTeq_pthwy1

./create_FIDEAL_ne16L72_prelim_SAI.sh 1 180 '_180day_newTeq_pthwy2' ~/repos/climate_analysis/CLDERA/SAI/sai_case_builder_E3SM/configs/user_nl_cam_SEne16_HSW_builtin_SAI_newTeq_pthwy2

./create_FIDEAL_ne16L72_prelim_SAI.sh 1 180 '_180day_newTeq_pthwy12' ~/repos/climate_analysis/CLDERA/SAI/sai_case_builder_E3SM/configs/user_nl_cam_SEne16_HSW_builtin_SAI_newTeq_pthwy12

./create_FIDEAL_ne16L72_prelim_SAI.sh 1 180 '_180day_newTeq_pthwy123' ~/repos/climate_analysis/CLDERA/SAI/sai_case_builder_E3SM/configs/user_nl_cam_SEne16_HSW_builtin_SAI_newTeq_pthwy123

