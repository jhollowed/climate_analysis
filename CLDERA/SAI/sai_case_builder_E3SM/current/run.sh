#!/bin/bash

# --- hsw plume runs with builtin SAI
#BUILD_FLAG=$1
#TOT_RUN_LENGTH=$2
#SUFFIX=$3
#CUSTOMNL=$4

# --- for heating ensemble
#./create_FIDEAL_ne16L72_heatingEns_e90st80.sh 1 181 '_heatingEns_e90st80' ~/repos/climate_analysis/CLDERA/SAI/sai_case_builder_E3SM/current/configs/user_nl_cam_SEne16_HSW_SAI_heatingEns
#./create_FIDEAL_ne16L72_heatingEns_e90st80.sh 1 181 '_heatingEns_e90st80_passive' ~/repos/climate_analysis/CLDERA/SAI/sai_case_builder_E3SM/current/configs/user_nl_cam_SEne16_HSW_SAI_heatingEns_passive

# --- for spinup
#./create_FIDEAL_ne16L72_spinup_e90st80.sh 1 720 '_e90st80_spinup' ~/repos/climate_analysis/CLDERA/SAI/sai_case_builder_E3SM/current/configs/user_nl_cam_SEne16_HSW_SAI_spinup
# ---- TEMP TO FIX TRUNCATED SPINUP
#./create_FIDEAL_ne16L72_spinup_e90st80.sh 1 361 '_e90st80_spinup2' ~/repos/climate_analysis/CLDERA/SAI/sai_case_builder_E3SM/current/configs/user_nl_cam_SEne16_HSW_SAI_spinup

# --- for tuning
#./create_FIDEAL_ne16L72_SAI.sh 1 100 '_passive' ~/repos/climate_analysis/CLDERA/SAI/sai_case_builder_E3SM/current/configs/user_nl_cam_SEne16_HSW_SAI_passive
#./create_FIDEAL_ne16L72_SAI.sh 1 100 '_allActive' ~/repos/climate_analysis/CLDERA/SAI/sai_case_builder_E3SM/current/configs/user_nl_cam_SEne16_HSW_SAI_allActive

# --- all tracers passive
#./create_FIDEAL_ne16L72_SAI.sh 1 360 '_passive' ~/repos/climate_analysis/CLDERA/SAI/sai_case_builder_E3SM/current/configs/user_nl_cam_SEne16_HSW_SAI_passive

# --- all pathways active
./create_FIDEAL_ne16L72_SAI.sh 1 360 '_allActive' ~/repos/climate_analysis/CLDERA/SAI/sai_case_builder_E3SM/current/configs/user_nl_cam_SEne16_HSW_SAI_allActive

# --- all pathways active, injection delayed to day 15
#./create_FIDEAL_ne16L72_SAI.sh 1 180 '_allActive_delay15days' ~/repos/climate_analysis/CLDERA/SAI/sai_case_builder_E3SM/current/configs/user_nl_cam_SEne16_HSW_SAI_allActive_delay15days


