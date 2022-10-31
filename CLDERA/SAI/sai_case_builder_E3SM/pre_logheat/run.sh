#!/bin/bash

# --- hsw plume runs with builtin SAI
#BUILD_FLAG=$1
#TOT_RUN_LENGTH=$2
#SUFFIX=$3
#CUSTOMNL=$4


# --- all tracers passive
#./create_FIDEAL_ne16L72_SAI.sh 1 180 '_passive' ~/repos/climate_analysis/CLDERA/SAI/sai_case_builder_E3SM/current/configs/user_nl_cam_SEne16_HSW_SAI_passive

# --- all pathways active
./create_FIDEAL_ne16L72_SAI.sh 1 180 '_allActive' ~/repos/climate_analysis/CLDERA/SAI/sai_case_builder_E3SM/current/configs/user_nl_cam_SEne16_HSW_SAI_allActive

# --- all pathways active, injection delayed to day 15
#./create_FIDEAL_ne16L72_SAI.sh 1 180 '_allActive_delay15days' ~/repos/climate_analysis/CLDERA/SAI/sai_case_builder_E3SM/current/configs/user_nl_cam_SEne16_HSW_SAI_allActive_delay15days

# --- all tracers passive, np4
#./create_FIDEAL_ne16np4L72_SAI.sh 1 180 '_passive_np4' ~/repos/climate_analysis/CLDERA/SAI/sai_case_builder_E3SM/current/configs/user_nl_cam_SEne16_HSW_SAI_passive

# --- all tracers active, np4
#./create_FIDEAL_ne16np4L72_SAI.sh 1 180 '_allActive_np4' ~/repos/climate_analysis/CLDERA/SAI/sai_case_builder_E3SM/current/configs/user_nl_cam_SEne16_HSW_SAI_allActive

# --- all pathways active, injection delayed to day 15, np4
#./create_FIDEAL_ne16np4L72_SAI.sh 1 180 '_allActive_np4_delay15days' ~/repos/climate_analysis/CLDERA/SAI/sai_case_builder_E3SM/current/configs/user_nl_cam_SEne16_HSW_SAI_allActive_delay15days

