#!/bin/bash

# --- cldera hsw sai
#BUILD_FLAG=$1
#TOT_RUN_LENGTH=$2
#SUFFIX=$3
#CUSTOMNL=$4

# --- all pathways active
#./create_FIDEAL_ne16L72_SAI.sh 1 150
#./create_FIDEAL_ne16L72_SAI.sh 1 900 _180delay
#./create_FIDEAL_ne16L72_SAI.sh 1 1200 _180delay_ens05__cf
#./create_FIDEAL_ne16L72_SAI.sh 1 3600 _mc

# --- for T pertlim investigations
#./create_FIDEAL_ne16L72_SAI.sh 1 300 _pertT
#./create_FIDEAL_ne16L72_SAI.sh 1 1200 _ptens

# --- for T pertlim injection ensemble
#./create_FIDEAL_ne16L72_SAI.sh 1 1200 _ptens

# --- testing ttend outputs...
#./create_FIDEAL_ne16L72_SAI.sh 1 3 _TEST

# --- testing ttend outputs...
#./create_FIDEAL_ne16L72_SAI.sh 1 90 _PERLMUTTER_INJECTION_TEST
./create_FIDEAL_ne16L72_SAI.sh 1 300 _PERLMUTTER_INJECTION_TEST
