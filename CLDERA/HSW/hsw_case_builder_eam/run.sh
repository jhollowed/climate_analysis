#!/bin/bash

# --------- for spinup
TOT_RUN_LENGTH_YEARS=10
TOT_RUN_LENGTH=$(expr $TOT_RUN_LENGTH_YEARS \* 360)
./create_FIDEAL_case_ne16L72.sh 1 $TOT_RUN_LENGTH "_10year_spinup" /global/homes/j/jhollo/repos/climate_analysis/CLDERA/HSW/hsw_case_builder_eam/configs/user_nl_eam_10yearSpinup

# ---------- for scaling tests
# run with 192, 384, 768, 1536, 3072, 6144
#PECOUNT=$1
#./create_FIDEAL_case_ne16L72.sh 1 30 "_pecount_${PECOUNT}" /global/homes/j/jhollo/repos/climate_analysis/CLDERA/HSW/hsw_case_builder_eam/configs/user_nl_eam_10yearSpinup_newTeq $PECOUNT

#./create_FIDEAL_case_ne16L72.sh 1 30 "_30day" /global/homes/j/jhollo/repos/climate_analysis/CLDERA/HSW/hsw_case_builder_eam/configs/user_nl_eam_daily_animation

#./create_FIDEAL_case_ne16L72.sh 1 10 "_10day_testTeq" /global/homes/j/jhollo/repos/climate_analysis/CLDERA/HSW/hsw_case_builder_eam/configs/user_nl_eam_check_Teq
