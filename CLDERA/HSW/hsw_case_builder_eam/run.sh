#!/bin/bash

# --------- for quick test runs
#TOT_RUN_LENGTH=2
#./create_FIDEAL_case_ne16L72.sh 1 $TOT_RUN_LENGTH "_testrun" /global/homes/j/jhollo/repos/climate_analysis/CLDERA/HSW/hsw_case_builder_eam/configs/user_nl_eam_testrun

# --------- for doing autocorrelation test
#TOT_RUN_LENGTH_YEARS=15
#TOT_RUN_LENGTH=$(expr $TOT_RUN_LENGTH_YEARS \* 360)
#./create_FIDEAL_case_ne16L72.sh 1 $TOT_RUN_LENGTH "_autoCorr_15year" /global/homes/j/jhollo/repos/climate_analysis/CLDERA/HSW/hsw_case_builder_eam/configs/user_nl_eam_autoCorr_15year

# --------- for doing HSW-HS modifications comparison runs for HSW paper appendix
#TOT_RUN_LENGTH_YEARS=8
#TOT_RUN_LENGTH=$(expr $TOT_RUN_LENGTH_YEARS \* 360)
#./create_FIDEAL_case_ne16L72.sh 1 $TOT_RUN_LENGTH "_HsHswMods_comparison_HS_8year" /global/homes/j/jhollo/repos/climate_analysis/CLDERA/HSW/hsw_case_builder_eam/configs/user_nl_eam_HsHswMods_comparison_HS
#./create_FIDEAL_case_ne16L72.sh 1 $TOT_RUN_LENGTH "_HsHswMods_comparison_HSWMods_8year" /global/homes/j/jhollo/repos/climate_analysis/CLDERA/HSW/hsw_case_builder_eam/configs/user_nl_eam_HsHswMods_comparison_HSWMods
#./create_FIDEAL_case_ne16L72.sh 1 $TOT_RUN_LENGTH "_HsHswMods_comparison_HSW_8year" /global/homes/j/jhollo/repos/climate_analysis/CLDERA/HSW/hsw_case_builder_eam/configs/user_nl_eam_HsHswMods_comparison_HSW

# --------- for doing HSW 5-year run with daily output for TEM for HSW paper BDC appendix
TOT_RUN_LENGTH_YEARS=5
TOT_RUN_LENGTH=$(expr $TOT_RUN_LENGTH_YEARS \* 360)
./create_FIDEAL_case_ne16L72.sh 1 $TOT_RUN_LENGTH "_HSW_TEM_daily_5year" /global/homes/j/jhollo/repos/climate_analysis/CLDERA/HSW/hsw_case_builder_eam/configs/user_nl_eam_HSW_TEM_daily


# --------- for generating monthly ICs
#TOT_RUN_LENGTH_YEARS=4
#TOT_RUN_LENGTH=$(expr $TOT_RUN_LENGTH_YEARS \* 360)
# --- ne16
#./create_FIDEAL_case_ne16L72.sh 1 $TOT_RUN_LENGTH "_3year_ensICs" /global/homes/j/jhollo/repos/climate_analysis/CLDERA/HSW/hsw_case_builder_eam/configs/user_nl_eam_3yearEnsICS
# --- ne30
#./create_FIDEAL_case_ne30L72.sh 1 $TOT_RUN_LENGTH "_3year_ensICs" /global/homes/j/jhollo/repos/climate_analysis/CLDERA/HSW/hsw_case_builder_eam/configs/user_nl_eam_3yearEnsICS_ne30

# --------- for spinup
#TOT_RUN_LENGTH_YEARS=10
#TOT_RUN_LENGTH=$(expr $TOT_RUN_LENGTH_YEARS \* 360)
# --- ne16
#./create_FIDEAL_case_ne16L72.sh 1 $TOT_RUN_LENGTH "_10year_spinup" /global/homes/j/jhollo/repos/climate_analysis/CLDERA/HSW/hsw_case_builder_eam/configs/user_nl_eam_10yearSpinup
# --- ne30
#./create_FIDEAL_case_ne30L72.sh 1 $TOT_RUN_LENGTH "_10year_spinup" /global/homes/j/jhollo/repos/climate_analysis/CLDERA/HSW/hsw_case_builder_eam/configs/user_nl_eam_10yearSpinup

# ---------- for scaling tests
# run with 192, 384, 768, 1536, 3072, 6144
#PECOUNT=$1
#./create_FIDEAL_case_ne16L72.sh 1 30 "_pecount_${PECOUNT}" /global/homes/j/jhollo/repos/climate_analysis/CLDERA/HSW/hsw_case_builder_eam/configs/user_nl_eam_10yearSpinup_newTeq $PECOUNT

#./create_FIDEAL_case_ne16L72.sh 1 30 "_30day" /global/homes/j/jhollo/repos/climate_analysis/CLDERA/HSW/hsw_case_builder_eam/configs/user_nl_eam_daily_animation

#./create_FIDEAL_case_ne16L72.sh 1 10 "_10day_testTeq" /global/homes/j/jhollo/repos/climate_analysis/CLDERA/HSW/hsw_case_builder_eam/configs/user_nl_eam_check_Teq
