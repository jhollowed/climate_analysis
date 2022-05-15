#!/bin/bash

#TOT_RUN_LENGTH_YEARS=30
#TOT_RUN_LENGTH=$(expr $TOT_RUN_LENGTH_YEARS \* 360)
#./create_FIDEAL_case_ne16L72.sh 1 $TOT_RUN_LENGTH "_30year_safetystring"

./create_FIDEAL_case_ne16L72.sh 1 30 "_30day" /global/homes/j/jhollo/repos/climate_analysis/CLDERA/HSW/hsw_case_builder_eam/configs/user_nl_eam_daily_animation
