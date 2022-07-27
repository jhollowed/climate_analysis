#!/bin/bash

set -e

# Created 2022-07-27 16:05:14

CASEDIR="/global/homes/j/jhollo/repos/climate_analysis/CLDERA/SAI/sai_case_builder_E3SM/cldera_sai_examples/cases/E3SM_ne16_L72_FIDEAL_SAICLDERA_TEST"

/global/homes/j/jhollo/E3SM/CLDERA-E3SM_SAIFork_BKP/cime/scripts/create_newcase --compset FIDEAL --res ne16pg2_ne16pg2 --case "${CASEDIR}" --pecount 384 --output-root /global/cscratch1/sd/jhollo/E3SM/E3SMv2_cases/sai_cases --machine cori-knl --compiler intel --project m4014

cd "${CASEDIR}"

./xmlchange DEBUG=FALSE,DOUT_S=FALSE,STOP_OPTION=ndays,STOP_N=10

./xmlchange --append --file env_build.xml --id CAM_CONFIG_OPTS --val -cldera_sai_trcs -verbose

./xmlchange JOB_WALLCLOCK_TIME=00:30:00

./xmlchange SAVE_TIMING=TRUE

./xmlchange JOB_QUEUE=debug

./case.setup

./case.build --clean

./case.submit

