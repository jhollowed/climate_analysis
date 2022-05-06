#!/bin/bash

set -e

pause(){
    while read -r -t 0.001; do :; done # dump the buffer
    read -n1 -rsp $'Press any key to continue or Ctrl+C to exit...\n'
}

# 2/9/22
# script for performing E3SM runs with enabled AOA clock tracers

MACHINE=cori-knl
COMPILER=intel
PROJECT=m4014
COMPSET=FIDEAL
GRID=ne16_ne16
NLEV=72

BUILD_FLAG=$1

# ----- for debug queueing -----
PECOUNT=$2
QUEUE=debug
WALLCLOCK=00:30:00

MY_E3SM_ROOT="/global/homes/j/jhollo/E3SM/E3SM_FIDEAL"      # ------------  important!
OUTROOT="/global/cscratch1/sd/jhollo/E3SM/E3SMv2_cases/hsw_validate_cases"
CONFIGS="/global/homes/j/jhollo/repos/climate_analysis/CLDERA/HSW/hsw_case_builder/configs"
MODEL="${MY_E3SM_ROOT}/cime/scripts/create_newcase"

CASES="/global/homes/j/jhollo/repos/climate_analysis/CLDERA/HSW/hsw_case_builder/cases"
CASENAME="E3SM_${RES}_L72_${COMPSET}_pecount${PECOUNT}"
CASE=${CASES}/${CASENAME}

STOP_N=30
NL=${CONFIGS}/user_nl_eam


# -------------- Define case, Build model ---------------

if [[ ! -d "$CASE"  ||  $BUILD_FLAG != "0" ]]; then

    # ---------- do cleaning if requested
    printf "\n\n========== CLEANING ==========\n"
    if [ "$BUILD_FLAG" != "0" ]; then
        read -p "Press Y/y to CLEAN case ${CASE}, INCLUDING OUTPUT" -n 1 -r
	    echo    # (optional) move to a new line
	    if [[ $REPLY =~ ^[Yy]$ ]]; then
			rm -rfv ${CASE}
			rm -rfv ${OUTROOT}/${CASENAME}/run
		fi 
    fi

    printf "\n\n========== CREATING CASE ==========\n"
    $MODEL --compset $COMPSET --res $GRID --case $CASE --pecount $PECOUNT \
           --output-root $OUTROOT --machine $MACHINE --compiler $COMPILER --project $PROJECT
    
    # ---------- configure case
    # nadv_11 = 1 for passive clock tracer
    cd $CASE
    ./xmlchange DEBUG=FALSE,DOUT_S=FALSE,STOP_OPTION=ndays,STOP_N=$STOP_N
    ./xmlchange --append --file env_build.xml --id CAM_CONFIG_OPTS --val "-age_of_air_trcs "
    ./xmlchange JOB_WALLCLOCK_TIME=$WALLCLOCK
    ./xmlchange SAVE_TIMING=TRUE
    
    printf "\n\n========== COPYING NAMELISTS,SOURCEMODS ==========\n"
    # ---------- copy namelist settings
    cp --verbose $NL ./user_nl_eam
    # ---------- copy sourcemods
    cp --verbose ${CONFIGS}/SourceMods/* ./SourceMods/src.eam/
    
    printf "\n\n========== CASE SETUP ==========\n"
    ./case.setup
    
    # ---------- build, submit
    printf "\n\n========== BUILDING, SUBMITTING JOB ==========\n"
    ./case.build 2>&1 | tee ./log.case.buid
    ./case.submit 2>&1 | tee ./log.case.submit
else
    printf "\n\n========== CASE EXISTS; ABORTING ==========\n"
fi
