#!/bin/bash

set -e

pause(){
    while read -r -t 0.001; do :; done # dump the buffer
    read -n1 -rsp $'Press any key to continue or Ctrl+C to exit...\n'
}

# 6/27/23
# script for performing 1-degree HSW runs

MACHINE=pm-cpu
COMPILER=intel
PROJECT=m4014
COMPSET=FIDEAL
GRID=ne30pg2_ne30pg2
RES=ne30

BUILD_FLAG=$1
TOT_RUN_LENGTH=$2
SUFFIX=$3
CUSTOMNL=$4
STOP_N=$TOT_RUN_LENGTH

PECOUNT=1350   # quarter number of cubedsphere elements in ne30
QUEUE=regular
WALLCLOCK=08:00:00

MY_E3SM_ROOT="/global/homes/j/jhollo/E3SM/CLDERA-E3SM_PV"
MODEL="${MY_E3SM_ROOT}/cime/scripts/create_newcase"
CASES="/global/homes/j/jhollo/repos/climate_analysis/CLDERA/HSW/hsw_case_builder_eam/cases"
CASENAME="E3SM_${RES}_L72_${COMPSET}${SUFFIX}"
CASE=${CASES}/${CASENAME}
OUTROOT="/pscratch/sd/j/jhollo/E3SM/E3SMv2_cases/hsw_cases"
RUNDIR="${OUTROOT}/${CASENAME}/run"

CONFIGS="/global/homes/j/jhollo/repos/climate_analysis/CLDERA/HSW/hsw_case_builder_eam/configs"
if [ -z "$CUSTOMNL" ]; then
    NL=${CONFIGS}/user_nl_eam_10yearSpinup
else
    NL=$CUSTOMNL
fi

BOLD="\033[1m"
YELLOW="\033[38;5;11m"
RESET="\033[0m"

# -------------- Define case, Build model ---------------

if [[ ! -d "$CASE"  ||  $BUILD_FLAG != "0" ]]; then

    # ---------- do cleaning if requested
    printf "\n\n========== CLEANING ==========\n"
    if [[ "$BUILD_FLAG" != "0" && -d ${CASE} ]]; then
        echo -e ${CASES}"/"${BOLD}${YELLOW}${CASENAME}${RESET}
        echo -e ${BOLD}${YELLOW}${RUNDIR}${RESET}
        read -p $'\n'"Press Y/y to CLEAN case and output above" -n 1 -r
	    echo
	    if [[ $REPLY =~ ^[Yy]$ ]]; then
			rm -rfv ${CASE}
			rm -rfv ${RUNDIR}
		fi 
    fi

    printf "\n\n========== CREATING CASE ==========\n"
    $MODEL --compset $COMPSET --res $GRID --case $CASE --pecount $PECOUNT \
           --output-root $OUTROOT --machine $MACHINE --project $PROJECT
    
    # ---------- configure case
    cd $CASE
    ./xmlchange DEBUG=FALSE,DOUT_S=FALSE,STOP_OPTION=ndays,STOP_N=$STOP_N
    ./xmlchange JOB_WALLCLOCK_TIME=$WALLCLOCK
    ./xmlchange SAVE_TIMING=TRUE
    ./xmlchange JOB_QUEUE=$QUEUE

    printf "\n\n========== COPYING NAMELISTS,SOURCEMODS ==========\n"
    # ---------- copy namelist settings
    cp --verbose $NL ./user_nl_eam
    
    printf "\n\n========== CASE SETUP ==========\n"
    ./case.setup
    
    # ---------- build, submit
    printf "\n\n========== BUILDING, SUBMITTING JOB ==========\n"
    ./case.build 2>&1 | tee ./log.case.buid
    ./case.submit 2>&1 | tee ./log.case.submit
else
    printf "\n\n========== CASE EXISTS; ABORTING ==========\n"
fi
