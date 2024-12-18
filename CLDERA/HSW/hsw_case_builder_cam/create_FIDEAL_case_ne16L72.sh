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
COMPSET=FHS94
GRID=ne16_ne16_mg17
RES=ne16
NLEV=72

BUILD_FLAG=$1
TOT_RUN_LENGTH_YEARS=$2
PREFIX=$3

# ----- for debug queueing -----
PECOUNT=768   # half number of cubedsphere elements in ne16
QUEUE=debug
WALLCLOCK=00:30:00

MY_CESM_ROOT="/global/homes/j/jhollo/CESM"
MODEL="${MY_CESM_ROOT}/cime/scripts/create_newcase"

CONFIGS="/global/homes/j/jhollo/repos/climate_analysis/CLDERA/HSW/hsw_case_builder_cam/configs"
NL=${CONFIGS}/user_nl_eam

CASES="/global/homes/j/jhollo/repos/climate_analysis/CLDERA/HSW/hsw_case_builder_cam/cases"
CASENAME="E3SM_${RES}_L72_${COMPSET}${PREFIX}"
CASE=${CASES}/${CASENAME}

OUTROOT="/global/cscratch1/sd/jhollo/CAM/hsw_validate_cases"
RUNDIR="${OUTROOT}/${CASENAME}/run"

# total run length in days
TOT_RUN_LENGTH=$(expr $TOT_RUN_LENGTH_YEARS \* 360)
# two years per run (sub-30 min for debug queue at ne16 with pecount 768)
STOP_N=720

# check if resubmissions should be done
DO_RESUBS=false
if [ "$TOT_RUN_LENGTH" -gt "$STOP_N" ]; then
    DO_RESUBS=true
fi

# if so, compute number of resubs to effectively round up total length of simulation 
# by one unit of STOP_N
RESUBMIT=0
if $DO_RESUBS; then
    RESUBMIT=$(expr $(expr $TOT_RUN_LENGTH - 1) / $STOP_N)
                           # the divide expr rounds down, so RESUBMIT * STOP_N would give 
                           # TOT_RUN_LENGTH - STOP_N
                           # this is correct, since RESUBMIT should be the total number of desired
                           # submissions, minus one for the initial run
    RESUBMIT=$(expr $RESUBMIT - 1) 
                           # since we are starting from E3SM IC after one run of length STOP_N
fi

# -------------- Define case, Build model ---------------

if [[ ! -d "$CASE"  ||  $BUILD_FLAG != "0" ]]; then

    # ---------- do cleaning if requested
    printf "\n\n========== CLEANING ==========\n"
    if [[ "$BUILD_FLAG" != "0" && -d ${CASE} ]]; then
        read -p "\nPress Y/y to CLEAN case and output \n${CASE}\n${RUNDIR}" -n 1 -r
	    echo
	    if [[ $REPLY =~ ^[Yy]$ ]]; then
			rm -rfv ${CASE}
			rm -rfv ${RUNDIR}
		fi 
    fi

    printf "\n\n========== CREATING CASE ==========\n"
    $MODEL --compset $COMPSET --res $GRID --case $CASE --pecount $PECOUNT \
           --output-root $OUTROOT --machine $MACHINE --compiler $COMPILER --project $PROJECT \
    
    # ---------- configure case
    # nadv_11 = 1 for passive clock tracer
    cd $CASE
    ./xmlchange DEBUG=FALSE,DOUT_S=FALSE,STOP_OPTION=ndays,STOP_N=$STOP_N
    ./xmlchange --append --file env_build.xml --id CAM_CONFIG_OPTS --val "-age_of_air_trcs "
    ./xmlchange JOB_WALLCLOCK_TIME=$WALLCLOCK
    ./xmlchange SAVE_TIMING=TRUE
    ./xmlchange JOB_QUEUE=$QUEUE

    # ---------- configure run restarts
    # setting RESUBMIT>0 automatically sets CONTINUE_RUN=TRUE for all runs
    # subsequent the initial run
    # see: https://www.cesm.ucar.edu/events/tutorials/2018/files/Practical2-shields.pdf
    if $DO_RESUBS; then
        ./xmlchange RESUBMIT=$RESUBMIT
        ./xmlchange REST_OPTION=ndays
        ./xmlchange REST_N=$STOP_N
    fi
    
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
