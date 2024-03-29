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

MY_E3SM_ROOT="/global/homes/j/jhollo/E3SM/E3SM"
OUTROOT="/global/cscratch1/sd/jhollo/E3SM/E3SMv2_cases/aoa_cases"
CASES="/global/homes/j/jhollo/E3SM/E3SMv2_cases/aoa/cases"
CONFIGS="/global/homes/j/jhollo/E3SM/E3SMv2_cases/aoa/configs"
MODEL="${MY_E3SM_ROOT}/cime/scripts/create_newcase"

#RES="ne4pg2_ne4pg2"
RES="ne4_ne4"
CASENAME="AQPcase_ne4_L72_AOA"
CASE=${CASES}/${CASENAME}
#NCDATA="cam_vcoords_L72_E3SM.nc"
STOP_N=30
PECOUNT=288
NLEV=72
BUILD_FLAG=1

# -------------- Define case, Build model ---------------

if [[ ! -d "$CASE"  ||  $BUILD_FLAG != "0" ]]; then

    # ---------- do cleaning if requested
    printf "\n\n========== CLEANING ==========\n"
    if [ "$BUILD_FLAG" != "0" ]; then
        rm -rfv ${CASE}
        rm -rfv ${OUTROOT}/${CASENAME}/run
    fi

    printf "\n\n========== CREATING CASE ==========\n"
    $MODEL --compset F-EAM-AQP1 --res ne4_ne4 --case $CASE --pecount $PECOUNT \
           --output-root $OUTROOT --machine $MACHINE --compiler $COMPILER --project $PROJECT
    
    # ---------- configure case
    # nadv_11 = 1 for passive clock tracer
    cd $CASE
    ./xmlchange DEBUG=FALSE,DOUT_S=FALSE,STOP_OPTION=ndays,STOP_N=$STOP_N
    ./xmlchange --append --file env_build.xml --id CAM_CONFIG_OPTS --val "-age_of_air_trcs "
    #./xmlchange --append --file env_build.xml --id CAM_CONFIG_OPTS --val "--nlev=$NLEV "
    ./xmlchange JOB_WALLCLOCK_TIME=01:00:00
    ./xmlchange SAVE_TIMING=TRUE
    
    printf "\n\n========== COPYING NAMELISTS,SOURCEMODS ==========\n"
    # ---------- copy namelist settings
    cp --verbose ${CONFIGS}/user_nl_eam_aoa ./user_nl_eam
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
