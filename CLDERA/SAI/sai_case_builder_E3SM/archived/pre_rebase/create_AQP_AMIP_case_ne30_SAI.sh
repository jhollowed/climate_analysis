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
COMPSET=$1
SPINUP=$2
BUILD_FLAG=$3
RES=$4

MY_E3SM_ROOT="/global/homes/j/jhollo/E3SM/E3SM"
OUTROOT="/global/cscratch1/sd/jhollo/E3SM/E3SMv2_cases/sai_cases"
CASES="/global/homes/j/jhollo/repos/climate_analysis/CLDERA/SAI/sai_case_builder_E3SM/cases"
CONFIGS="/global/homes/j/jhollo/repos/climate_analysis/CLDERA/SAI/sai_case_builder_E3SM/configs"
MODEL="${MY_E3SM_ROOT}/cime/scripts/create_newcase"

if [ "$RES"=="ne30" ]; then   GRID="ne30_ne30"
elif [ "$RES"=="ne16" ]; then GRID="ne16_ne16"; fi


PECOUNT=288
NLEV=72

if [ "$COMPSET" == "aqp" ]; then
    COMPSETNAME='F-EAM-AQP1'
elif [ "$COMPSET" == "amip" ]; then
    COMPSETNAME='F2010-CICE'
fi

if [ "$SPINUP" == "1" ]; then
    # a run without the plume for writing init files every month
    STOP_N=30
    RESUBMIT=11
    NL=${CONFIGS}/user_nl_cam_aoa_SE_spinup
    CASENAME="E3SM_case_${RES}_L72_SAI_${COMPSET}"
    WALLCLOCK=04:00:00
elif [ "$SPINUP" == "0" ]; then
    # 2 month plume run
    #STOP_N=60 
    #WALLCLOCK=06:00:00
    # 1 month plume run
    STOP_N=30 
    WALLCLOCK=03:00:00
    NL=${CONFIGS}/user_nl_cam_aoa_SE
    #CASENAME="E3SM_case_${RES}_L72_SAI_${COMPSET}_juneclimo"
    CASENAME="E3SM_case_${RES}_L72_SAI_${COMPSET}_juneclimo_inclmass"
fi
CASE=${CASES}/${CASENAME}

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
    $MODEL --compset $COMPSETNAME --res $GRID --case $CASE --pecount $PECOUNT \
           --output-root $OUTROOT --machine $MACHINE --compiler $COMPILER --project $PROJECT
    
    # ---------- configure case
    # nadv_11 = 1 for passive clock tracer
    cd $CASE
    ./xmlchange DEBUG=FALSE,DOUT_S=FALSE,STOP_OPTION=ndays,STOP_N=$STOP_N
    ./xmlchange --append --file env_build.xml --id CAM_CONFIG_OPTS --val "-age_of_air_trcs "
    ./xmlchange JOB_WALLCLOCK_TIME=$WALLCLOCK
    ./xmlchange SAVE_TIMING=TRUE
    if [ "$SPINUP" == "1" ]; then
        # turn on resubmits each month if spinup run
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
