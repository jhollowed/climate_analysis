#!/bin/bash

set -e

pause(){
    while read -r -t 0.001; do :; done # dump the buffer
    read -n1 -rsp $'Press any key to continue or Ctrl+C to exit...\n'
}

# 2/9/22
# script for performing WHS forced runs with AOA tracers via CAM modifications by Gupta+

DYCORE=$1
RES=$2
NLEV=$3
BUILD_FLAG=$4  # either 0,1,or 2; 
               # 0 throws error if case exists, 
               # 1 will rebuild (without clean-all) and submit,e.g. for namelist changes 
               # 2 will nuke the case and create from scratch
MOD_TYPE=$5    # Options:
               # 0: vanilla HS
               # 1: WHS mod
               # 2: WHS mod + RF sponge
               # 3: WHS mod + RF sponge + topography
               # 4: WHS mod + topography
SFX=$6         # optional suffix to be added to end of case name (useful for ensemble members)

MODEL=/glade/u/home/jhollowed/CAM/CAM_dev/cime/scripts/create_newcase
OUT=/glade/scratch/jhollowed/CAM/cases_589/project4
CASES=/glade/u/home/jhollowed/CAM/cases_589/project4/cases
NL=/glade/u/home/jhollowed/CAM/cases_589/project4/configs
VGRIDS=/glade/u/home/cjablono/CESM_vertical_grids
DATA=/glade/u/home/jhollowed/CAM/inputdata
CONFIGS=/glade/u/home/jhollowed/CAM/cases/aoa_cases/project3/configs
SRCMODS=${CONFIGS}/SourceMods/src.cam.$5

# =============== configure based on dycore, resolution ===============
if [ "$RES" == "C24" ]; then
    if [ "$DYCORE" == "FV3" ]; then 
        GRID="C24_C24_mg17"
        LAB="C24"
        PE=144
    fi 
elif [ "$RES" == "C48" ]; then
    
    if [ "$DYCORE" == "FV3" ]; then 
        GRID="C48_C48_mg17"
        LAB="C48"
        PE=144
        # use namelist settings for either builtin FV3 RF or custom CJ RF
        if [ "$MOD_TYPE" == "4" ]; then
            NL_EXT="_FV3RF"
        elif [ "$MOD_TYPE" == "3" ]; then
            NL_EXT="_CJRF"
        fi
    elif [ "$DYCORE" == "SE" ]; then 
        GRID="ne16_ne16_mg17"
        LAB="ne16"
        PE=288
        NL_EXT=""
    fi
fi

if [ "$NLEV" == "72" ]; then NCDATA="cam_vcoords_L72_E3SM.nc"; fi
if [ "$NLEV" == "93" ]; then NCDATA="cam_vcoords_L93_dz500m_high_top_86km.nc"; fi


STOP_N=720        # total simulated time will be STOP_N * (RESUBMIT+1)
                  # here 2 years --> RESUBMIT=14 for 30 years
RESUBMIT=14

CASENAME=${DYCORE}_${LAB}L${NLEV}_mod${MOD_TYPE}${SFX}
CASE=${CASES}/${CASENAME}


# =============== Create case ===============
if [[ ! -d "$CASE"  ||  $BUILD_FLAG != "0" ]]; then
    
    # ---------- do cleaning if requested
    printf "\n\n========== CLEANING ==========\n"
    if [ "$BUILD_FLAG" != "0" ]; then
        rm -rfv ${CASE}
        rm -rfv ${OUT}/${CASENAME}/run
    fi
    if [ "$BUILD_FLAG" == "2" ]; then
        rm -rfv ${OUT}/${CASENAME}
    fi

    printf "\n\n========== CREATING CASE ==========\n"
    $MODEL --compset FHS94 --res ${GRID} --case $CASE --run-unsupported --project UMIC0087 \
           --pecount $PE --output-root $OUT

    # ---------- configure case
    # nadv_11 = 1 for passive clock tracer
    cd $CASE
    ./xmlchange DEBUG=FALSE,DOUT_S=FALSE,STOP_OPTION=ndays,STOP_N=$STOP_N
    ./xmlchange --file env_build.xml --id CAM_CONFIG_OPTS --val "-phys held_suarez -analytic_ic"
    ./xmlchange --append --file env_build.xml --id CAM_CONFIG_OPTS --val "-age_of_air_trcs -nadv_tt=1"
    ./xmlchange --append --file env_build.xml --id CAM_CONFIG_OPTS --val "--nlev=$NLEV "
    ./xmlchange JOB_WALLCLOCK_TIME=02:00:00
    ./xmlchange SAVE_TIMING=TRUE
    
    # ---------- configure run restarts
    # setting RESUBMIT>0 automatically sets CONTINUE_RUN=TRUE for all runs 
    # subsequent the initial run
    # see: https://www.cesm.ucar.edu/events/tutorials/2018/files/Practical2-shields.pdf
    ./xmlchange RESUBMIT=$RESUBMIT
    ./xmlchange REST_OPTION=ndays
    ./xmlchange REST_N=$STOP_N
    

    printf "\n\n========== COPYING NAMELISTS ==========\n"
    # ---------- copy namelist settings, append vertical levels
    cp --verbose ${NL}/user_nl_cam_aoa_${DYCORE}${NL_EXT} ./user_nl_cam
    sed -i '$a NCDATA = '"\"${VGRIDS}/${NCDATA}\""'' ./user_nl_cam
    
    
    printf "\n\n========== COPYING SOURCEMODS ==========\n"
    # ---------- copy source mods
    cp --verbose ${SRCMODS}/* ./SourceMods/src.cam/
    
    
    printf "\n\n========== CASE SETUP ==========\n"
    ./case.setup
    
    # ---------- build, submit
    printf "\n\n========== BUILDING, SUBMITTING JOB ==========\n"
    #pause
    #./case.build 2>&1 | tee ./log.case.buid 
    qcmd -A UMIC0087 -- ./case.build 2>&1 | tee ./log.case.buid
    ./case.submit 2>&1 | tee ./log.case.submit
else
    printf "\n\n========== CASE EXISTS; ABORTING ==========\n"
fi
