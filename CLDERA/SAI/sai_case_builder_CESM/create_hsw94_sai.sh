#!/bin/bash

set -e

pause(){
    while read -r -t 0.001; do :; done # dump the buffer
    read -n1 -rsp $'Press any key to continue or Ctrl+C to exit...\n'
}

# 2/9/22
# script for performing WHS forced runs with AOA tracers via CAM modifications by Gupta+

BUILD_FLAG=$1  # either 0,1,or 2; 
               # 0 throws error if case exists, 
               # 1 will rebuild (without clean-all) and submit,e.g. for namelist changes 
               # 2 will nuke the case and create from scratch
FIX=$2      # turn on mass fix?
SPINUP=$3   # is spinup?
STOP_N=$4

TAUPHYS=$5  # injection timescale to match phys timestep?
NSPLIT=$6   # nsplit
NODIFF=$7   # turn off diffusion?

MODEL=/glade/u/home/jhollowed/CAM/CAM_dev/cime/scripts/create_newcase
OUT=/glade/scratch/jhollowed/CAM/cases/sai_runs
CASES=/glade/u/home/jhollowed/repos/climate_analysis/CLDERA/SAI/sai_case_builder_CESM/cases
CONFIGS=/glade/u/home/jhollowed/repos/climate_analysis/CLDERA/SAI/sai_case_builder_CESM/configs
VGRIDS=/glade/u/home/cjablono/CESM_vertical_grids
DATA=/glade/u/home/jhollowed/CAM/inputdata
# THIS NEEDED FOR EXTENDED SPINUP; UPDATE THIS TO POINT TO 2 YR OUTPUT
SPUN=/glade/scratch/jhollowed/CAM/cases/sai_runs/SE_ne16L72_whs_sai_spinup_FIRSTYEAR/run/SE_ne16L72_whs_sai_spinup.cam.i.0001-12-27-00000.nc

GRID="ne16_ne16_mg17"
COMPSET="FHS94"
LEVFILE="cam_vcoords_L72_E3SM.nc"
PROJECT="UMIC0087"
PE=288


if [ "$SPINUP" == '1' ]; then
    CASENAME=SE_ne16L72_whs_sai_spinup
    # DISABLING TEMPORARILY TO EXTEND SPINUP
    #IC='-analytic_ic'
    IC=''
fi
if [ "$SPINUP" == '0' ]; then
    CASENAME=SE_ne16L72_whs_sai_fix${FIX}_tau${TAUPHYS}_nsplit${NSPLIT}_nodiff${NODIFF}
    IC=''
fi
CASE=${CASES}/${CASENAME}

if [ "$FIX" == '1' ]; then
    SRCMODS=${CONFIGS}/SourceMods_massCorrect
fi
if [ "$FIX" == '0' ]; then
    SRCMODS=${CONFIGS}/SourceMods
fi


# =============== Create case ===============
if [[ ! -d "$CASE"  ||  $BUILD_FLAG != "0" ]]; then
    
    # ---------- do cleaning if requested
    printf "\n\n========== CLEANING ==========\n"
    if [ "$BUILD_FLAG" != "0" ]; then
        read -p "Press Y/y to CLEAN case ${CASE}, INCLUDING OUTPUT" -n 1 -r
        echo
        if [[ $REPLY =~ ^[Yy]$ ]]; then
            rm -rfv ${CASE}
            rm -rfv ${OUT}/${CASENAME}/run
        fi
    fi
    if [ "$BUILD_FLAG" == "2" ]; then
        rm -rfv ${OUT}/${CASENAME}
    fi

    printf "\n\n========== CREATING CASE ==========\n"
    $MODEL --compset $COMPSET --res $GRID --case $CASE --run-unsupported --project $PROJECT \
           --pecount $PE --output-root $OUT

    # ---------- configure case
    # nadv_11 = 1 for passive clock tracer
    cd $CASE
    ./xmlchange DEBUG=FALSE,DOUT_S=FALSE,STOP_OPTION=ndays,STOP_N=$STOP_N
    ./xmlchange --file env_build.xml --id CAM_CONFIG_OPTS --val "-phys held_suarez $IC"
    ./xmlchange --append --file env_build.xml --id CAM_CONFIG_OPTS --val "-age_of_air_trcs "
    ./xmlchange --append --file env_build.xml --id CAM_CONFIG_OPTS --val "--nlev=72 "
    ./xmlchange JOB_WALLCLOCK_TIME=02:00:00
    ./xmlchange SAVE_TIMING=TRUE
    
    printf "\n\n========== COPYING NAMELISTS ==========\n"
    # ---------- copy namelist settings, append vertical levels
    if [ "$SPINUP" == '1' ]; then
        cp --verbose ${CONFIGS}/user_nl_cam_aoa_SE_spinup ./user_nl_cam
        # TEMPORAY for continuing spinup for extra year
        #sed -i '$a NCDATA = '"\"${VGRIDS}/${LEVFILE}\""'' ./user_nl_cam
        sed -i '$a NCDATA = '"\"${SPUN}\""'' ./user_nl_cam
    fi
    if [ "$NODIFF" == '1' ]; then
        cp --verbose ${CONFIGS}/user_nl_cam_aoa_SE_massrun ./user_nl_cam
        sed -i '$a NCDATA = '"\"${SPUN}\""'' ./user_nl_cam
    fi
    if [ "$SPINUP" == '0' ]; then
        if [ "$NODIFF" == '0' ]; then
            cp --verbose ${CONFIGS}/user_nl_cam_aoa_SE ./user_nl_cam
            sed -i '$a NCDATA = '"\"${SPUN}\""'' ./user_nl_cam
        fi
    fi
    
    sed -i '$a se_nsplit = '"${NSPLIT}"'' ./user_nl_cam
    
    printf "\n\n========== COPYING SOURCEMODS ==========\n"
    # ---------- copy source mods
    cp --verbose ${SRCMODS}/* ./SourceMods/src.cam/
    if [ "$TAUPHYS" == '1' ]; then
        if [ "$FIX" == '1' ]; then
            sed -i '346s/.*/    tau    = 1.0\/dt/' ./SourceMods/src.cam/aoa_tracers.F90
        fi
        if [ "$FIX" == '0' ]; then
            sed -i '352s/.*/    tau    = 1.0\/dt/' ./SourceMods/src.cam/aoa_tracers.F90
        fi
    fi
    if [ "$NODIFF" == '1' ]; then
        if [ "$FIX" == '0' ]; then
            sed -i '349s/.*/    k_so2    = 0.0/' ./SourceMods/src.cam/aoa_tracers.F90
            sed -i '350s/.*/    k_ash    = 0.0/' ./SourceMods/src.cam/aoa_tracers.F90
        fi
        if [ "$FIX" == '1' ]; then
            sed -i '355s/.*/    k_so2    = 0.0/' ./SourceMods/src.cam/aoa_tracers.F90
            sed -i '356s/.*/    k_ash    = 0.0/' ./SourceMods/src.cam/aoa_tracers.F90
        fi
    fi
    
    printf "\n\n========== CASE SETUP ==========\n"
    ./case.setup
    
    # ---------- build, submit
    printf "\n\n========== BUILDING, SUBMITTING JOB ==========\n"
    ./case.build 2>&1 | tee ./log.case.buid
    ./case.submit 2>&1 | tee ./log.case.submit
else
    printf "\n\n========== CASE EXISTS; ABORTING ==========\n"
fi
