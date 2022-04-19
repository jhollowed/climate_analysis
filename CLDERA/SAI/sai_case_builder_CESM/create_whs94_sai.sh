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
TAUPHYS=$5  # injection timescale to match phys timestep?
FIX=$6      # turn on mass fix?
NSPLIT=$7   # nsplit
SPINUP=$8   # is spinup?
NODIFF=$9   # turn off diffusion?

MODEL=/glade/u/home/jhollowed/CAM/CAM_dev/cime/scripts/create_newcase
OUT=/glade/scratch/jhollowed/CAM/cases/sai_runs
CASES=/glade/u/home/jhollowed/repos/climate_analysis/CLDERA/SAI/sai_case_builder_CESM/cases
CONFIGS=/glade/u/home/jhollowed/repos/climate_analysis/CLDERA/SAI/sai_case_builder_CESM/configs
VGRIDS=/glade/u/home/cjablono/CESM_vertical_grids
DATA=/glade/u/home/jhollowed/CAM/inputdata
SPUN=/glade/scratch/jhollowed/CAM/cases/sai_runs/SE_ne16L72_whs_saiv2_spinup/run/SE_ne16L72_whs_saiv2_spinup.cam.i.0001-12-27-00000.nc

# =============== configure based on dycore, resolution ===============
if [ "$RES" == "C24" ]; then
    if [ "$DYCORE" == "FV3" ]; then 
        GRID="C24_C24_mg17"
        LAB="C24"
        PE=144
        #ATM_NCPL=48
        ATM_DOMAIN_MESH="/glade/p/cesmdata/cseg/inputdata/share/scripgrids/C24_SCRIP_desc.181018.nc"
    fi 
elif [ "$RES" == "C48" ]; then
    if [ "$DYCORE" == "FV3" ]; then 
        GRID="C48_C48_mg17"
        LAB="C48"
        PE=144
        ATM_DOMAIN_MESH="/glade/p/cesmdata/cseg/inputdata/share/scripgrids/C48_SCRIP_desc.181018.nc"
    elif [ "$DYCORE" == "SE" ]; then 
        GRID="ne16_ne16_mg17"
        LAB="ne16"
        PE=288
        #ATM_NCPL=96
    fi
fi

if [ "$NLEV" == "72" ]; then NCDATA="cam_vcoords_L72_E3SM.nc"; fi
if [ "$NLEV" == "93" ]; then NCDATA="cam_vcoords_L93_dz500m_high_top_86km.nc"; fi



if [ "$SPINUP" == '1' ]; then
    CASENAME=${DYCORE}_${LAB}L${NLEV}_whs_saiv2_spinup
    STOP_N=360
fi
if [ "$SPINUP" == '0' ]; then
    CASENAME=${DYCORE}_${LAB}L${NLEV}_whs_saiv2_fix${FIX}_tau${TAUPHYS}_qsplit${NSPLIT}
    STOP_N=60
    if [ "$NODIFF" == '1' ]; then
        CASENAME=${DYCORE}_${LAB}L${NLEV}_whs_saiv2_fix${FIX}_tau${TAUPHYS}_qsplit${NSPLIT}_NODIFF
        STOP_N=4
    fi
fi
CASE=${CASES}/${CASENAME}

if [ "$FIX" == '1' ]; then
    SRCMODS=${CONFIGS}/SourceMods_massCorrect
    TAULINE=348
fi
if [ "$FIX" == '0' ]; then
    SRCMODS=${CONFIGS}/SourceMods
    TAULINE=342
fi


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
    ./xmlchange --append --file env_build.xml --id CAM_CONFIG_OPTS --val "-age_of_air_trcs "
    ./xmlchange --append --file env_build.xml --id CAM_CONFIG_OPTS --val "--nlev=$NLEV "
    ./xmlchange JOB_WALLCLOCK_TIME=02:00:00
    ./xmlchange SAVE_TIMING=TRUE
    
    printf "\n\n========== COPYING NAMELISTS ==========\n"
    # ---------- copy namelist settings, append vertical levels
    if [ "$SPINUP" == '1' ]; then
        cp --verbose ${CONFIGS}/user_nl_cam_aoa_${DYCORE}_spinup ./user_nl_cam
        sed -i '$a NCDATA = '"\"${VGRIDS}/${NCDATA}\""'' ./user_nl_cam
    fi
    if [ "$NODIFF" == '1' ]; then
        cp --verbose ${CONFIGS}/user_nl_cam_aoa_${DYCORE}_massrun ./user_nl_cam
        sed -i '$a NCDATA = '"\"${SPUN}\""'' ./user_nl_cam
    fi
    if [ "$SPINUP" == '0' ]; then
        if [ "$NODIFF" == '0' ]; then
            cp --verbose ${CONFIGS}/user_nl_cam_aoa_${DYCORE} ./user_nl_cam
            sed -i '$a NCDATA = '"\"${SPUN}\""'' ./user_nl_cam
        fi
    fi
    
    sed -i '$a se_nsplit = '"${NSPLIT}"'' ./user_nl_cam
    
    printf "\n\n========== COPYING SOURCEMODS ==========\n"
    # ---------- copy source mods
    cp --verbose ${SRCMODS}/* ./SourceMods/src.cam/
    if [ "$TAUPHYS" == '1' ]; then
        if [ "$FIX" == '1' ]; then
            sed -i '348s/.*/    tau    = 1.0\/dt/' ./SourceMods/src.cam/aoa_tracers.F90
        fi
        if [ "$FIX" == '0' ]; then
            sed -i '342s/.*/    tau    = 1.0\/dt/' ./SourceMods/src.cam/aoa_tracers.F90
        fi
    fi
    if [ "$NODIFF" == '1' ]; then
        if [ "$FIX" == '0' ]; then
            sed -i '345s/.*/    k_so2    = 0.0/' ./SourceMods/src.cam/aoa_tracers.F90
            sed -i '346s/.*/    k_ash    = 0.0/' ./SourceMods/src.cam/aoa_tracers.F90
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
