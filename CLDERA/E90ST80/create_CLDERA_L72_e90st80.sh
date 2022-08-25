#!/bin/bash

# --------------------------------------------------------------------------------------
# Joe Hollowed 8/4/22
# Script for creating E3SM cases with E90, ST80 tracer advection
# 
# Usage:
# ./create_CLDERA_L72_e90st80.sh [TRACER IMPLEMENTATION] [RUN LENGTH] [GRID DEFINITION] [COMPSET]
# 
# Arguments
# ---------
#   [TRACER IMPLEMENTATION] : integer, optional
#       Which E90, ST80 tracer implementation to use
#       1 : the modular (JH) implementation
#       2 : the WACCM-like (AH) implementation
#       Defaults to option 1
#   [RUN LENGTH] : float, optional
#       The total run lenth in units of days.
#       Defaults to 10
#   [GRID DEFINITION] : string, optional
#       Full name of the desired horizontal grid. See available options via:
#       >> E3SM/cime/scripts/query_config --grids
#       Defaults to ne30pg2_ne30pg2 (~1 degree or 100km SE subed sphere grid with 2x2 physics), 
#       described here: https://e3sm.org/new-physgrid-and-dycore-methods-speed-up-eam-by-2x/
#   [COMPSET] : string, optional
#        Name of desired compset. See available options via
#        >> E3SM/cime/scripts/query_config --compsets
#        Defaults to F1850 (AMIP with MPAS sea ice)
#
#   For F1850, the MPAS Sea Ice model is used, which requires the following horizontal grid at
#   e.g. ne30: 
#   ne30pg2_EC30to60E2r2
#   If needing to use a different resolution for which an MPAS grid is not available, e.g.
#   ne4pg2_ne4pg2, must use F1850-CICE as the compset
#
# Note: Make sure to update personal paths to your own locations as noted below
#
# --------------------------------------------------------------------------------------

set -e
wd=$(cd $(dirname $0) && pwd)



# ========== CHANGE THESE VALUES! ====================

# This is the parent directory where new case directories created by this script will be placed
CASE_ROOT="/global/homes/j/jhollo/repos/climate_analysis/CLDERA/E90ST80/cases"
# This is the parent directory where new case outputs and build files will be placed
OUT_ROOT="/global/cscratch1/sd/jhollo/E3SM/E3SMv2_cases/e90st80_cases"

# =====================================================

printf "\n\n========== SETTING VARIABLES ==========\n"

MACHINE=cori-knl
COMPILER=intel
PROJECT=m4014
TRACER_IMPL=$1
TOT_RUN_LENGTH=$2
GRID=$3
COMPSET=$4

# ------ set argument defaults ------
if [ -z "$TRACER_IMPL" ]; then
    TRACER_IMPL=1
fi
if [ -z "$TOT_RUN_LENGTH" ]; then
    TOT_RUN_LENGTH=10
fi
if [ -z "$GRID" ]; then
    GRID=ne30pg2_EC30to60E2r2
fi
if [ -z "$COMPSET" ]; then
    COMPSET=F1850
fi
printf "Grid: ${GRID}\n"
printf "Compset: ${COMPSET}\n"
printf "Run length in days: ${TOT_RUN_LENGTH}\n"

# ------ configure tracer implementation ------
if [ "$TRACER_IMPL" -eq 1 ]; then
    # point to Joe's branch
    CONFIG_APPEND_VAL="-cldera_passive_trcs "
    MY_E3SM_ROOT="/global/homes/j/jhollo/E3SM/CLDERA-E3SM_passive"
elif [ "$TRACER_IMPL" -eq 2 ]; then
    # point to Allen's branch
    CONFIG_APPEND_VAL="-usr_mech_infile /global/homes/a/allenhu/chem_mech_bothtracers.in "
    MY_E3SM_ROOT="/global/homes/j/jhollo/E3SM/CLDERA-E3SM_AllenE90ST80"
fi
printf "Using tracer implementation ${TRACER_IMPL}\n"
printf "with config option ${CONFIG_APPEND_VAL}\n"
printf "on branch ${MY_E3SM_ROOT}\n"

# ------ for debug queueing ------
PECOUNT=240         # this must match an available mpassi pelayout, seen here:
                    # https://web.lcrc.anl.gov/public/e3sm/inputdata/ice/mpas-seaice/EC30to60E2r2/
QUEUE=debug         # run on debug queue for faster queue times
WALLCLOCK=00:30:00  # max debug run time
printf "${PECOUNT} procs on queue ${QUEUE} for ${WALLCLOCK}\n"

# ------ point to CIME case creation script ------
MODEL="${MY_E3SM_ROOT}/cime/scripts/create_newcase"

# ------ create case, build model ------
CASENAME="E3SM_${COMPSET}_${GRID}_L72_e90st80_impl${TRACER_IMPL}"
CASE=${CASE_ROOT}/${CASENAME}
RUNDIR="${OUT_ROOT}/${CASENAME}/run"

# ------ automatically compute number of resubmits ------
# if total run length 30 days, then require resubmits to avoid running over 30 min on debug queue
# (this value may need to be tuned for ne30pg2)
MAX_STOP_N=31
DO_RESUBS=false
if [ "$TOT_RUN_LENGTH" -gt "$MAX_STOP_N" ]; then
    STOP_N=$MAX_STOP_N
    DO_RESUBS=true
else
    STOP_N=$TOT_RUN_LENGTH
    DO_RESUBS=false
fi

# if do resubs, compute number of resubs to effectively round up total length of simulation 
# by one unit of STOP_N
RESUBMIT=0
if $DO_RESUBS; then
    RESUBMIT=$(expr $(expr $TOT_RUN_LENGTH - 1) / $STOP_N)
               # the divide expr rounds down, so RESUBMIT * STOP_N would give TOT_RUN_LENGTH - STOP_N
               # this is correct, since RESUBMIT should be the total number of desired submissions, 
               # minus one for the initial run
fi
printf "With max run length of ${MAX_STOP_N} days, \
run length of ${TOT_RUN_LENGTH} requires ${RESUBMIT} resubmissions\n"



printf "\n\n========== CREATING CASE ==========\n"
$MODEL --compset $COMPSET --res $GRID --case $CASE --pecount $PECOUNT \
       --output-root $OUT_ROOT --machine $MACHINE --compiler $COMPILER --project $PROJECT

# ---------- configure case
cd $CASE
./xmlchange DEBUG=FALSE,DOUT_S=FALSE,STOP_OPTION=ndays,STOP_N=$TOT_RUN_LENGTH
./xmlchange JOB_WALLCLOCK_TIME=$WALLCLOCK
./xmlchange SAVE_TIMING=TRUE
./xmlchange JOB_QUEUE=$QUEUE
./xmlchange --append --file env_build.xml --id CAM_CONFIG_OPTS --val "$CONFIG_APPEND_VAL"
    
# ---------- configure run restarts
# setting RESUBMIT>0 automatically sets CONTINUE_RUN=TRUE for all runs subsequent the initial run
# see: https://www.cesm.ucar.edu/events/tutorials/2018/files/Practical2-shields.pdf
if $DO_RESUBS; then
    ./xmlchange RESUBMIT=$RESUBMIT
fi



printf "\n\n========== POPULATING NAMELIST SETTINGS ==========\n"
cat << EOF >> ./user_nl_eam
empty_htapes     = .TRUE.              ! output only the varibales listed below
! output a few quantities for frequently on their instantaneous values for SSW quantification
fincl1 = 'U','V','T','Z3','OMEGA','PS','E90','ST80_25'

avgflag_pertape   = 'A'             ! hist file 1 is avg
NHTFRQ            = -24             ! output frequency every day
MFILT             = 360             ! allow 360 time samples per hist file
NDENS             = 2               ! single-precision for each hist file

inithist='ENDOFRUN'

! turn off any random perturbations
pertlim = 0
EOF
cat ./user_nl_eam


printf "\n\n========== CASE SETUP ==========\n"
./case.setup



# ---------- build, submit
printf "\n\n========== BUILDING, SUBMITTING JOB ==========\n"
./case.build 2>&1 | tee ./log.case.buid
./case.submit 2>&1 | tee ./log.case.submit
