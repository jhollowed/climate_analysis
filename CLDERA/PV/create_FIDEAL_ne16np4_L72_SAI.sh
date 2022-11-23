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
RES=ne16np4
PECOUNT=384   # quarter the number of cubedsphere elements in ne16

BUILD_FLAG=$1
TOT_RUN_LENGTH=$2

# ----- for debug queueing -----
QUEUE=debug
WALLCLOCK=00:30:00

wd="/global/homes/j/jhollo/repos/climate_analysis/CLDERA/PV/"
MY_E3SM_ROOT="/global/homes/j/jhollo/E3SM/CLDERA-E3SM_PV"
MODEL="${MY_E3SM_ROOT}/cime/scripts/create_newcase"

CASES="${wd}/cases"
CASENAME="E3SM_${RES}_L72_${COMPSET}_PV"
CASE=${CASES}/${CASENAME}

OUTROOT="/global/cscratch1/sd/jhollo/E3SM/E3SMv2_cases/pv_cases"
RUNDIR="${OUTROOT}/${CASENAME}/run"

# if total run length >1/2 year per run (~15 min on Cori with 768 ranks), then require resubmits
MAX_STOP_N=30
DO_RESUBS=false
if [ "$TOT_RUN_LENGTH" -gt "$MAX_STOP_N" ]; then
    STOP_N=$MAX_STOP_N
    DO_RESUBS=true
else
    STOP_N=$TOT_RUN_LENGTH
fi

# if do resubs, compute number of resubs to effectively round up total length of simulation 
# by one unit of STOP_N
RESUBMIT=0
if $DO_RESUBS; then
    RESUBMIT=$(expr $(expr $TOT_RUN_LENGTH - 1) / $STOP_N)
                           # the divide expr rounds down, so RESUBMIT * STOP_N would give 
                           # TOT_RUN_LENGTH - STOP_N
                           # this is correct, since RESUBMIT should be the total number of desired
                           # submissions, minus one for the initial run
fi

echo "STOP_N = $STOP_N"
echo "TOT_RUN_LENGTH = $TOT_RUN_LENGTH"
echo "DO_RESUBS = $DO_RESUBS"
echo "RESUBMIT = $RESUBMIT"

BOLD="\033[1m"
YELLOW="\033[38;5;11m"
RESET="\033[0m"

# -------------- Define case, Build model ---------------

if [[ ! -d "$CASE"  ||  $BUILD_FLAG != "0" ]]; then
    :
else
    printf "\n\n========== CASE EXISTS; ABORTING ==========\n"
fi

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
       --output-root $OUTROOT --machine $MACHINE --compiler $COMPILER --project $PROJECT

# ---------- configure case
cd $CASE
./xmlchange DEBUG=FALSE,DOUT_S=FALSE,STOP_OPTION=ndays,STOP_N=$STOP_N
./xmlchange JOB_WALLCLOCK_TIME=$WALLCLOCK
./xmlchange SAVE_TIMING=TRUE
./xmlchange JOB_QUEUE=$QUEUE

# ---------- configure run restarts
# setting RESUBMIT>0 automatically sets CONTINUE_RUN=TRUE for all runs
# subsequent the initial run
# see: https://www.cesm.ucar.edu/events/tutorials/2018/files/Practical2-shields.pdf
if $DO_RESUBS; then
    ./xmlchange RESUBMIT=$RESUBMIT
fi

printf "\n\n========== POPULATING NAMELIST SETTINGS ==========\n"
cat << EOF >> ./user_nl_eam
empty_htapes     = .TRUE.              ! output only the varibales listed below
avgflag_pertape  = 'I'                 ! hist file 1 is instantaneous

fincl1 = 'U','V','T','PS','PV','DYN_U','DYN_V','DYN_T','DYN_PS','DYN_PV'

NHTFRQ           = -24      ! output frequency every 6 hours
MFILT            = 360      ! allow 360 time samples per hist file (90 days)
NDENS            = 1        ! single-precision for each hist file

inithist='NONE'

! uses new 5-year spinup with extended stratopause
NCDATA="/project/projectdirs/m4014/data/HSW/initial_conditions/netcdf/E3SM_ne16_L72_FIDEAL_10year_spinup.eam.i.0005-01-01-00000.nc.newCoordNames"

! don't let analytic ICs overwrite input from NCDATA
ideal_phys_analytic_ic = .false.

! select HSW idealized forcing (nl added by Ben)
ideal_phys_option = 'held-suarez-williamson'

! turn off any random perturbations (not needed for HS in SE)
pertlim = 0
EOF
cat ./user_nl_eam

printf "\n\n========== CASE SETUP ==========\n"
./case.setup

# ---------- build, submit
printf "\n\n========== BUILDING, SUBMITTING JOB ==========\n"
./case.build 2>&1 | tee ./log.case.buid
./case.submit 2>&1 | tee ./log.case.submit
