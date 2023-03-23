#!/bin/bash

set -e
wd=$(cd $(dirname $0) && pwd)

MACHINE=cori-knl
COMPILER=intel
PROJECT=m4014
COMPSET=FIDEAL
GRID=ne16pg2_ne16pg2
RES=ne16

BUILD_FLAG=$1
TOT_RUN_LENGTH=$2

# ----- for debug queueing -----
PECOUNT=768   # half the number of cubedsphere elements in ne16
QUEUE=debug
WALLCLOCK=00:30:00

wd="/global/homes/j/jhollo/repos/climate_analysis/CLDERA/PV"
MY_E3SM_ROOT="/global/homes/j/jhollo/E3SM/CLDERA-E3SM_PV"
MODEL="${MY_E3SM_ROOT}/cime/scripts/create_newcase"

CASES="${wd}/cases"
CASENAME="E3SM_${RES}_L72_${COMPSET}_dynTrac"
CASE=${CASES}/${CASENAME}

OUTROOT="/global/cscratch1/sd/jhollo/E3SM/E3SMv2_cases/pv_cases"
RUNDIR="${OUTROOT}/${CASENAME}/run"

DO_RESUBS=false
STOP_N=$TOT_RUN_LENGTH
if [[ "$QUEUE" == "debug" ]]; then
    # if total run length >300 days per run (~20 min on Cori with 768 ranks), then require resubmits
    MAX_STOP_N=300
    if [ "$TOT_RUN_LENGTH" -gt "$MAX_STOP_N" ]; then
        STOP_N=$MAX_STOP_N
        DO_RESUBS=true
    fi
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
./xmlchange --append --file env_build.xml --id CAM_CONFIG_OPTS --val "-cldera_dynamic_trcs "
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

! output a few quantities for frequently on their instantaneous values for SSW quantification
fincl1 = 'PV','PV_TRCR','PT_TRCR'

avgflag_pertape  = 'A'               ! hist file 1 is avg, 2 is instant, 3 is avg
NHTFRQ            = -6               ! output frequency every day
MFILT             = 720              ! allow 720 time samples per hist file
NDENS             = 2                ! single-precision for each hist file

inithist='ENDOFRUN'

! uses IC from HSW 5-year spinup
NCDATA="/global/cscratch1/sd/jhollo/E3SM/E3SMv2_cases/hsw_cases/E3SM_ne16_L72_FIDEAL_10year_spinup/run/E3SM_ne16_L72_FIDEAL_10year_spinup.eam.i.0005-01-01-00000.nc.newCoordNames"

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
