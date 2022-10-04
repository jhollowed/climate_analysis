#!/bin/bash

set -e
wd=$(cd $(dirname $0) && pwd)


# =====================================================

printf "\n\n========== SETTING VARIABLES ==========\n"

MACHINE=cori-knl
COMPILER=intel
PROJECT=m4014
COMPSET=F1850
GRID=ne4pg2_oQU480
BUILD_FLAG=$1
TOT_RUN_LENGTH=$2

printf "Grid: ${GRID}\n"
printf "Compset: ${COMPSET}\n"
printf "Run length in days: ${TOT_RUN_LENGTH}\n"

# ------ configure tracers ------
CONFIG_APPEND_VAL="-cldera_passive_trcs "
MY_E3SM_ROOT="/global/homes/j/jhollo/E3SM/CLDERA-E3SM_passive"
MODEL="${MY_E3SM_ROOT}/cime/scripts/create_newcase"

CASES="${wd}/cases"
CASENAME="E3SM_${GRID}_L72_${COMPSET}_e90st80"
CASE=${CASES}/${CASENAME}
RUNDIR="${OUT_ROOT}/${CASENAME}/run"
OUT_ROOT="/global/cscratch1/sd/jhollo/E3SM/E3SMv2_cases/e90st80_cases"

# ------ for debug queueing ------
PECOUNT=S           # this must match an available mpassi pelayout if using F1850, seen here:
                    # https://web.lcrc.anl.gov/public/e3sm/inputdata/ice/mpas-seaice/EC30to60E2r2/
                    # https://web.lcrc.anl.gov/public/e3sm/inputdata/ice/mpas-cice/oQU480/
QUEUE=debug         # run on debug queue for faster queue times
WALLCLOCK=00:30:00  # max debug run time
printf "${PECOUNT} procs on queue ${QUEUE} for ${WALLCLOCK}\n"

# ------ point to CIME case creation script ------

# ------ automatically compute number of resubmits for debug ------
DO_RESUBS=false
if [ "$QUEUE" = "debug" ]; then
    # if total run length 180 days, then require resubmits to avoid running over 30 min on debug queue
    # (this value will need to be tuned for ne30pg2 to avoid timeouts; 
    #  was originally tuned to ne16 FIDEAL w/ pecount=384)
    MAX_STOP_N=181
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
        # the divide expr rounds down, so RESUBMIT * STOP_N would give TOT_RUN_LENGTH - STOP_N
        # this is correct, since RESUBMIT should be the total number of desired submissions, 
        # minus one for the initial run
    fi
    printf "With max run length of ${MAX_STOP_N} days, \
    run length of ${TOT_RUN_LENGTH} requires ${RESUBMIT} resubmissions\n"
fi

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
fincl1 = 'U','V','T','Z3','OMEGA','PS','E90j','ST80_25j','AOA'

avgflag_pertape   = 'A'             ! hist file 1 is avg
NHTFRQ            = -24             ! output frequency every day
MFILT             = 360             ! allow 360 time samples per hist file
NDENS             = 2               ! single-precision for each hist file

inithist='ENDOFRUN'
EOF
cat ./user_nl_eam

printf "\n\n========== CASE SETUP ==========\n"
./case.setup

# ---------- build, submit
printf "\n\n========== BUILDING, SUBMITTING JOB ==========\n"
./case.build 2>&1 | tee ./log.case.buid
./case.submit 2>&1 | tee ./log.case.submit
