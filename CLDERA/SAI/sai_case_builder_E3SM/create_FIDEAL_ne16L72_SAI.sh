#!/bin/bash

set -e

pause(){
    while read -r -t 0.001; do :; done # dump the buffer
    read -n1 -rsp $'Press any key to continue or Ctrl+C to exit...\n'
}

# 2/9/22
# script for performing E3SM runs with enabled AOA clock tracers

#MACHINE=cori-knl
MACHINE=pm-cpu
COMPILER=intel
PROJECT=m4014
COMPSET=FIDEAL
GRID=ne16pg2_ne16pg2
RES=ne16pg2

BUILD_FLAG=$1
TOT_RUN_LENGTH=$2
SUFFIX=$3

# ----- for debug queueing -----
#PECOUNT=768   # half the number of cubedsphere elements in ne16
PECOUNT=384   # quarter the number of cubedsphere elements in ne16
#QUEUE=debug
#WALLCLOCK=00:30:00
QUEUE=regular
WALLCLOCK=01:30:00

wd="/global/homes/j/jhollo/repos/climate_analysis/CLDERA/SAI/sai_case_builder_E3SM/"
MY_E3SM_ROOT="/global/homes/j/jhollo/E3SM/CLDERA-E3SM"
#MY_E3SM_ROOT="/global/homes/j/jhollo/E3SM/CLDERA-E3SM_SAI"
MODEL="${MY_E3SM_ROOT}/cime/scripts/create_newcase"

#CASES="${wd}/cases"
#CASES="${wd}/cases/pertlim_ens"
#CASES="${wd}/cases/pertlim_ic_ens"
CASES="${wd}/cases/pertlim_ic_ens/for_hsdwpaperfigs_12.7.23/root_case"
CASENAME="HSW_SAI_${RES}_L72_${TOT_RUN_LENGTH}day${SUFFIX}"
CASE=${CASES}/${CASENAME}

#OUTROOT="/global/cscratch1/sd/jhollo/E3SM/E3SMv2_cases/sai_cases"
#OUTROOT="/global/cscratch1/sd/jhollo/E3SM/E3SMv2_cases/sai_cases/pertlim_ens"
OUTROOT="/pscratch/sd/j/jhollo/E3SM/E3SMv2_cases/sai_cases/ic_ens/extra_LIMVAR_outputs_for_hsw_paper/root_case"
RUNDIR="${OUTROOT}/${CASENAME}/run"

DO_RESUBS=false
STOP_N=$TOT_RUN_LENGTH
if [[ "$QUEUE" == "debug" ]]; then
    # if total run length >600 days per run (360 days is ~10 min on PM with 384 ranks), 
    # then require resubmits
    MAX_STOP_N=600
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
       --output-root $OUTROOT --machine $MACHINE --project $PROJECT

# ---------- configure case
cd $CASE
./xmlchange DEBUG=FALSE,DOUT_S=FALSE,STOP_OPTION=ndays,STOP_N=$STOP_N
./xmlchange --append --file env_build.xml --id CAM_CONFIG_OPTS --val "-cldera_sai_trcs "
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

#if [[ $SUFFIX == "noinject" ]]; then
if [[ $SUFFIX == *"_cf" ]]; then
  # counterfactual run with no injection
  DOHEAT=".false."
  OUTVARS="'U','V','T','Z3','OMEGA','PS','T1000','T050','T025'"
  NHTFRQ="-6,-6,-48"
  AVGFLG="'A','I','A'"
  MFILT="720,720,720"
  INITHIST="'ENDOFRUN'"
  DELAY="180"
elif [[ $SUFFIX == *"_pertT" ]]; then
  # counterfactual run with no injection, small T pertubations
  DOHEAT=".false."
  OUTVARS="'U','V','T','Z3','OMEGA','PS','T1000','T050','T025'"
  NHTFRQ="-6,-24"
  AVGFLG="'A','A'"
  MFILT="720,720"
  INITHIST="'DAILY'"
  DELAY="180"
elif [[ $SUFFIX == *"_ptens" ]]; then
  # ic ensemble of small T pertubations
  DOHEAT=".true."
  OUTVARS="'U','V','T','Z3','OMEGA','PS','SO2','ASH','SULFATE','AIR_MASS','AOD','T1000','T050','T025'"
  NHTFRQ="-6,-99999,-24" # suppress instantaneous data
  AVGFLG="'A','I','A'"
  MFILT="720,720,720"
  INITHIST="'ENDOFRUN'"
  DELAY="90"             # this should be the number of days post-pertubation to trigger the injection
elif [[ $SUFFIX == "_mc" ]]; then
  # long-time mean climate run; yearly avgs only
  DOHEAT=".false."
  OUTVARS="'U','V','T','Z3','OMEGA','PS','T1000','T050','T025'"
  NHTFRQ="-8640,-99999,-99999"
  AVGFLG="'A','A','A'"
  MFILT="100,1,1"
  INITHIST="'ENDOFRUN'"
  DELAY="180"
elif [[ $SUFFIX == *"_TEST" ]]; then
  # for random testing...
  DOHEAT=".true."
  OUTVARS="'SAI_HEAT','ATTEN_LW','ATTEN_SW','AOD','SO2','SULFATE','ASH','T1000','T050'"
  NHTFRQ="-6, -48, -48"
  AVGFLG="'A', 'I', 'A'"
  MFILT="720,720,720"
  INITHIST="'ENDOFRUN'"
  DELAY="10"
elif [[ $SUFFIX == *"_extraRunsForHSWPaper" ]]; then
  # for rerunning high-variability root case on PM for HSW paper figures
  DOHEAT=".true."
  OUTVARS="'SAI_HEAT','ATTEN_LW','ATTEN_SW','TROP_P','TROP_Z'"
  NHTFRQ="-9999, -9999, -48"
  AVGFLG="'A', 'I', 'A'"
  MFILT="720,720,720"
  INITHIST="'ENDOFRUN'"
  DELAY="180"
elif [[ $SUFFIX == *"_extraRunsForHSWPaper_LIMVAR" ]]; then
  # for rerunning high-variability root case on PM for HSW paper figures
  DOHEAT=".true."
  OUTVARS="'AOD','U'"
  NHTFRQ="-9999, -9999, -48"
  AVGFLG="'A', 'I', 'A'"
  MFILT="720,720,900"
  INITHIST="'ENDOFRUN'"
  DELAY="90"
else
  # ensemble members
  DOHEAT=".true."
  OUTVARS="'U','V','T','Z3','OMEGA','PS','SO2','ASH','SULFATE','AIR_MASS','AOD','T1000','T050','T025'"
  NHTFRQ="-6,-6,-48"
  AVGFLG="'A','I','A'"
  MFILT="720,720,720"
  INITHIST="'ENDOFRUN'"
  DELAY="180"
fi

printf "\n\n========== POPULATING NAMELIST SETTINGS ==========\n"
cat << EOF >> ./user_nl_eam
empty_htapes     = .TRUE.              ! output only the varibales listed below

! --- for prod runs
NHTFRQ            = ${NHTFRQ}      ! output frequency
avgflag_pertape   = ${AVGFLG}      ! history flags
MFILT             = ${MFILT}       ! history write frequency
                                   ! 720 --> 6-month files at 6hr sampling,
                                   !         4-year files at  2day sampling
cldera_sai_t0     = ${DELAY}       ! injection delay, in days
fincl1 = ${OUTVARS}
fincl2 = ${OUTVARS}
fincl3 = ${OUTVARS}

! --- For tuning 
!NHTFRQ            = -168            ! output frequency every 7 days
!avgflag_pertape   = 'A'             ! hist file 1 is avg, 2 is instant, 3 is avg
!NDENS             = 2               ! single-precision for each hist file
!MFILT             = 1000             ! time samples per hist file
!cldera_sai_t0 = 5                   ! injection starts at 5 days
!fincl1 = 'T','SO2','ASH','SULFATE','SAI_HEAT','AOD','T1000','T975','T050','T025','ATTEN_LW','ATTEN_SW','I_SW','I_LW','COL_SAI','AREA'

inithist=${INITHIST}

! uses IC from HSW 5-year spinup
!NCDATA="/pscratch/sd/j/jhollo/E3SM/E3SMv2_cases/hsw_cases/E3SM_ne16_L72_FIDEAL_10year_spinup/run/E3SM_ne16_L72_FIDEAL_10year_spinup.eam.i.0005-01-01-00000.nc.newCoordNames"
! uses IC of ens05 (FOR LOW-VAR)
NCDATA="/pscratch/sd/j/jhollo/E3SM/E3SMv2_cases/hsw_cases/E3SM_ne16_L72_FIDEAL_3year_ensICs/run/E3SM_ne16_L72_FIDEAL_3year_ensICs.eam.i.0002-05-01-00000.nc.newCoordNames"

! don't let analytic ICs overwrite input from NCDATA
ideal_phys_analytic_ic = .false.

! select HSW idealized forcing (nl added by Ben)
ideal_phys_option = 'held-suarez-williamson'

! turn off any random perturbations (not needed for HS in SE)
pertlim = 0

! activate built-in analytic init of cldera sai tracers
cldera_sai_read_from_ic_file = .false.

! toggle pathways
cldera_sai_formSulfate = .true.
cldera_sai_stratHeating = ${DOHEAT}
cldera_sai_surfCooling = ${DOHEAT}

! heating tuning parameters normalization
cldera_sai_z_cool   = 100
cldera_sai_zeta     = 0.004
cldera_sai_blw_so2  = 0.01 
cldera_sai_blw_sulf = 29
cldera_sai_blw_ash  = 0.00001 
cldera_sai_bsw_so2  = 400
cldera_sai_bsw_sulf = 1900
cldera_sai_bsw_ash  = 400

! injection duration, vertial position
cldera_sai_duration = 24
cldera_sai_z0 = 14

EOF
cat ./user_nl_eam

printf "\n\n========== CASE SETUP ==========\n"
./case.setup

# ---------- build, submit
printf "\n\n========== BUILDING, SUBMITTING JOB ==========\n"
./case.build 2>&1 | tee ./log.case.buid
./case.submit 2>&1 | tee ./log.case.submit
