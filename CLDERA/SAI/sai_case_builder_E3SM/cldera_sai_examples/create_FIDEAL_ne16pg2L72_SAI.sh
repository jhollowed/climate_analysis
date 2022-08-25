#!/bin/bash

# --------------------------------------------------------------------------------------
# Joe Hollowed 8/4/22
# Script for creating E3SM cases with HSW++ climatology and tracer injection
# 
# Usage:
# ./create_HSW_SAI_case.sh [RUN LENGTH] [GRID DEFINITION] [CASE SUFFIX] [NAMELIST FILE PATH]
# 
# Arguments
# ---------
#   [RUN LENGTH] : float
#       The total run lenth in units of days.
#       Defaults to 30
#   [GRID DEFINITION] : string, optional
#       Full name of the desired ATM grid. See available options via:
#       >> E3SM/cime/scripts/query_config --grids
#       Defaults to ne16pg2_ne16pg2 (~2 degree or 200km SE subed sphere grid with 2x2 physics), 
#       described here: https://e3sm.org/new-physgrid-and-dycore-methods-speed-up-eam-by-2x/
#   [CASE SUFFIX] : string, optional
#       Substring to append to the end of the case name for identification
#       Defaults to an empty string ''
#   [NAMELIST FILE PATH] : string, optional
#       Location of namelist file to copy into case directorty. 
#       Defaults to ./user_nl_eam, where '.' is the current working directory
#
# Note: Make sure to update personal paths to your own locations as noted below
#
# --------------------------------------------------------------------------------------


set -e
wd=$(cd $(dirname $0) && pwd)

# ========== CHANGE THESE TO PERSONAL PATHS ==========

# This points to the repo clone you wish to use; currently defaulting to my HSW++ branch at
# https://github.com/sandialabs/CLDERA-E3SM/tree/jhollowed/eam/cldera-sai-module
#MY_E3SM_ROOT="/global/homes/j/jhollo/E3SM/CLDERA-E3SM_SAIBranch"
MY_E3SM_ROOT="/global/homes/j/jhollo/E3SM/CLDERA-E3SM_SAIBranch_PV"
# This is the parent directory where new case directories created by this script will be placed
CASE_ROOT="/global/homes/j/jhollo/repos/climate_analysis/CLDERA/SAI/sai_case_builder_E3SM/cldera_sai_examples/cases"
# This is the parent directory where new case outputs and build files will be placed
#OUT_ROOT="/global/cscratch1/sd/jhollo/E3SM/E3SMv2_cases/sai_cases/cldera_template_cases"

# =====================================================

MACHINE=cori-knl
COMPILER=intel
PROJECT=m4014
COMPSET=FIDEAL
TOT_RUN_LENGTH=$1
GRID=$2
SUFFIX=$3
NL=$4

# ------ set argument defaults ------
if [ -z "$GRID" ]; then
    #GRID=ne16pg2_ne16pg2
    GRID=ne16_ne16
fi
if [ -z "$TOT_RUN_LENGTH" ]; then
    TOT_RUN_LENGTH=30
fi
if [ -z "$SUFFIX" ]; then
    SUFFIX=""
else
    SUFFIX="_${SUFFIX}"
fi
if [ -z "$NL" ]; then
    NL="${wd}/user_nl_eam"
fi

# ------ for debug queueing ------
PECOUNT=384   # quarter the number of cubedsphere elements in ne16
QUEUE=debug
WALLCLOCK=00:30:00

# ------ point to CIME case creation script ------
MODEL="${MY_E3SM_ROOT}/cime/scripts/create_newcase"

# ------ create case, build model ------
CASENAME="E3SM_${GRID}_L72_${COMPSET}_SAI${SUFFIX}"
CASE=${CASE_ROOT}/${CASENAME}
RUNDIR="${OUT_ROOT}/${CASENAME}/run"
OUT_ROOT="${CASE}/OUTPUT"

printf "\n\n========== CREATING CASE ==========\n"
$MODEL --compset $COMPSET --res $GRID --case $CASE --pecount $PECOUNT \
       --output-root $OUT_ROOT --machine $MACHINE --compiler $COMPILER --project $PROJECT

# ---------- configure case
cd $CASE
./xmlchange DEBUG=FALSE,DOUT_S=FALSE,STOP_OPTION=ndays,STOP_N=$TOT_RUN_LENGTH
./xmlchange JOB_WALLCLOCK_TIME=$WALLCLOCK
./xmlchange SAVE_TIMING=TRUE
./xmlchange JOB_QUEUE=$QUEUE

# Here we update the ATM model's "config options", which are a set of options different from
# the case configration options set directly by xmlchange, and different from namelist settings.
# See qualitative description here (E3SM will not agree with all details on this page):
# https://ncar.github.io/CAM/doc/build/html/users_guide/customizing-compsets.html
# 
# 'cldera_sai_trcs' is a flag which activates the idealized injection. If this is ommited, 
# the tracer fields will never be initialized
# 'nlev' gives the number or vertical levels, which triggers the model to read a specific 
# file for the hybrid coefficients if a vertical discretization for this nlev exists. Should
# default to 72 for most grids on most compsets; be explicit to be safe
./xmlchange --append --file env_build.xml --id CAM_CONFIG_OPTS --val "-cldera_sai_trcs -nlev 72"

# Here we set the physics timestep, in number of timesteps per simulated day. Currently
# set to the ne16pg2_ne16pg2 default of ATM_NCPL=48, giving a 30 minute physics timestep.
# The dynamics timestep is typically shorter than this; see comments in provided 'user_nl_eam',
# and some info here:
# https://acme-climate.atlassian.net/wiki/spaces/DOC/pages/1044644202/EAM+s+HOMME+dycore
./xmlchange ATM_NCPL=48

printf "\n\n========== COPYING NAMELISTS ==========\n"
cp --verbose $NL ./user_nl_eam

printf "\n\n========== CASE SETUP ==========\n"
./case.setup

# ---------- build, submit
printf "\n\n========== BUILDING, SUBMITTING JOB ==========\n"
./case.build 2>&1 | tee ./log.case.buid
./case.submit 2>&1 | tee ./log.case.submit
