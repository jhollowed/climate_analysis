#!/bin/bash


# --- hsw plume runs with builtin SAI
# usage:
# BUILD_FLAG=$1      (0 = abort if case exists
#                     1 = remove run files but don't build clean
#                     2 = rm case directory and create from scratch)
# TOT_RUN_LENGTH=$2  (in days)
# SUFFIX=$3          (substring to tag onto the end of case name, optional)
# CUSTOMNL=$4        (namelist file to use if not default, optional)

#./create_FIDEAL_ne16pg2L72_SAI.sh 1 180
./create_FIDEAL_ne16pg2L72_SAI.sh 1 10 'CLDERA_TEST'

