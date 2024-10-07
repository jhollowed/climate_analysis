#!/bin/bash

# specify ensemble member number
ens=$1
mass=$2
t1=0
t2=36

# this should be a string contining the steps to perform, e.g.
# 12345678910
# or
# 34567
steps=$3

# whether or not to overwrite the data for the processing steps specified
overwrite=$4

# whether or not to do a "dry run"
dry=$5


# step 1
if [[ $steps == *"1"* ]]; then
    echo
    echo "===================="
    echo "====== STEP 1 ======"
    echo "===================="
    echo
    python ./limvar_zinterp.py $ens 0 $mass $t1 $t2 $overwrite $dry # monthly
    python ./limvar_zinterp.py $ens 1 $mass $t1 $t2 $overwrite $dry # daily
fi

# step 2
if [[ $steps == *"2"* ]]; then
    echo
    echo "===================="
    echo "====== STEP 2 ======"
    echo "===================="
    echo
    python ./limvar_TEM.py $ens $mass $t1 $t2 $overwrite $dry
fi

# step 3
if [[ $steps == *"3"* ]]; then
    echo
    echo "===================="
    echo "====== STEP 3 ======"
    echo "===================="
    echo
    python ./limvar_zonalmeans.py $ens 0 $mass $t1 $t2 $overwrite $dry # monthly
    python ./limvar_zonalmeans.py $ens 0 $mass $t1 $t2 $overwrite $dry 1 # monthly tropopause
    python ./limvar_zonalmeans.py $ens 1 $mass $t1 $t2 $overwrite $dry # daily
    python ./limvar_zonalmeans.py $ens 1 $mass $t1 $t2 $overwrite $dry 1 # daily tropopause
fi

#step 4
if [[ $steps == *"4"* ]]; then
    echo
    echo "===================="
    echo "====== STEP 4 ======"
    echo "===================="
    echo
    python ./limvar_timeConcat.py $ens 0 $mass $t1 $t2 $overwrite $dry # monthly
    python ./limvar_timeConcat.py $ens 1 $mass $t1 $t2 $overwrite $dry # daily
fi

#step 5
if [[ $steps == *"5"* ]]; then
    echo
    echo "===================="
    echo "====== STEP 5 ======"
    echo "===================="
    echo
    python ./limvar_monthly_seasonal_means.py $ens $mass $t1 $t2 $overwrite $dry
fi

# step 6
if [[ $steps == *"6"* ]]; then
    echo
    echo "===================="
    echo "====== STEP 6 ======"
    echo "===================="
    echo
    python ./limvar_TEM_budget.py $ens 0 $mass $t1 $t2 $overwrite $dry
    #python ./limvar_TEM_budget.py $ens 1 $mass $t1 $t2 $overwrite $dry
fi

# step 7
if [[ $steps == *"7"* ]]; then
    echo
    echo "===================="
    echo "====== STEP 7 ======"
    echo "===================="
    echo
    python ./limvar_daily_to_10daily.py $ens $mass $t1 $t2 $overwrite $dry
fi
