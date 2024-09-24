#!/bin/bash

mass=$1
t1=0
t2=36
overwrite=1

# this should be a string contining the steps to perform, e.g.
# 12345678910
# or
# 34567
steps=$2


# step 1
if [[ $steps == *"1[^0]"* ]]; then
    echo
    echo "===================="
    echo "====== STEP 1 ======"
    echo "===================="
    echo
    for ens in $(seq 1 10); do
        python ./limvar_zinterp.py $ens 0 $mass $t1 $t2 $overwrite 0 # monthly
        python ./limvar_zinterp.py $ens 1 $mass $t1 $t2 $overwrite 0 # daily
    done
fi

# step 2
if [[ $steps == *"2"* ]]; then
    echo
    echo "===================="
    echo "====== STEP 2 ======"
    echo "===================="
    echo
    for ens in $(seq 1 10); do
        python ./limvar_TEM.py $ens $mass $t1 $t2 $overwrite 0
    done
fi

# step 3
if [[ $steps == *"3"* ]]; then
    echo
    echo "===================="
    echo "====== STEP 3 ======"
    echo "===================="
    echo
    for ens in $(seq 1 10); do
        #python ./limvar_zonalmeans.py $ens 0 $mass $t1 $t2 $overwrite 0 0 # monthly
        #python ./limvar_zonalmeans.py $ens 1 $mass $t1 $t2 $overwrite 0 0 # daily
        python ./limvar_zonalmeans.py $ens 0 $mass $t1 $t2 0 0 1 # monthly
        python ./limvar_zonalmeans.py $ens 1 $mass $t1 $t2 0 0 1 # daily
    done
fi

#step 4
if [[ $steps == *"4"* ]]; then
    echo
    echo "===================="
    echo "====== STEP 4 ======"
    echo "===================="
    echo
    for ens in $(seq 1 10); do
        python ./limvar_timeConcat.py $ens 0 $mass 0 $t1 $t2 $overwrite 0 # monthly
        python ./limvar_timeConcat.py $ens 1 $mass $t1 $t2 $overwrite 0      # daily
    done
fi

#step 5
if [[ $steps == *"5"* ]]; then
    echo
    echo "===================="
    echo "====== STEP 5 ======"
    echo "===================="
    echo
    for ens in $(seq 1 10); do
        python ./limvar_monthly_seasonal_means.py $ens $mass $t1 $t2 $overwrite 0
    done
fi

# step 6
if [[ $steps == *"6"* ]]; then
    echo
    echo "===================="
    echo "====== STEP 6 ======"
    echo "===================="
    echo
    for ens in $(seq 1 10); do
        python ./limvar_latband_means.py $ens $mass $t1 $t2 $overwrite 0
    done
fi

# step 7
if [[ $steps == *"7"* ]]; then
    echo
    echo "===================="
    echo "====== STEP 7 ======"
    echo "===================="
    echo
    for ens in $(seq 1 10); do
        python ./limvar_latband_10daily.py $ens $mass $t1 $t2 $overwrite 0
    done
fi

# step 8
if [[ $steps == *"8"* ]]; then
    echo
    echo "===================="
    echo "====== STEP 8 ======"
    echo "===================="
    echo
    for ens in $(seq 1 10); do
        python ./limvar_TEM_budget_monthly.py $ens $mass $t1 $t2 $overwrite 0
    done
fi

# step 9
if [[ $steps == *"9"* ]]; then
    echo
    echo "===================="
    echo "====== STEP 9 ======"
    echo "===================="
    echo
    for ens in $(seq 1 10); do
        python ./limvar_TEM_budget_latbands.py $ens $mass $t1 $t2 $overwrite 0
    done
fi

# step 10
if [[ $steps == *"10"* ]]; then
    echo
    echo "===================="
    echo "====== STEP 10 ======"
    echo "===================="
    echo
    for ens in $(seq 1 10); do
        python ./limvar_TEM_budget_latbands_10daily.py 1 $mass $t1 $t2 $overwrite 0
    done
fi
