#!/bin/bash

mass=$1 #10 or 0
#ens=$2
#./run_pipeline_for_member_and_massmag.sh $ens $mass 1234567 0 0
#./run_pipeline_for_member_and_massmag.sh $ens $mass 6 1 0


batch=$2

if [ "$batch" -eq 1 ]; then

    for ens in $(seq 1 9); do
        ./run_pipeline_for_member_and_massmag.sh $ens $mass 67 1 0
    done
fi

if [ "$batch" -eq 2 ]; then

    for ens in $(seq 10 13); do
        ./run_pipeline_for_member_and_massmag.sh $ens $mass 67 1 0
    done
    for ens in $(seq 91 95); do
        ./run_pipeline_for_member_and_massmag.sh $ens $mass 67 1 0
    done

fi
