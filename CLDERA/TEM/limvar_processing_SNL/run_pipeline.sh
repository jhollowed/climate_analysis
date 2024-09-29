#!/bin/bash

mass=$1 #10 or 0
#ens=$2
#./run_pipeline_for_member_and_massmag.sh $ens $mass 1234567 0 0
#./run_pipeline_for_member_and_massmag.sh $ens $mass 6 1 0

for ens in $(seq 1 13); do
    ./run_pipeline_for_member_and_massmag.sh $ens $mass 67 1 0
done
for ens in $(seq 91 95); do
    ./run_pipeline_for_member_and_massmag.sh $ens $mass 67 1 0
done
