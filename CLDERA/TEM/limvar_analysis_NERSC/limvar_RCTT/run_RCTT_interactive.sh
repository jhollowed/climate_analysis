#!/bin/bash

ens=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15)
mass=$1

for ensi in "${ens[@]}"; do 
    echo "submitting mass=${massi}, ens=${ensi}"
    mpiexec -n 32 python RCTT_on_combined_limvar_fullvar.py $ensi $mass
done
