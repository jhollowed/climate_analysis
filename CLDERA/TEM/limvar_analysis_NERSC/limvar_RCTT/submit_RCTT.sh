#!/bin/bash

ens=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15)
mass=(0 10)

for massi in "${mass[@]}"; do
    for ensi in "${ens[@]}"; do 
        echo "submitting mass=${massi}, ens=${ensi}"
        sbatch ./run_RCTT.sh $ensi $massi
    done
done
