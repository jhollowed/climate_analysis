#!/bin/bash

# compute ensemble mean
# $1 = hist num 0,1,or 2 
# for hist num, all output is on one file, so skip for loop over year, month

if [ "$1" -lt "2" ]; then
    for YEAR in 0001 0002 0003; do
        for MONTH in 01 06 12; do
            find . -name "*ens0*eam.h${1}.${YEAR}-${MONTH}*00.nc" -print0 | xargs -r0 ncea -o ens_stats/ens_mean.h${1}.${YEAR}-${MONTH}.nc -D 7
        done
    done
else
    find . -name "*ens0**eam.h${1}.*00.nc" -print0 | xargs -r0 ncea -o ens_stats/ens_mean.h${1}.nc -D 7
fi
