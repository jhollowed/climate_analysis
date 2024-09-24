#!/bin/bash

#Nens  = int(sys.argv[1])
#cfb   = bool(int(sys.argv[2]))
#dry   = bool(int(sys.argv[3]))

cf=$1
ens=$2
# ensemble number is assigned based on the skybridge login node number if not passed
if [ -z "$2" ] 
then
    ens=$(hostname | tr -d -c 0-9)
fi

dry=0

python ./CLDERA_limvar_zonalMeans_monthly_timeConcat.py $ens $cf $dry

