#!/bin/bash

#Nens  = int(sys.argv[1])
#tmini = int(sys.argv[2])
#tmaxi = int(sys.argv[3])
#cfb   = bool(int(sys.argv[4]))
#dry   = bool(int(sys.argv[5]))

cf=$1
ens=$2
# ensemble number is assigned based on the skybridge login node number if not passed
if [ -z "$2" ] 
then
    ens=$(hostname | tr -d -c 0-9)
fi

tmini=0
tmaxi=36
dry=0

python ./CLDERA_limvar_zinterp.py $ens $tmini $tmaxi $cf $dry

