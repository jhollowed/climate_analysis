#!/bin/bash

CASES="/glade/scratch/jhollowed/CAM/cases/aoa_runs/project3"

for d in ${CASES}/*/; do
    if [[ $d =~ "FV3" ]]; then
        printf '\nRemapping %s \n' "$d"
        if [[ $d =~ "C24" ]]; then
            srcgrid=C24
            dstgrid=4x4
        fi
        if [[ $d =~ "C48" ]]; then
            scrgrid=C48
            dstgrid=2x2
        fi
        griddir=/glade/u/home/cjablono/ncl/regrid/
        #file=$(ls $d/run/*h0* | head -n 1)

        for file in $d/run/*h[0-1]*.nc; do
            srcinitfile=${file%.*}.nc
            dstinitfile=${file%.*}.regrid.${dstgrid}.nc
            
            printf '\nremapping file %s to %s\n' "$srcinitfile" "$dstinitfile" 
            ncremap -m ${griddir}/map_${srcgrid}_to_${dstgrid}_aave.nc -i ${srcinitfile} -o ${dstinitfile}
            printf "\n\n"
        done
    fi
done

