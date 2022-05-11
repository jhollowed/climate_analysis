#!/bin/bash
set -e

CASES="/glade/scratch/jhollowed/CAM/cases_589/project4"

#for d in ${CASES}/*mod${1}__ens[0-9]${2}/; do
for d in ${CASES}/*mod[0-9]/; do
    if [[ $d =~ "FV3" ]]; then
        printf '\nRemapping %s \n' "$d"
        if [[ $d =~ "C24" ]]; then
            printf 'grid is C24'
            srcgrid=C24
            dstgrid=4x4
        fi
        if [[ $d =~ "C48" ]]; then
            printf 'grid is C48'
            srcgrid=C48
            dstgrid=2x2
        fi
        griddir=/glade/u/home/cjablono/ncl/regrid
        #file=$(ls $d/run/*h0* | head -n 1)

        for file in $d/run/*h[0-1]*000.nc; do
            srcinitfile=${file%.*}.nc
            dstinitfile=${file%.*}.regrid.${dstgrid}.nc
            
            printf '\nremapping file %s to %s\n' "$srcinitfile" "$dstinitfile" 
            ncremap -m ${griddir}/map_${srcgrid}_to_${dstgrid}_aave.nc -i ${srcinitfile} -o ${dstinitfile}
            printf "\n\n"
        done
    fi
done
exit 0

# ======= BKP


CASES="/glade/scratch/jhollowed/CAM/cases_589/project4"

for d in ${CASES}/*/; do
    if [[ $d =~ "FV3" ]]; then
        printf '\nRemapping %s \n' "$d"
        if [[ $d =~ "C24" ]]; then
            printf 'grid is C24'
            srcgrid=C24
            dstgrid=4x4
        fi
        if [[ $d =~ "C48" ]]; then
            printf 'grid is C48'
            srcgrid=C48
            dstgrid=2x2
        fi
        griddir=/glade/u/home/cjablono/ncl/regrid
        #file=$(ls $d/run/*h0* | head -n 1)

        for file in $d/run/*h[0-1]*000.nc; do
            srcinitfile=${file%.*}.nc
            dstinitfile=${file%.*}.regrid.${dstgrid}.nc
            
            printf '\nremapping file %s to %s\n' "$srcinitfile" "$dstinitfile" 
            ncremap -m ${griddir}/map_${srcgrid}_to_${dstgrid}_aave.nc -i ${srcinitfile} -o ${dstinitfile}
            printf "\n\n"
        done
    fi
done
