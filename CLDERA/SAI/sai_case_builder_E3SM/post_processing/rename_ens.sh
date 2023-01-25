#!/bin/bash

for i in {1..5}
do
    if [ ${i} -eq 5 ]; then PFX=''; else PFX='0'; fi
    ENS_NUM=$(($i+5))
    ENS_NUM=${PFX}${ENS_NUM}
    echo ens${ENS_NUM}

    # rename matches in files
    # https://superuser.com/questions/428493/how-can-i-do-a-recursive-find-and-replace-from-the-command-line
    find ./*ens${ENS_NUM} -type f \( -iname \*.nml -o -iname \*.txt -o -iname \*.atm -o -iname \*.json -o -iname \*.drv -o -iname CASEROOT -o -iname drv_in \) -print0 | xargs -0 sed -i'' -e "s/__ens${ENS_NUM}/_ens0${i}/g"
    find ./*ens${ENS_NUM} -type f \( -iname \*.nml -o -iname \*.txt -o -iname \*.atm -o -iname \*.json -o -iname \*.drv -o -iname CASEROOT -o -iname drv_in \) -print0 | xargs -0 sed -i'' -e "s/_ens${ENS_NUM}/_ens0${i}/g"

    # rename file names
    mv ./HSW_SAI_ne16pg2_L72_900day_180delay__ens${ENS_NUM} HSW_SAI_ne16pg2_L72_900day_180delay_ens0${i}
    for file in ./*ens0${i}/run/* ; do mv $file ${file//__ens${ENS_NUM}/_ens0${i}} ; done
    for file in ./*ens0${i}/run/* ; do mv $file ${file//_ens${ENS_NUM}/_ens0${i}} ; done
    
    # change netcdf attributes
    for file in ./*ens0${i}/run/*eam.h[0-2].*.nc ; do
        ncatted -O -a case,global,m,c,"HSW_SAI_ne16pg2_L72_900day_180delay__ens0${i}" $file
    done
    for file in ./*ens0${i}/run/*eam.i.*.nc ; do
        ncatted -O -a case,global,m,c,"HSW_SAI_ne16pg2_L72_900day_180delay__ens0${i}" $file
    done
    for file in ./*ens0${i}/run/*eam.r.*.nc ; do
        ncatted -O -a caseid,global,m,c,"HSW_SAI_ne16pg2_L72_900day_180delay__ens0${i}" $file
    done
done
