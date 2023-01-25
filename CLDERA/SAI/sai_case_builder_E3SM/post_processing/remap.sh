#!/bin/bash

source /global/homes/j/jhollo/util/nco_utils.sh

##############################################################
# remapping
###############################################################

# for making map file if not exist...
#ncremap -s ne16pg2_scrip_c20191218.nc -g 91x180_SCRIP.20170401.nc -m map_ne16pg2_to_91x180.nc

map_dir=/global/cscratch1/sd/jhollo/E3SM/maps
regrid_str='91x180'
data_dir='ens_stats'

#for alg in 'bilinear' 'aave'; do
for alg in 'bilinear'; do

    if [[ $1 == *"np4"* ]]; then
        map=${map_dir}/map_ne16np4_to_91x180_${alg}.nc
    else
        map=${map_dir}/map_ne16pg2_to_91x180_${alg}.nc
    fi
   
    for input_file in ${data_dir}/$1.nc; do
        #echo "native: min(T) = `ncmin T ${input_file}`"
        remapped_file=${input_file%.nc}.regrid.${regrid_str}_${alg}.nc
        echo "output dest = " ${remapped_file}
        ncremap -m ${map} -i ${input_file} -o ${remapped_file}
        #echo "${alg}: min(T) = `ncmin T ${remapped_file}`"
    done
done
