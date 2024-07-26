#!/bin/bash

map_dir=/pscratch/sd/j/jhollo/E3SM/maps
data_dir=/global/homes/j/jhollo/repos/PyTEMDiags/PyTEMDiags/tests/data
#input_file=${data_dir}/f2.nc
input_file=${data_dir}/Y20.nc
regrid_str='180x360'
data_dir='.'
histnum=1

alg='aave'
map=${map_dir}/map_ne30pg2_to_180x360_${alg}.nc
remapped_file=${input_file%.nc}.regrid.${regrid_str}_${alg}.nc
echo "output dest = " ${remapped_file}
ncremap -m ${map} -i ${input_file} -o ${remapped_file} -a ${alg}

alg='bilinear'
map=${map_dir}/map_ne30pg2_to_180x360_${alg}.nc
remapped_file=${input_file%.nc}.regrid.${regrid_str}_${alg}.nc
echo "output dest = " ${remapped_file}
ncremap -m ${map} -i ${input_file} -o ${remapped_file} -a ${alg}

