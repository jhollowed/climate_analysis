#/bin/bash
set -e
grid_root=.
map_root=.
source_grid_file=${grid_root}/ne30pg2_scrip_20200209.nc
#source_grid_file=${grid_root}/ne30pg2_scrip.nc
target_grid_file=${grid_root}/cmip6_180x360_scrip.20181001.nc

alg='bilinear'
res=30
source_grid=${res}pg2
target_grid='180x360'

# Create source grid file
if [ ! -e ${source_grid_file} ]; then
    GenerateCSMesh --alt --res ${res} --file ${grid_root}/ne${res}.g
    GenerateVolumetricMesh --in ${grid_root}/ne${res}.g --out ${grid_root}/ne${res}pg2.g --np 2 --uniform
    ConvertExodusToSCRIP --in ${grid_root}/ne${res}pg2.g --out ${grid_root}/ne${res}pg2.scrip.nc
    source_grid_file=${grid_root}/ne${res}pg2.scrip.nc
fi

# Create maps
map_file=${map_root}/map_${source_grid}_to_${target_grid}_${alg}.nc
ncremap -a ${alg} -s ${source_grid_file} -g ${target_grid_file} -m ${map_file}
