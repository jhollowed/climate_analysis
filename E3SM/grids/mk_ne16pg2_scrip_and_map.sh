#/bin/bash
set -e
grid_root='/global/cscratch1/sd/jhollo/E3SM/grids'
map_root='/global/cscratch1/sd/jhollo/E3SM/maps'
source_grid_file=${grid_root}/ne16pg2.scrip.nc
target_grid_file=${grid_root}/91x180_SCRIP.20170401.nc

res=16
source_grid=${res}pg2
target_grid='91x180'

# Create source grid file
if [ ! -e ${source_grid_file} ]; then
    GenerateCSMesh --alt --res ${res} --file ${grid_root}/ne${res}.g
    GenerateVolumetricMesh --in ${grid_root}/ne${res}.g --out ${grid_root}/ne${res}pg2.g --np 2 --uniform
    ConvertExodusToSCRIP --in ${grid_root}/ne${res}pg2.g --out ${grid_root}/ne${res}pg2.scrip.nc
fi

# Create maps
for alg in 'bilinear' 'aave'; do
    map_file=${map_root}/map_${source_grid}_to_${target_grid}_${alg}.nc
    ncremap -a ${alg} -s ${source_grid_file} -g ${target_grid_file} -m ${map_file}
done
