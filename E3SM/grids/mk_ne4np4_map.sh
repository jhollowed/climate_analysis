#/bin/bash
set -e
grid_root='/global/cscratch1/sd/jhollo/E3SM/grids'
map_root='/global/cscratch1/sd/jhollo/E3SM/maps'
source_grid_file='/global/cscratch1/sd/jhollo/E3SM/grids/ne4np4-pentagons_c100308.nc'
target_grid_file=${grid_root}/25x48_SCRIP.20170401.nc

source_grid='ne4np4'
target_grid='25x48'

# Create maps
for alg in 'bilinear' 'aave'; do
    map_file=${map_root}/map_${source_grid}_to_${target_grid}_${alg}.nc
    ncremap -a ${alg} -s ${source_grid_file} -g ${target_grid_file} -m ${map_file}
done
