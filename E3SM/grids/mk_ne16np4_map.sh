#/bin/bash
set -e
grid_root='/global/cscratch1/sd/jhollo/E3SM/grids'
map_root='/global/cscratch1/sd/jhollo/E3SM/maps'
source_grid_file='/global/cscratch1/sd/jhollo/E3SM/grids/ne16np4_110512_pentagons.nc'
target_grid_file=${grid_root}/91x180_SCRIP.20170401.nc

source_grid=ne16np4
target_grid='91x180'

# Create maps
for alg in 'bilinear' 'aave'; do
    map_file=${map_root}/map_${source_grid}_to_${target_grid}_${alg}.nc
    ncremap -a ${alg} -s ${source_grid_file} -g ${target_grid_file} -m ${map_file}
done
