#/bin/bash
set -e
grid_root=.
map_root=.
source_grid_file=${grid_root}/ne30np4_pentagons.091226.nc
target_grid_file=${grid_root}/cmip6_180x360_scrip.20181001.nc

alg='bilinear'
res=30
source_grid=${res}np4
target_grid='180x360'

# Create maps
map_file=${map_root}/map_${source_grid}_to_${target_grid}_${alg}.nc
ncremap -a ${alg} -s ${source_grid_file} -g ${target_grid_file} -m ${map_file}
