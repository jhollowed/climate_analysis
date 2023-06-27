#!/bin/bash

##############################################################
# remapping
###############################################################

alg="bilinear"

#map="map_30pg2_to_180x360_bilinear.nc"
#input_file="fullvar_2ndHistorical_1850-2014/v2.LR.WCYCL20TR.0211.trc.pmcpu.eam.h0.2014-12.nc"

map="map_30np4_to_180x360_bilinear.nc"
input_file="fullvar_2ndHistorical_1850-2014/yearly_ic_avg.nc"

remapped_file=${input_file%.nc}.regrid.${regrid_str}_${alg}.nc
ncremap -m ${map} -i ${input_file} -o ${remapped_file}
