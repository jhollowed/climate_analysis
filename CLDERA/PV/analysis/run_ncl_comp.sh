#!/bin/bash

ncl 'res_name="ne16np4"' 'compset_name="FIDEAL"' 'grid_name="dyn"' 'remap_res="91x180"' cldera_PV_cdf_or_quadratic.ncl
ncl 'res_name="ne16np4"' 'compset_name="FIDEAL"' 'grid_name="phys"' 'remap_res="91x180"' cldera_PV_cdf_or_quadratic.ncl

ncl 'res_name="ne16pg2"' 'compset_name="FIDEAL"' 'grid_name="dyn"' 'remap_res="91x180"' cldera_PV_cdf_or_quadratic.ncl
ncl 'res_name="ne16pg2"' 'compset_name="FIDEAL"' 'grid_name="phys"' 'remap_res="91x180"' cldera_PV_cdf_or_quadratic.ncl

#ncl 'res_name="ne4pg2"' 'compset_name="F1850"' 'grid_name="dyn"' 'remap_res="25x48"' cldera_PV_cdf_or_quadratic.ncl
#ncl 'res_name="ne4pg2"' 'compset_name="F1850"' 'grid_name="phys"' 'remap_res="25x48"' cldera_PV_cdf_or_quadratic.ncl
