#!/bin/bash

var=$1
data_in_loc="/pscratch/sd/w/wagmanbe/E3SM_simulations/CLDERA/v2.LR.WCYCL20TR.0211.trc.pmcpu/archive/atm/hist"
data_out_loc="/pscratch/sd/j/jhollo/E3SM/historical/CLDERA_2ndHistorical_1950-2014_monthly"

output_file="${data_out_loc}/histoircal_h0_${var}_1850-2014.nc"
echo "working on $var"
ncrcat -D 2 -O -v $var ${data_in_loc}/*h0*.nc $output_file
