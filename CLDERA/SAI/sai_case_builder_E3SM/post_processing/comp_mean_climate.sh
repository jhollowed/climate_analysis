#!/bin/bash

# obtain a single averaged "mean climate" by averaging over yearly outputs f# obtain a single averaged "mean climate" by averaging over yearly outputs from a 10-year HSW runrom a 10-year HSW run

DAT="./*mc/run/HSW_SAI_ne16pg2_L72_3600day_mc.eam.h0.0001-01-01-00000.nc"
ncra -d time,0,11 $DAT mean_climate.nc -D 7
