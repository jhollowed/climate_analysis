#!/bin/bash

cp held_suarez_1994.F90 aman_mods/src.cam
cp aman.physpkg.F90 aman_mods/src.cam/physpkg.F90
cp aman.tracers.F90 aman_mods/src.cam/tracers.F90
cp aoa.aoa_tracers.F90 aman_mods/src.cam/aoa_tracers.F90

cp held_suarez_1994.F90 aoa_mods/src.cam
cp aoa.physpkg.F90 aoa_mods/src.cam/physpkg.F90
cp aoa.aoa_tracers.F90 aoa_mods/src.cam/aoa_tracers.F90
