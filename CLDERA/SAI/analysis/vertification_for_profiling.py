# Joe Hollowed
# University of Michigan 2023
# 
# Performing comparisons between one of my released enesemble members, and independent runs
# performed by Andrew Steyer. As of 5/4/23, there are concerns that my runs cannot be reproduced; 
# this script compares the datasets in a few cummary statistics

import pdb
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

AS_run = '/global/cscratch1/sd/jhollo/E3SM/E3SMv2_cases/sai_cases/verification_for_profiling/E3SM_ne16pg2_ne16pg2_L72_FIDEAL_SAI.eam.h0.0001-01-01-00000.nc'
JH_run = '/project/projectdirs/m4014/data/HSW/outputs/release_030123/netcdf/low_var_ens/ens05/HSW_SAI_ne16pg2_L72_1200day_90delay__ens05.eam.h0.0001-01-01-00000.nc'
JW_run = '/project/projectdirs/m4014/data/HSW/outputs/release_030123/netcdf/low_var_mass_ens/ens05_mass1.00X/E3SM_ne16pg2_ne16pg2_L72_FIDEAL_SAI__ens05_mass1.00X.eam.h0.0001-01-01-00000.nc'
mean_climate = '/global/cfs/cdirs/m4014/data/HSW/outputs/release_011423/netcdf/ens_stats/mean_climate.nc'

print('opening datasets')
d1 = xr.open_dataset(AS_run)
d2 = xr.open_dataset(JH_run)
d3 = xr.open_dataset(JW_run)
lat, lon, lev, ncol = d1['lat'], d1['lon'], d1['lev'], d1['ncol']
ndcur, nscur = d1['ndcur'], d1['nscur']
time = ndcur + nscur/(3600*24)

xr.testing.assert_allclose(lat, d2['lat'])
xr.testing.assert_allclose(lon, d2['lon'])
xr.testing.assert_allclose(lev, d2['lev'])
xr.testing.assert_allclose(ncol, d2['ncol'])
xr.testing.assert_allclose(ndcur, d2['ndcur'])
xr.testing.assert_allclose(nscur, d2['nscur'])

xr.testing.assert_allclose(lat, d3['lat'])
xr.testing.assert_allclose(lon, d3['lon'])
xr.testing.assert_allclose(lev, d3['lev'])
xr.testing.assert_allclose(ncol, d3['ncol'])
xr.testing.assert_allclose(ndcur, d3['ndcur'])
xr.testing.assert_allclose(nscur, d3['nscur'])

check_andrews_run = 1
check_jerrys_run = 1

# ---- check time of injection
tmpdir = '/global/homes/j/jhollo/repos/climate_analysis/CLDERA/SAI/analysis/tmp'
if(1):
    
    print('checking injection time')
    m1 = (d1['SO2'] != 0).sum(dim=('ncol', 'lev'))
    m2 = (d2['SO2'] != 0).sum(dim=('ncol', 'lev'))    
    m3 = (d3['SO2'] != 0).sum(dim=('ncol', 'lev'))
    
    t01 = float(time[np.where(m1.values != 0)[0][0]])
    t02 = float(time[np.where(m2.values != 0)[0][0]])
    t03 = float(time[np.where(m3.values != 0)[0][0]])
    print('injection time for dataset d1: {}'.format(t01))
    print('injection time for dataset d2: {}'.format(t02))
    print('injection time for dataset d3: {}'.format(t03))

# ===============================================
# ============= ANDREWS RUN =====================

if(check_andrews_run):
    # compare global T mean at all levels
    tmp = '/global/homes/j/jhollo/repos/climate_analysis/CLDERA/SAI/analysis/tmp/Tdiff_andrew.nc'
    tmp_srf1 = '/global/homes/j/jhollo/repos/climate_analysis/CLDERA/SAI/analysis/tmp/T1000_joe.nc'
    tmp_srf2 = '/global/homes/j/jhollo/repos/climate_analysis/CLDERA/SAI/analysis/tmp/T1000_andrew.nc'
    try:
        Tdiff = xr.open_dataset(tmp)['T']
        T1000_joe = xr.open_dataset(tmp_srf1)['T1000']
        T1000_joe = xr.open_dataset(tmp_srf1)['T1000']
        print('read Tdiff')
    except FileNotFoundError:
        print('reading T1')
        T1_gm = d1['T'].mean('ncol')
        print('reading T2')
        T2_gm = d2['T'].mean('ncol')
        print('taking diff')
        Tdiff = T1_gm - T2_gm
        Tdiff.to_netcdf(tmp)
        print('wrote Tdiff')

    print('plotting Tdiff')
    fig = plt.figure()
    ax = fig.add_subplot(111)

    LEV, TIME = np.meshgrid(lev, time)
    cax_T = ax.contourf(TIME, LEV, Tdiff, levels=20, cmap=plt.cm.rainbow)
    cbar = fig.colorbar(cax_T)
    cbar.set_label('global-mean T difference\n((Andrew\'s run) - (low_var_ens ens05))')

    ax.set_ylim([5, 1000])
    ax.set_xlim([80, 180])
    ax.set_yscale('log')
    ax.invert_yaxis()
    ax.set_xlabel('time [days]')
    ax.set_ylabel('lev [hPa]')

    # compare global SO2 mean at all levels
    tmp = '/global/homes/j/jhollo/repos/climate_analysis/CLDERA/SAI/analysis/tmp/SO2diff_andrew.nc'
    try:
        SO2diff = xr.open_dataset(tmp)['SO2']
        print('read SO2diff')
    except FileNotFoundError:
        print('reading SO21')
        SO21_gm = d1['SO2'].mean('ncol')
        print('reading SO22')
        SO22_gm = d2['SO2'].mean('ncol')
        print('taking diff')
        SO2diff = SO21_gm - SO22_gm
        SO2diff.to_netcdf(tmp)
        print('wrote SO2diff')

    print('plotting SO2diff')
    fig = plt.figure()
    ax = fig.add_subplot(111)

    cax_SO2 = ax.contourf(TIME, LEV, SO2diff, levels=20, cmap=plt.cm.rainbow)
    cbar = fig.colorbar(cax_SO2)
    cbar.set_label('global-mean SO2 difference\n((Andrew\'s run) - (low_var_ens ens05))')

    ax.set_ylim([5, 1000])
    ax.set_xlim([80, 180])
    ax.set_yscale('log')
    ax.invert_yaxis()
    ax.set_xlabel('time [days]')
    ax.set_ylabel('lev [hPa]')


# ==================================================
# ============= JERRYS ENS RUN =====================

if(check_jerrys_run):
    # compare global T mean at all levels
    tmp = '/global/homes/j/jhollo/repos/climate_analysis/CLDERA/SAI/analysis/tmp/Tdiff_jerry.nc'
    try:
        Tdiff = xr.open_dataset(tmp)['T']
        print('read Tdiff')
    except FileNotFoundError:
        print('reading T1')
        T1_gm = d1['T'].mean('ncol')
        print('reading T2')
        T2_gm = d3['T'].mean('ncol')
        print('taking diff')
        Tdiff = T1_gm - T2_gm
        Tdiff.to_netcdf(tmp)
        print('wrote Tdiff')

    print('plotting Tdiff')
    fig = plt.figure()
    ax = fig.add_subplot(111)

    LEV, TIME = np.meshgrid(lev, time)
    cax = ax.contourf(TIME, LEV, Tdiff, levels=cax_T.levels, cmap=plt.cm.rainbow)
    cbar = fig.colorbar(cax)
    cbar.set_label('global-mean T difference\n((Jerry\'s ens05_mass1.00X) - (low_var_ens ens05))')

    ax.set_ylim([5, 1000])
    ax.set_xlim([80, 180])
    ax.set_yscale('log')
    ax.invert_yaxis()
    ax.set_xlabel('time [days]')
    ax.set_ylabel('lev [hPa]')

    # compare global SO2 mean at all levels
    tmp = '/global/homes/j/jhollo/repos/climate_analysis/CLDERA/SAI/analysis/tmp/SO2diff_jerry.nc'
    try:
        SO2diff = xr.open_dataset(tmp)['SO2']
        print('read SO2diff')
    except FileNotFoundError:
        print('reading SO21')
        SO21_gm = d1['SO2'].mean('ncol')
        print('reading SO22')
        SO22_gm = d3['SO2'].mean('ncol')
        print('taking diff')
        SO2diff = SO21_gm - SO22_gm
        SO2diff.to_netcdf(tmp)
        print('wrote SO2diff')

    print('plotting SO2diff')
    fig = plt.figure()
    ax = fig.add_subplot(111)

    cax = ax.contourf(TIME, LEV, SO2diff, levels=cax_SO2.levels, cmap=plt.cm.rainbow)
    cbar = fig.colorbar(cax)
    cbar.set_label('global-mean SO2 difference\n((Jerry\'s ens05_mass1.00X) - (low_var_ens ens05))')

    ax.set_ylim([5, 1000])
    ax.set_xlim([80, 180])
    ax.set_yscale('log')
    ax.invert_yaxis()
    ax.set_xlabel('time [days]')
    ax.set_ylabel('lev [hPa]')

#plt.show()
