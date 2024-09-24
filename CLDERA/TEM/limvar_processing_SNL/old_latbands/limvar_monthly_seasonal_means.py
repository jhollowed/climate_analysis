'''
Joe Hollowed 
University of Michigan 2024

This script computes monthly and seasonal means of limvar TEM data, and limvar zonal means. 
Some limvar zonal means will already exist from the h0 history files. These variables will 
be combined with variables which only exist as daily data, monthly averaged.
Command line arguments are:

Usage
-----
python ./limvar_zinterp.py [Nens] [massMag] [tmini] [tmaxi] [overwrite] [dry]

Parameters
----------
Nens : int
    the ensemble member specified as an integer
massMag : int
    the mass magnitude of the ensemble. Valid values are:
    0, 1, 3, 5, 7, 10, 13, 15
    if massMag = 0, this sepcifies the counterfactual ensemble
tmini : int
    index of the first history file to include in the averaging
tmaxi : int
    index of the last history file to include in the averaging
overwrite : int
    whether or not to overwrite the averaged data if it already exists
    at the specified outpit file. 0 = False, 1 = True, 
    If False, then the script skips the averaging for any file that 
    already exists
dry : int
    whether or not to do a "dry run" of the function, printing information
    about the execution but skipping the averaging. 0 = False, 1 = True

Note that tmini, tmaxi specify only the *index* of the identified time files. 
Each file may not cover the same amount of model time. 
'''

import sys
import pdb
import glob
import dask
import copy
import cftime
import os.path
import numpy as np
import xarray as xr
from datetime import datetime

monthlyloc = '/gpfs/cldera/data/e3sm/e3sm-cldera/CLDERA-historical/limvar_src_tag/tem_processed/interp_to_pressure_zonalMeans'
dailyloc = '/ascldap/users/jphollo/data/limvar/limvar_daily'
monthlyoutloc   = '/ascldap/users/jphollo/data/limvar/limvar_monthly' 
seasonaloutloc   = '/ascldap/users/jphollo/data/limvar/limvar_seasonal' 

Nens      = int(sys.argv[1])
massMag   = int(sys.argv[2])
tmini     = int(sys.argv[3])
tmaxi     = int(sys.argv[4])
overwrite = bool(int(sys.argv[5]))
dry       = bool(int(sys.argv[6]))
print('args: {}'.format(sys.argv))

cfb = massMag == 0   # counterfactual flag


# ----------------------------------------------------------------

dpm = {'noleap': [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       '365_day': [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       'standard': [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       'gregorian': [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       'proleptic_gregorian': [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       'all_leap': [0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       '366_day': [0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       '360_day': [0, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30]}


def leap_year(year, calendar='standard'):
    '''
    Determine if year is a leap year
    '''
    leap = False
    if ((calendar in ['standard', 'gregorian',
        'proleptic_gregorian', 'julian']) and
        (year % 4 == 0)):
        leap = True
        if ((calendar == 'proleptic_gregorian') and
            (year % 100 == 0) and
            (year % 400 != 0)):
            leap = False
        elif ((calendar in ['standard', 'gregorian']) and
                 (year % 100 == 0) and (year % 400 != 0) and
                 (year < 1583)):
            leap = False
    return leap

def get_dpm(time, calendar='standard'):
    '''
    return a array of days per month corresponding to the months provided in `months`
    '''
    month_length = np.zeros(len(time), dtype=int)

    cal_days = dpm[calendar]

    for i, (month, year) in enumerate(zip(time.month, time.year)):
        month_length[i] = cal_days[month]
        if leap_year(year, calendar=calendar) and month == 2:
            month_length[i] += 1
    return month_length

def change_djf_year(ds):
    '''
    Changes the year of DJF seasons to be the yeart o that January, rather than that December
    '''
    time = ds.time.values
    for i in range(len(time)):
        if(time[i].month == 12):
            time[i] = cftime.DatetimeNoLeap(time[i].year+1, time[i].month, time[i].day, 
                                            time[i].hour, time[i].minute, time[i].second)
    ds['time'] = time
    return ds

def season_mean(ds, calendar='standard'):
    '''
    Takes seasonal means as 3-month intervals starting each December
    '''
    
    month_length = xr.DataArray(
        get_dpm(ds.time.to_index(), calendar=calendar),
        coords=[ds.time],
        name='month_length'
    )
    
    # compute monthly mean 
    result = ((ds.drop_vars('time_bnds') * month_length).resample(time='QS-DEC').sum() / 
              month_length.resample(time='QS-DEC').sum())
    result = change_djf_year(result)

    # build new time bounds
    grouped_time_bnds = ds['time_bnds'].resample(time='QS-DEC')
    new_time_bnds = []
    for i in range(len(result.time)):
        tb = ds['time_bnds'].isel(time=list(grouped_time_bnds.groups.values())[i])
        new_time_bnds.append([tb.values[0][0], tb.values[-1][-1]])

    new_time_bnds = xr.DataArray(new_time_bnds, dims=('time', 'nbnd'))
    new_time_bnds.assign_coords(time=result.time)
    result['time_bnds'] = new_time_bnds

    return result


# ----------------------------------------------------------------


# get model monthly mean time bounds
if(not cfb):
    zmh0 = sorted(glob.glob('{}/*{}Tg*ens{}.eam.h0*zonalmeans.nc'.format(
                                monthlyloc, massMag, Nens)))
else:
    zmh0 = sorted(glob.glob('{}/*ens{}.cf.eam.h0*zonalmeans.nc'.format(
                                monthlyloc, Nens)))
zmh0 = xr.concat([xr.open_dataset(ds) for ds in zmh0], dim='time')
h0_time_bounds = zmh0['time_bnds']

def monthly_mean(ds, dsh0_time_bnds=h0_time_bounds):
    '''
    Take the monthly mean of daily data

    Parameters
    ----------
    ds : xarray DataArray
        the data to average
    dsh0 : xarray DataArray
        monthly-mean data time bounds as output from the model
    '''
    # previously, we were taking monthly averages simply as 
    # ds.resample(time='1MS').mean('time')
    # However, the result of this calculation does not match the monthly mean data
    # as output from the model. It looks like the date range used in the model is different
    # (shifted by one day, running from the second of month A, to the first of month B,
    #  inclusive). So, we need to do this manually.
    
    # first compute monthly mean the old way as a template
    ds_monmean = ds.resample(time='1MS').mean('time')
    ds_monmean = ds_monmean.transpose('time', 'lat', 'plev')
    
    # drop vars that are not on (lat, plev, time)
    if('P0' in ds.data_vars): 
        P0 = ds['P0']
        ds = ds.drop_vars('P0')
    if('time_bnds' in ds.data_vars): 
        time_bnds = ds['time_bnds']
        ds = ds.drop_vars('time_bnds')

    if(len(dsh0_time_bnds) != len(ds_monmean.time)):
        raise RuntimeError('input data and model monthly data don\'t cover same amunt of time')

    for i in range(len(dsh0_time_bnds)):
        start, end = dsh0_time_bnds[i].values[0], dsh0_time_bnds[i].values[1]
        delta =  (end - start).days
        x = ds.sel(time=slice(start, end)).isel(time=slice(1, delta+1)).mean('time')
        for var in x.data_vars:
            ds_monmean[var][i,:,:] = x[var]
    
    dsh0_time_bnds['time']  = ds_monmean.time
    ds_monmean['time_bnds'] = dsh0_time_bnds
    if('P0' in ds.data_vars): ds_monmean['P0'] = P0
    return ds_monmean
    

# ----------------------------------------------------------------
    
skip_tem = 0 # flag to skip averaging the TEM data only averaging the zonal-mean model variables

if(not skip_tem):
    # ---- First average TEM data
    print('======== locating TEM data for ens{}, {} Tg...'.format(Nens, massMag))
    # loop over tracers
    for qi in [0, 1, 2]:
        qstr  = ['','_TRACER-AOA','_TRACER-E90j'][int(qi)]
        print('====== working on tracer {}...'.format(qi))
        if(not cfb):
            enstem = sorted(glob.glob('{}/*{}Tg*ens{}.eam.h1*TEM*L45{}.nc'.format(
                                        dailyloc, massMag, Nens, qstr)))[0]
        else:
            enstem = sorted(glob.glob('{}/*ens{}.cf.eam.h1*TEM*L45{}.nc'.format(
                                        dailyloc, Nens, qstr)))[0]
       
        print('time concatenated file: {}'.format(enstem))
        name = enstem.split('/')[-1].split('.nc')[0]
        name = name.split('199')[0] + name.split('000_')[-1]

        
        print('---- computing monthly mean...')
        if(not dry): 
            
            print('reading data...')
            tem = xr.open_dataset(enstem)
            
            # skip if exists
            outfile = '{}/{}_monthlymean.nc'.format(monthlyoutloc, name)
            existing = glob.glob(outfile)
            if(len(existing) > 0 and not overwrite): 
                print('monthly file exists; skipping...')
            else:
                print('averaging data...')
                tem = monthly_mean(tem)
                print('writing out...')
                tem.to_netcdf(outfile)
        
            print('---- computing seasonal mean...')
            
            # skip if exists
            outfile = '{}/{}_seasonalmean.nc'.format(seasonaloutloc, name)
            existing = glob.glob(outfile)
            if(len(existing) > 0 and not overwrite):
                print('seasonal file exists; skipping...')
            else:
                print('averaging data...')
                tem = season_mean(tem)
                print('writing out...')
                tem.to_netcdf(outfile)


# ---- Now do the same for the zonal mean data
print('======== locating zonal mean data for ens{}, {} Tg...'.format(Nens, massMag))

if(not cfb):
    enszm   = sorted(glob.glob('{}/*{}Tg*ens{}.eam.h1*zonalmeans.nc'.format(
                                dailyloc, massMag, Nens)))[0]
else:
    enszm   = sorted(glob.glob('{}/*ens{}.cf.eam.h1*zonalmeans.nc'.format(
                                dailyloc, Nens)))[0]

print('time concatenated file: {}'.format(enszm))
name = enszm.split('/')[-1].split('.nc')[0]
name = name.split('199')[0] + name.split('000_')[-1]


print('---- computing monthly mean...')
if(not dry): 
    
    print('reading data...')
    zm   = xr.open_dataset(enszm)
    
    # skip if exists
    outfile = '{}/{}_monthlymean.nc'.format(monthlyoutloc, name)
    existing = glob.glob(outfile)
    if(len(existing) > 0 and not overwrite): 
        print('monthly file exists; skipping...')
    else:
        print('averaging data...')
        zm = monthly_mean(zm)
        
        print('combining with existing monthly data...')
        # do sanity check on U...
        min_err = (zm.isel(time=0) - zmh0.isel(time=0)).min()
        mean_err = (zm.isel(time=0) - zmh0.isel(time=0)).mean()
        max_err = (zm.isel(time=0) - zmh0.isel(time=0)).max()
        print('min, mean, max U err: {}, {}, {}'.format(np.atleast_1d(min_err['U'].values)[0], 
                                                        np.atleast_1d(mean_err['U'].values)[0], 
                                                        np.atleast_1d(max_err['U'].values)[0]))

        # add variables from monthly data
        add_vars = [d for d in zmh0.data_vars if d not in zm.data_vars]
        for var in add_vars:
            add_var = zmh0[var]
            add_var['time'] = zm.time
            zm[var] = add_var

        print('writing out...')
        zm.to_netcdf(outfile)

    print('---- computing seasonal mean...')
    
    # skip if exists
    outfile = '{}/{}_seasonalmean.nc'.format(seasonaloutloc, name)
    existing = glob.glob(outfile)
    if(len(existing) > 0 and not overwrite):
        print('seasonal file exists; skipping...')
    else:
        print('averaging data...')
        zm = season_mean(zm)
        print('writing out...')
        zm.to_netcdf(outfile)
