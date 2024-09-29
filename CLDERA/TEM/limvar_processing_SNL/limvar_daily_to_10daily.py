'''
Joe Hollowed 
University of Michigan 2024

This script computes 10-day averages of limvar zonal means, limvar TEM data, and limvar TEM budgets.

Command line arguments are:

Usage
-----
python ./limvar_zinterp.py [Nens] [massMag] [tmini] [tmaxi] [overwrite] [dry] [time_window_overlap]

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
time_window_overlap : int
    number of days to overlap the 10-day averaging windows by 
    time_window_overlap = 0 
        will produce consecutive windows with no overlap, and a 10x reduction in data size
    time_window_overlap = 5
        will produce 10-day means ever 5 days, with 5-days on each side of the average being
        shared by the neighboring averages, and a 5x reduction in datasize.

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

#dailyloc = '/ascldap/users/jphollo/data/limvar/limvar_daily'
dailyloc = '/gpfs/cldera/data/e3sm/e3sm-cldera/CLDERA-historical/limvar_src_tag/tem_processed/timeConcat_zonalMeans'
outloc   = '/ascldap/users/jphollo/data/limvar/limvar_10daily' 

print('args: {}'.format(sys.argv))
Nens      = int(sys.argv[1])
massMag   = int(sys.argv[2])
tmini     = int(sys.argv[3])
tmaxi     = int(sys.argv[4])
overwrite = bool(int(sys.argv[5]))
dry       = bool(int(sys.argv[6]))
time_window_overlap = 5

cfb = massMag == 0   # counterfactual flag
ovr = time_window_overlap # shorthand

# if Nens was >90, this will denote the non-source-tagged ensemble
if(Nens < 90):
    no_src_tag = False
    massStr = '.{}Tg.'.format(massMag)
else:
    massStr = ''
    no_src_tag = True
# the non-source-tagged ensemble only exists for 10 Tg and 0 Tg
if(massMag != 0 and massMag != 10 and no_src_tag):
    raise RuntimeError('must have massMag = 10 or 0 for the non-source tagged ensemble')

#  CHECK THESE FLAGS BEFORE RUNNING
SKIP_TEM    = True
SKIP_ZM     = True
SKIP_BUDGET = False

# ----------------------------------------------------------------

# ---- First average TEM data
print('======== locating TEM data for ens{}, {} Tg...'.format(Nens, massMag))
# loop over tracers
for qi in [0, 1, 2]:
    if(SKIP_TEM): continue
    if(qi > 0 and no_src_tag): continue # no tracers for the non-source tagged ensemble
    qstr  = ['','_TRACER-AOA','_TRACER-E90j'][int(qi)]
    print('====== working on tracer {}...'.format(qi))
    if(not cfb):
        enstem = sorted(glob.glob('{}/*{}*ens{}.eam.h1*_TEM_*L45{}.nc'.format(
                                    dailyloc, massStr, Nens, qstr)))[0]
    else:
        enstem = sorted(glob.glob('{}/*ens{}.cf.eam.h1*_TEM_*L45{}.nc'.format(
                                    dailyloc, Nens, qstr)))[0]
   
    print('time concatenated file: {}'.format(enstem))
    name = enstem.split('/')[-1].split('.nc')[0]
    name = name.split('199')[0] + name.split('000_')[-1]
    outfile = '{}/{}_10daily.nc'.format(outloc, name)

    if(os.path.isfile(outfile) and not overwrite):
        print('file exists and overwrite=False; skipping')
        continue
    
    print('---- computing 10-day average...')
    if(not dry):  
        # skip if exists
        existing = glob.glob(outfile)
        if(len(existing) > 0 and not overwrite): 
            print('10-day average file exists; skipping...')
        else:
            print('reading data...')
            time_avg = xr.open_dataset(enstem)
            if('time_bnds' in time_avg.data_vars): time_avg = time_avg.drop_vars('time_bnds')
            print('averaging data...')
            # average data in 10-day segments 
            #time_avg = time_avg.resample(time='10D').mean('time')
            time_avg = time_avg.rolling(time=10, min_periods=4, center=True).mean()
            time_avg = time_avg.isel(time=slice(9-ovr, None, 10-ovr))
            print('writing out...')
            time_avg.to_netcdf(outfile)
        
# ---- Now do the same for the zonal mean data

if(not SKIP_ZM):
    print('======== locating zonal mean data for ens{}, {} Tg...'.format(Nens, massMag))

    if(not cfb):
        enszm   = sorted(glob.glob('{}/*{}*ens{}.eam.h1*zonalmeans.nc'.format(
                                    dailyloc, massStr, Nens)))[0]
    else:
        enszm   = sorted(glob.glob('{}/*ens{}.cf.eam.h1*zonalmeans.nc'.format(
                                    dailyloc, Nens)))[0]

    print('time concatenated file: {}'.format(enszm))
    name = enszm.split('/')[-1].split('.nc')[0]
    name = name.split('199')[0] + name.split('000_')[-1]
    outfile = '{}/{}_10daily.nc'.format(outloc, name)
        
    if(os.path.isfile(outfile) and not overwrite):
        print('file exists and overwrite=False; skipping')
    else:
        print('---- computing 10-day average...')
        if(not dry): 
            # skip if exists
            existing = glob.glob(outfile)
            if(len(existing) > 0 and not overwrite): 
                print('10-day average file exists; skipping...')
            else:
                print('reading data...')
                time_avg = xr.open_dataset(enszm)
                if('time_bnds' in time_avg.data_vars): time_avg = time_avg.drop_vars('time_bnds')
                print('averaging data...')
                #time_avg = time_avg.resample(time='10D').mean('time')
                time_avg = time_avg.rolling(time=10, min_periods=4, center=True).mean()
                time_avg = time_avg.isel(time=slice(9-ovr, None, 10-ovr))
                print('writing out...')
                time_avg.to_netcdf(outfile)


# ---- Now do the same for the TEM budget data
print('======== locating TEM budget data for ens{}, {} Tg...'.format(Nens, massMag))
# loop over tracers
for qi in [0, 1, 2]:
    if(SKIP_BUDGET): continue
    if(qi > 0): continue # TEMPORARY FOR UTEND FIXES, SKIPPING RUNNING TRACERS FOR NOW, REMOVE LATER
    if(qi > 0 and no_src_tag): continue # no tracers for the non-source tagged ensemble
    qstr  = ['','_TRACER-AOA','_TRACER-E90j'][int(qi)]
    print('====== working on tracer {}...'.format(qi))
    if(not cfb):
        enstem = sorted(glob.glob('{}/*{}*ens{}.eam.h1*_TEMBudget{}.nc'.format(
                                    dailyloc, massStr, Nens, qstr)))[0]
    else:
        enstem = sorted(glob.glob('{}/*ens{}.cf.eam.h1*_TEMBudget{}.nc'.format(
                                    dailyloc, Nens, qstr)))[0]
   
    print('time concatenated file: {}'.format(enstem))
    name = enstem.split('/')[-1].split('.nc')[0]
    name = name.split('199')[0] + name.split('000_')[-1]
    outfile = '{}/{}_10daily.nc'.format(outloc, name)
    
    if(os.path.isfile(outfile) and not overwrite):
        print('file exists and overwrite=False; skipping')
        continue
    
    print('---- computing 10-day average...')
    if(not dry):  
        # skip if exists
        existing = glob.glob(outfile)
        if(len(existing) > 0 and not overwrite): 
            print('10-day average file exists; skipping...')
        else:
            print('reading data...')
            time_avg = xr.open_dataset(enstem)
            if('time_bnds' in time_avg.data_vars): time_avg = time_avg.drop_vars('time_bnds')
            print('averaging data...')
            # average data in 10-day segments 
            #time_avg = time_avg.resample(time='10D').mean('time')
            time_avg = time_avg.rolling(time=10, min_periods=4, center=True).mean()
            time_avg = time_avg.isel(time=slice(9-ovr, None, 10-ovr))
            print('writing out...')
            time_avg.to_netcdf(outfile)
