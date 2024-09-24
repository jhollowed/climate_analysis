'''
Joe Hollowed 
University of Michigan 2024

This script computes meridional-band means of limvar TEM data, and limvar zonal means, and
averages the result in time over 10-day intervals.
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

dailyloc = '/ascldap/users/jphollo/data/limvar/limvar_daily'
bandloc   = '/ascldap/users/jphollo/data/limvar/limvar_latbands' 
outloc   = '/ascldap/users/jphollo/data/limvar/limvar_latbands_10daily' 

print('args: {}'.format(sys.argv))
Nens      = int(sys.argv[1])
massMag   = int(sys.argv[2])
tmini     = int(sys.argv[3])
tmaxi     = int(sys.argv[4])
overwrite = bool(int(sys.argv[5]))
dry       = bool(int(sys.argv[6]))

cfb = massMag == 0   # counterfactual flag

lat_band_names = ['SHpole', 'SHmid', 'tropics', 'eq', 'NHmid', 'NHpole']


# ----------------------------------------------------------------
    
for i in range(len(lat_band_names)):
    print('\n============== working on lat band {}'.format(lat_band_names[i]))

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
        band_file = '{}/{}_latband_{}.nc'.format(bandloc, name, lat_band_names[i])
        outfile = '{}/{}_latband_{}_10daily.nc'.format(outloc, name, lat_band_names[i])
        
        print('---- computing 10-day average for band {}...'.format(lat_band_names[i]))
        if(not dry): 
            
            print('reading data...')
            band_avg = xr.open_dataset(band_file)
            
            # skip if exists
            existing = glob.glob(outfile)
            if(len(existing) > 0 and not overwrite): 
                print('{} lat band 10-day average file exists; skipping...'.format(
                                                                 lat_band_names[i]))
            else:
                print('averaging data...')
                # average data in 10-day segments, 
                band_avg = band_avg.resample(time='10D').mean('time')
                print('writing out...')
                band_avg.to_netcdf(outfile)
            
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
    band_file = '{}/{}_latband_{}.nc'.format(bandloc, name, lat_band_names[i])
    outfile = '{}/{}_latband_{}_10daily.nc'.format(outloc, name, lat_band_names[i])

    print('---- computing 10-day average for band {}...'.format(lat_band_names[i]))
    if(not dry): 
        
        print('reading data...')
        band_avg = xr.open_dataset(band_file)
        
        # skip if exists
        existing = glob.glob(outfile)
        if(len(existing) > 0 and not overwrite): 
            print('{} lat band 10-day average file exists; skipping...'.format(
                                                             lat_band_names[i]))
        else:
            print('averaging data...')
            band_avg = band_avg.resample(time='10D').mean('time')
            
            print('writing out...')
            band_avg.to_netcdf(outfile)

    # ADD CODE TO ALSO WRITEOUT THE 10-DAY AVERAGE OF THIS TIME SERIES
