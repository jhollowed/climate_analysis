'''
Joe Hollowed 
University of Michigan 2024

This script computes meridional-band  means of limvar TEM data, and limvar zonal means.
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
outloc   = '/ascldap/users/jphollo/data/limvar/limvar_latbands' 

Nens      = int(sys.argv[1])
massMag   = int(sys.argv[2])
tmini     = int(sys.argv[3])
tmaxi     = int(sys.argv[4])
overwrite = bool(int(sys.argv[5]))
dry       = bool(int(sys.argv[6]))
print('args: {}'.format(sys.argv))

cfb = massMag == 0   # counterfactual flag

lat_bands = [slice(-90, -60), slice(-60, -20), 
             slice(-20, 20), slice(-5, 5),
             slice(20, 60), slice(60, 90)]
lat_band_names = ['SHpole', 'SHmid', 'tropics', 'eq', 'NHmid', 'NHpole']

def lat_band_mean(data, band_slice):

    try:
        data = data.drop_vars('time_bnds')
    except ValueError:
        pass
    
    data_weights = np.cos(np.deg2rad(data.lat))
    data_weights.name = 'weights'
    
    data_lat_band = data.sel(lat=band_slice)
    data_lat_band = data_lat_band.weighted(data_weights)
    data_lat_band = data_lat_band.mean('lat')

    return data_lat_band


# ----------------------------------------------------------------
    
for i in range(len(lat_bands)):
    print('\n============== working on lat band {}'.format(lat_band_names[i]))

    # ---- First average TEM data
    print('======== locating TEM data for ens{}, {} Tg...'.format(Nens, massMag))
    # loop over tracers
    for qi in [0, 1, 2]:
        qstr  = ['','_TRACER-AOA','_TRACER-E90j'][int(qi)]
        print('====== working on tracer {}...'.format(qi))
        if(not cfb):
            enstem = sorted(glob.glob('{}/*{}Tg*ens{}*.eam.h1*TEM*L45{}.nc'.format(
                                        dailyloc, massMag, Nens, qstr)))[0]
        else:
            enstem = sorted(glob.glob('{}/*ens{}*.cf.eam.h1*TEM*L45{}.nc'.format(
                                        dailyloc, Nens, qstr)))[0]
       
        print('time concatenated file: {}'.format(enstem))
        name = enstem.split('/')[-1].split('.nc')[0]
        name = name.split('199')[0] + name.split('000_')[-1]
        
        print('---- computing {} band mean...'.format(lat_band_names[i]))
        if(not dry): 
            
            print('reading data...')
            tem = xr.open_dataset(enstem)
            
            # skip if exists
            outfile = '{}/{}_latband_{}.nc'.format(outloc, name, lat_band_names[i])
            existing = glob.glob(outfile)
            if(len(existing) > 0 and not overwrite): 
                print('{} lat band file exists; skipping...'.format(lat_band_names[i]))
            else:
                print('averaging data...')
                tem = lat_band_mean(tem, lat_bands[i])
                print('writing out...')
                tem.to_netcdf(outfile)
        

    # ---- Now do the same for the zonal mean data
    print('======== locating zonal mean data for ens{}, {} Tg...'.format(Nens, massMag))

    if(not cfb):
        enszm   = sorted(glob.glob('{}/*{}Tg*ens{}*.eam.h1*zonalmeans.nc'.format(
                                    dailyloc, massMag, Nens)))[0]
    else:
        enszm   = sorted(glob.glob('{}/*ens{}*.cf.eam.h1*zonalmeans.nc'.format(
                                    dailyloc, Nens, qstr)))[0]

    print('time concatenated file: {}'.format(enstem))
    name = enszm.split('/')[-1].split('.nc')[0]
    name = name.split('199')[0] + name.split('000_')[-1]


    print('---- computing {} band mean...'.format(lat_band_names[i]))
    if(not dry): 
        
        print('reading data...')
        zm   = xr.open_dataset(enszm)
        
        # skip if exists
        outfile = '{}/{}_latband_{}.nc'.format(outloc, name, lat_band_names[i])
        existing = glob.glob(outfile)
        if(len(existing) > 0 and not overwrite): 
            print('{} lat band file exists; skipping...'.format(lat_band_names[i]))
        else:
            print('averaging data...')
            zm = lat_band_mean(zm, lat_bands[i])
            
            print('writing out...')
            zm.to_netcdf(outfile)
