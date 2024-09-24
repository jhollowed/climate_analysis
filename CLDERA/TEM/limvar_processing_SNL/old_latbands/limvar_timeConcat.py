'''
Joe Hollowed 
University of Michigan 2024

This script concatenates reduced limvar ensemble member datasets in time to form
single netcdf files per ensemble member. Specifically, the TEM output and zonal means.
Command line arguments are:

Usage
-----
python ./limvar_zinterp.py [Nens] [histnum] [massMag] [tmini] [tmaxi] [overwrite] [dry]

Parameters
----------
Nens : int
    the ensemble member specified as an integer
histnum : int
    history file number to use
    0 = monthly data
    1 = daily data
massMag : int
    the mass magnitude of the ensemble. Valid values are:
    0, 1, 3, 5, 7, 10, 13, 15
    if massMag = 0, this sepcifies the counterfactual ensemble
tmini : int
    index of the first history file to include in the concatenation
tmaxi : int
    index of the last history file to include in the concatenation
overwrite : int
    whether or not to overwrite the concatenated data if it already exists
    at the specified outpit file. 0 = False, 1 = True, 
    If False, then the script skips the concatenation for any file that 
    already exists
dry : int
    whether or not to do a "dry run" of the function, printing information
    about the execution but skipping the concatenation. 0 = False, 1 = True

Note that tmini, tmaxi specify only the *index* of the identified time files. 
Each file may not cover the same amount of model time, depending on the histnum. 
However, for the limvar ensembles, both the monthly (histnum=0) and daily (histnum=1) 
history files include a single month of time, and thus a given combination of 
(tmini,tmaxi) will specify the same period of simulated time for either choice of 
histnum. 
'''

import sys
import pdb
import glob
import dask
import copy
import os.path
import numpy as np
import xarray as xr
from datetime import datetime

zmloc = '/gpfs/cldera/data/e3sm/e3sm-cldera/CLDERA-historical/limvar_src_tag/tem_processed/interp_to_pressure_zonalMeans'
temloc = '/gpfs/cldera/data/e3sm/e3sm-cldera/CLDERA-historical/limvar_src_tag/tem_processed/tem_output'
outloc = '/ascldap/users/jphollo/data/limvar/limvar_daily' 

Nens      = int(sys.argv[1])
histnum   = int(sys.argv[2])
massMag   = int(sys.argv[3])
tmini     = int(sys.argv[4])
tmaxi     = int(sys.argv[5])
overwrite = bool(int(sys.argv[6]))
dry       = bool(int(sys.argv[7]))
print('args: {}'.format(sys.argv))

cfb = massMag == 0   # counterfactual flag

skip_tem = 0 # flag to skip TEM data and concatenate only the zonal-mean model variables

# ----------------------------------------------------------------

if(histnum == 1 and not skip_tem):

    # ---- First concatenate TEM data, if operating on daily data
    print('locating TEM data for ens{}, {} Tg...'.format(Nens, massMag))
    # loop over tracers
    for qi in [0, 1, 2]:
        qstr  = ['','_TRACER-AOA','_TRACER-E90j'][int(qi)]
        print('working on tracer {}...'.format(qi))
        if(not cfb):
            enshist = sorted(glob.glob('{}/*{}Tg*ens{}.eam.h{}*TEM*L45{}.nc'.format(
                                        temloc, massMag, Nens, histnum, qstr)))
        else:
            enshist = sorted(glob.glob('{}/*ens{}.cf.eam.h{}*TEM*L45{}.nc'.format(
                                        temloc, Nens, histnum, qstr)))

        print('first file: {}'.format(enshist[0]))
        print('last file: {}'.format(enshist[-1]))
        if(not dry): 
            print('reading data...')
            enshistdat = [xr.open_dataset(hist) for hist in enshist]
            print('combining data...')
            ensconcat = xr.concat(enshistdat, dim='time')
            print('writing out...')
            name = enshist[0].split('/')[-1]
            ensconcat.to_netcdf('{}/{}'.format(outloc, name))

if(histnum == 0 or histnum == 1):

    # ---- Now do the same for the zonal mean data
    print('locating zonal mean data for ens{}, {} Tg, h{}...'.format(Nens, massMag, histnum))
    if(not cfb):
        enshist = sorted(glob.glob('{}/*{}Tg*ens{}.eam.h{}*.nc'.format(
                                    zmloc, massMag, Nens, histnum)))
    else:
        enshist = sorted(glob.glob('{}/*ens{}.cf.eam.h{}*.nc'.format(
                                    zmloc, Nens, histnum)))

    print('first file: {}'.format(enshist[0]))
    print('last file: {}'.format(enshist[-1]))
    if(not dry):
        print('reading data...')
        enshistdat = [xr.open_dataset(hist) for hist in enshist]
        print('combining data...')
        ensconcat = xr.concat(enshistdat, dim='time')
        print('writing out...')
        name = enshist[0].split('/')[-1]
        ensconcat.to_netcdf('{}/{}'.format(outloc, name))
