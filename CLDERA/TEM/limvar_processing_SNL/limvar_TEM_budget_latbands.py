'''
Joe Hollowed 
University of Michigan 2024

This script computes the TEM budget for latitude-band averaged limvar data. That is, the
total tendency expected by the TEM terms plus the parameterized gravity wave forcing
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

dataloc = '/ascldap/users/jphollo/data/limvar/limvar_latbands' 
outloc  = '/ascldap/users/jphollo/data/limvar/limvar_latbands_tembudget'

Nens      = int(sys.argv[1])
massMag   = int(sys.argv[2])
tmini     = int(sys.argv[3])
tmaxi     = int(sys.argv[4])
overwrite = bool(int(sys.argv[5]))
dry       = bool(int(sys.argv[6]))
print('args: {}'.format(sys.argv))

cfb = massMag == 0   # counterfactual flag

bands = ['SHpole', 'SHmid', 'tropics', 'NHmid', 'NHpole']

# -----------------------------------------------------------------

for band in bands:

    print('+++++++++++ working on band {}'.format(band))

    for qi in [0, 1, 2]:
        qstr  = ['','_TRACER-AOA','_TRACER-E90j'][int(qi)]
        print('====== working on tracer {}...'.format(qi))
        if(not cfb):
            temfile = sorted(glob.glob('{}/*{}Tg*ens{}*.eam.h1*TEM*L45{}_latband_{}.nc'.format(
                                        dataloc, massMag, Nens, qstr, band)))[0]
            zmfile  = sorted(glob.glob('{}/*{}Tg*ens{}*.eam.h1*zonal*latband_{}.nc'.format(
                                        dataloc, massMag, Nens, band)))[0]
        else:
            temfile = sorted(glob.glob('{}/*ens{}*.cf.eam.h1*TEM*L45{}_latband_{}.nc'.format(
                                        dataloc, Nens, qstr, band)))[0]
            zmfile  = sorted(glob.glob('{}/*ens{}*.cf.eam.h1*zonal*latband_{}.nc'.format(
                                        dataloc, Nens, band)))[0]
        
        outfile = '{}/{}_TEMBudget{}.nc'.format(outloc, zmfile.split('/')[-1].split('.nc')[0], qstr)

        zm = xr.open_dataset(zmfile)
        tem = xr.open_dataset(temfile)
        orig_vars = list(zm.data_vars)

        if(qi == 0):
            
            # compute total residual velocity U tendency
            zm['UTRESVEL'] = tem['utendvtem'] + tem['utendwtem']
            zm['UTRESVEL'].attrs['long name'] = 'sum of zonal mean utendvtem, utendwtem'
            zm['UTRESVEL'].attrs['units'] = 'm/s'

            # compute total U tendency
            zm['UTTOTAL_NOGW'] = tem['utendepfd'] + tem['utendvtem'] + tem['utendwtem']
            zm['UTTOTAL_NOGW'].attrs['long name'] = 'sum of zonal mean UTRESVEL, utendepfd '\
                                                    '(not inckuding gravity waves)'
            zm['UTTOTAL_NOGW'].attrs['units'] = 'm/s2'

            # compute difference of TEM U tendency and measured U tendency
            zm['UTDIFF_NOGW'] = zm['UTTOTAL_NOGW'] - zm['UTEND']
            zm['UTDIFF_NOGW'].attrs['long name'] = 'difference of UTTOTAL_NOGW and UTEND'
            zm['UTDIFF_NOGW'].attrs['units'] = 'm/s2'

        else:
            
            # compute total residual velocity tracer tendency
            zm['QTRESVEL'] = tem['qtendvtem'] + tem['qtendwtem']
            zm['QTRESVEL'].attrs['long name'] = 'sum of zonal mean qtendvtem, qtendwtem'
            zm['QTRESVEL'].attrs['units'] = '1/s'
            
            # compute total tracer tendency
            zm['QTTOTAL'] = tem['qtendetfd'] + tem['qtendvtem'] + tem['qtendwtem']
            zm['QTTOTAL'].attrs['long name'] = 'sum of zonal mean qtendepfd, QTRESVEL'
            zm['QTTOTAL'].attrs['units'] = '1/s'
            
            # template for total sources/sinks
            zm['QTSOURCE'] = zm['AOA'] * 0 
            zm['QTSOURCE'].attrs['long name'] = 'tracer source'
            zm['QTSOURCE'].attrs['units'] = '1/s'
            zm['QTSINK'] = zm['AOA'] * 0 
            zm['QTSINK'].attrs['long name'] = 'tracer sink'
            zm['QTSINK'].attrs['units'] = '1/s'
            zm['QTSRCSNK'] = zm['AOA'] * 0 
            zm['QTSRCSNK'].attrs['long name'] = 'sum of zonal mean QTSOURCE, QTSINK'
            zm['QTSRCSNK'].attrs['units'] = '1/s'
            
            if(qi == 1):
                # compute AOA tendency by parameterized tracer source
                # this won't be exactly right below 700 hPa since we've interpolated 
                # to pure pressure levels...
                zm['QTSOURCE'] = zm['AOA']/zm['AOA'] # clock tracer above 700 hPa, 1 s/s
                zm['QTSOURCE'].loc[{'plev':slice(700, 10000)}] = 0
                # AOA has no sinks
                zm['QTSINK'] = zm['AOA'] * 0
                # sum sources, sinks
                zm['QTSRCSNK'] = zm['QTSOURCE'] + zm['QTSINK']
                
                # add this to the total tendency computed above
                zm['QTTOTAL'] = zm['QTTOTAL'] + zm['QTSRCSNK']
                zm['QTTOTAL'].attrs['long name'] = 'sum of zonal mean qtendepfd, QTRESVEL, QTSRCSNK'

                # compute difference of TEM AOA tendency and measured AOA tendency
                zm['QTDIFF'] = zm['QTTOTAL'] - zm['AOATEND']
                zm['QTDIFF'].attrs['long name'] = 'difference of QTTOTAL and AOATEND'
                zm['QTDIFF'].attrs['units'] = '1/s'

            if(qi == 2):
                # compute E90 tendency by parameterized tracer sources and sinks
                # E90 has no source
                zm['QTSOURCE'] = zm['E90j'] * 0
                # e-folding decay with timescale of 90 days
                tau = 90 * 24 * 60 * 60
                zm['QTSINK'] = -1/tau * zm['E90j']
                # sum sources, sinks
                zm['QTSRCSNK'] = zm['QTSOURCE'] + zm['QTSINK']
                
                # add this to the total tendency computed above
                zm['QTTOTAL'] = zm['QTTOTAL'] + zm['QTSRCSNK']
                zm['QTTOTAL'].attrs['long name'] = 'sum of zonal mean qtendepfd, QTRESVEL, QTSRCSNK'

                # compute difference of TEM E90 tendency and measured E90 tendency
                zm['QTDIFF'] = zm['QTTOTAL'] - zm['E90TEND']
                zm['QTDIFF'].attrs['long name'] = 'difference of QTTOTAL and E90TEND'
                zm['QTDIFF'].attrs['units'] = '1/s'

        print('writing out...')
        zm = zm.drop_vars(orig_vars)
        zm.to_netcdf(outfile)

















