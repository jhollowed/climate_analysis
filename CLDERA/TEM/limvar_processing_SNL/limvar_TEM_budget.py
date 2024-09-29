'''
Joe Hollowed 
University of Michigan 2024

This script computes the TEM budget for monthly mean limvar data. That is, the
total tendency expected by the TEM terms plus the parameterized gravity wave forcing
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
from scipy import integrate

Nens      = int(sys.argv[1])
histnum   = int(sys.argv[2])
massMag   = int(sys.argv[3])
tmini     = int(sys.argv[4])
tmaxi     = int(sys.argv[5])
overwrite = bool(int(sys.argv[6]))
dry       = bool(int(sys.argv[7]))
print('args: {}'.format(sys.argv))

cfb = massMag == 0   # counterfactual flag

# set data paths
if(histnum == 0):
    dataloc = '/ascldap/users/jphollo/data/limvar/limvar_monthly' 
    outloc  = '/ascldap/users/jphollo/data/limvar/limvar_monthly'
    sfx = '_monthlymean'
elif(histnum == 1):
    #dataloc = '/ascldap/users/jphollo/data/limvar/limvar_daily' 
    #outloc  = '/ascldap/users/jphollo/data/limvar/limvar_daily'
    dataloc = '/gpfs/cldera/data/e3sm/e3sm-cldera/CLDERA-historical/limvar_src_tag/tem_processed/timeConcat_zonalMeans'
    outloc = '/gpfs/cldera/data/e3sm/e3sm-cldera/CLDERA-historical/limvar_src_tag/tem_processed/timeConcat_zonalMeans'
    sfx = ''
else:
    raise RuntimeError('histnum must be 0 or 1')

spd = 24*60*60  # seconds per day
tol = 1e-5      # integrated wind sanity check tolerance

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


# -----------------------------------------------------------------

for qi in [0, 1, 2]:
    if(qi > 0): continue # TEMPORARY FOR UTEND FIXES, SKIPPING RUNNING TRACERS FOR NOW, REMOVE LATER
    if(qi > 0 and no_src_tag): continue # no tracers for the non-source tagged ensemble
    qstr  = ['','_TRACER-AOA','_TRACER-E90j'][int(qi)]
    print('====== working on tracer {}...'.format(qi))
    if(not cfb):
        temfile = sorted(glob.glob('{}/*{}*ens{}.eam.h1*TEM*L45{}{}.nc'.format(
                                    dataloc, massStr, Nens, qstr, sfx)))[0]
        zmfile  = sorted(glob.glob('{}/*{}*ens{}.eam.h1*zonal*{}.nc'.format(
                                    dataloc, massStr, Nens, sfx)))[0]
    else:
        temfile = sorted(glob.glob('{}/*ens{}.cf.eam.h1*TEM*L45{}{}.nc'.format(
                                    dataloc, Nens, qstr, sfx)))[0]
        zmfile  = sorted(glob.glob('{}/*ens{}.cf.eam.h1*zonal*{}.nc'.format(
                                    dataloc, Nens, sfx)))[0]
    
    outfile = '{}/{}_TEMBudget{}.nc'.format(outloc, zmfile.split('/')[-1].split('.nc')[0], qstr)
    print('outfile is {}'.format(outfile))

    if(os.path.isfile(outfile) and not overwrite):
        print('file exists and overwrite=False; skipping')
        continue

    zm = xr.open_dataset(zmfile)
    tem = xr.open_dataset(temfile)
    orig_vars = list(zm.data_vars)

    if(qi == 0):
       
        if(histnum == 0):
            # compute total gravity wave U tendency; these variables only exist for monthly data
            zm['UTGWTOTAL'] = zm['BUTGWSPEC'] + zm['UTGWSPEC'] + zm['UTGWORO']
            zm['UTGWTOTAL'].attrs['long name'] = 'sum of zonal mean BUTGWSPEC, UTGWSPEC, UTGWORO'
            zm['UTGWTOTAL'].attrs['units'] = 'm/s2'

        # compute total residual velocity U tendency
        zm['UTRESVEL'] = tem['utendvtem'] + tem['utendwtem']
        zm['UTRESVEL'].attrs['long name'] = 'sum of zonal mean utendvtem, utendwtem'
        zm['UTRESVEL'].attrs['units'] = 'm/s'

        # compute total U tendency
        if(histnum == 0):
            # monthly data; total includes gravity waves
            zm['UTTOTAL'] = zm['UTGWTOTAL'] + tem['utendepfd'] + tem['utendvtem'] + tem['utendwtem']
            zm['UTTOTAL'].attrs['long name'] = 'sum of zonal mean UTGWTOTAL, UTRESVEL, utendepfd'
            zm['UTTOTAL'].attrs['units'] = 'm/s2'
            # compute difference of TEM U tendency and measured U tendency
            zm['UTDIFF'] = zm['UTTOTAL'] - zm['UTEND']
            zm['UTDIFF'].attrs['long name'] = 'difference of UTTOTAL and UTEND'
            zm['UTDIFF'].attrs['units'] = 'm/s2'
        if(histnum == 1):
            # daily data; total does not include gravity waves
            zm['UTTOTAL_NOGW'] = tem['utendepfd'] + tem['utendvtem'] + tem['utendwtem']
            zm['UTTOTAL_NOGW'].attrs['long name'] = 'sum of zonal mean UTRESVEL, utendepfd'
            zm['UTTOTAL_NOGW'].attrs['units'] = 'm/s2'
            # compute difference of TEM U tendency and measured U tendency
            zm['UTDIFF'] = zm['UTTOTAL_NOGW'] - zm['UTEND']
            zm['UTDIFF'].attrs['long name'] = 'difference of UTTOTAL_NOGW and UTEND'
            zm['UTDIFF'].attrs['units'] = 'm/s2'
            
            # ----- get integrated tendencies by starting from U on each day

            # for integration methods which take N point and output N-1 points as a 
            # numpy array (scipy), this function converts them back to an xr dataarray
            # with dimensions and coordinates consistent with the data, and adds back
            # the missing point at the beginning of the time series. That first point
            # added will be all zeros, and it will be removed when later passed to
            # merge_U_shifted()
            def values_to_xr(values):
                out = xr.zeros_like(zm['U'])
                out[{'time':slice(1, None)}] = values
                return out
            # this function takes an integrated tendency, shfits it in time by one index, 
            # and adds it to the zonal-mean zonal wind, such that U at t=0 is retained as 
            # the initial condition, and the integrated tendency at t=-1 (which is 
            # probably nan) is dropped
            def merge_U_shifted(inttend):
                out = zm['U'].copy()
                out[{'time':slice(1, None)}] = inttend.isel(time=slice(None, -1)).values
                return out

            # get integrated U components
            # first sanity check with UTEND
            tmp = zm['U'] + zm['UTEND'] * spd
            zm['UTEND_INT'] = merge_U_shifted(tmp)
            if(np.max(np.abs(zm['UTEND_INT'] - zm['U'])) > tol):
                raise RuntimeError('UTEND integration failed sanity check')
            zm['UTEND_INT'].attrs['long name'] = 'integrated UTEND'
            zm['UTEND_INT'].attrs['units'] = 'm/s'
            
            tmp = zm['U'] + zm['UTRESVEL'] * spd 
            zm['UTRESVEL_INT'] = merge_U_shifted(tmp)
            zm['UTRESVEL_INT'].attrs['long name'] = 'integrated UTRESVEL'
            zm['UTRESVEL_INT'].attrs['units'] = 'm/s'
            
            tmp = zm['U'] + tem['utendepfd'] * spd 
            zm['UTEPFD_INT'] = merge_U_shifted(tmp)
            zm['UTEPFD_INT'].attrs['long name'] = 'integrated utendepfd'
            zm['UTEPFD_INT'].attrs['units'] = 'm/s'
            
            tmp = zm['U'] + zm['UTDIFF'] * spd 
            zm['UTDIFF_INT'] = merge_U_shifted(tmp)
            zm['UTDIFF_INT'][{'time':slice(1, None)}] = tmp.isel(time=slice(None, -1)).values
            zm['UTDIFF_INT'].attrs['long name'] = 'integrated UTDIFF'
            zm['UTDIFF_INT'].attrs['units'] = 'm/s' 
            
            # ----- get integrated tendencies by starting from U on day 0 only, summing tendencies
            u0 = zm['U'].isel(time=0, drop=True) # initial condition

            tmp = u0 + (zm['UTEND'] * spd).cumsum('time')
            zm['UTEND_INT2'] = merge_U_shifted(tmp)
            zm['UTEND_INT2'].attrs['long name'] = 'integrated UTEND, method 2'
            zm['UTEND_INT2'].attrs['units'] = 'm/s'
            
            tmp = u0 + (zm['UTRESVEL'] * spd).cumsum('time')
            zm['UTRESVEL_INT2'] = merge_U_shifted(tmp)
            zm['UTRESVEL_INT2'].attrs['long name'] = 'integrated UTRESVEL, method 2'
            zm['UTRESVEL_INT2'].attrs['units'] = 'm/s'
            
            tmp = u0 + (tem['utendepfd'] * spd).cumsum('time')
            zm['UTEPFD_INT2'] = merge_U_shifted(tmp)
            zm['UTEPFD_INT2'].attrs['long name'] = 'integrated utendepfd, method 2'
            zm['UTEPFD_INT2'].attrs['units'] = 'm/s'
            
            tmp = u0 + (zm['UTDIFF'] * spd).cumsum('time')
            zm['UTDIFF_INT2'] = merge_U_shifted(tmp)
            zm['UTDIFF_INT2'].attrs['long name'] = 'integrated UTDIFF, method 2'
            zm['UTDIFF_INT2'].attrs['units'] = 'm/s' 
            
            # ----- get integrated tendencies by starting from U on day 0 only, using trapezoid method
            trapz = integrate.cumulative_trapezoid

            tidx = zm['UTEND'].dims.index('time')
            tmp = u0 + values_to_xr(trapz(zm['UTEND'].values, dx=spd, axis=tidx))
            zm['UTEND_INT3'] = merge_U_shifted(tmp)
            zm['UTEND_INT3'].attrs['long name'] = 'integrated UTEND, method 3'
            zm['UTEND_INT3'].attrs['units'] = 'm/s'
            
            tidx = zm['UTRESVEL'].dims.index('time')
            tmp = u0 + values_to_xr(trapz(zm['UTRESVEL'].values, dx=spd, axis=tidx))
            zm['UTRESVEL_INT3'] = merge_U_shifted(tmp)
            zm['UTRESVEL_INT3'].attrs['long name'] = 'integrated UTRESVEL, method 3'
            zm['UTRESVEL_INT3'].attrs['units'] = 'm/s'
            
            tidx = tem['utendepfd'].dims.index('time')
            tmp = u0 + values_to_xr(trapz(tem['utendepfd'].values, dx=spd, axis=tidx))
            zm['UTEPFD_INT3'] = merge_U_shifted(tmp)
            zm['UTEPFD_INT3'].attrs['long name'] = 'integrated utendepfd, method 3'
            zm['UTEPFD_INT3'].attrs['units'] = 'm/s'
            
            tidx = zm['UTDIFF'].dims.index('time')
            tmp = u0 + values_to_xr(trapz(zm['UTDIFF'].values, dx=spd, axis=tidx))
            zm['UTDIFF_INT3'] = merge_U_shifted(tmp)
            zm['UTDIFF_INT3'].attrs['long name'] = 'integrated UTDIFF, method 3'
            zm['UTDIFF_INT3'].attrs['units'] = 'm/s' 
            
            # ----- get integrated tendencies by starting from U on day 0 only, using simpson method
            simp = integrate.cumulative_simpson

            tidx = zm['UTEND'].dims.index('time')
            tmp = u0 + values_to_xr(simp(zm['UTEND'].values, dx=spd, axis=tidx))
            zm['UTEND_INT4'] = merge_U_shifted(tmp)
            zm['UTEND_INT4'].attrs['long name'] = 'integrated UTEND, method 4'
            zm['UTEND_INT4'].attrs['units'] = 'm/s'
            
            tidx = zm['UTRESVEL'].dims.index('time')
            tmp = u0 + values_to_xr(simp(zm['UTRESVEL'].values, dx=spd, axis=tidx))
            zm['UTRESVEL_INT4'] = merge_U_shifted(tmp)
            zm['UTRESVEL_INT4'].attrs['long name'] = 'integrated UTRESVEL, method 4'
            zm['UTRESVEL_INT4'].attrs['units'] = 'm/s'
            
            tidx = tem['utendepfd'].dims.index('time')
            tmp = u0 + values_to_xr(simp(tem['utendepfd'].values, dx=spd, axis=tidx))
            zm['UTEPFD_INT4'] = merge_U_shifted(tmp)
            zm['UTEPFD_INT4'].attrs['long name'] = 'integrated utendepfd, method 4'
            zm['UTEPFD_INT4'].attrs['units'] = 'm/s'
            
            tidx = zm['UTDIFF'].dims.index('time')
            tmp = u0 + values_to_xr(simp(zm['UTDIFF'].values, dx=spd, axis=tidx))
            zm['UTDIFF_INT4'] = merge_U_shifted(tmp)
            zm['UTDIFF_INT4'].attrs['long name'] = 'integrated UTDIFF, method 4'
            zm['UTDIFF_INT4'].attrs['units'] = 'm/s' 

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
        # template for tendency difference
        zm['QTDIFF'] = zm['AOA'] * 0 
        zm['QTDIFF'].attrs['long name'] = 'difference of QTTOTAL and AOATEND'
        zm['QTDIFF'].attrs['units'] = '1/s'
        # templates for integrated tendencies
        zm['QTRESVEL_INT'] = zm['AOA'] * 0
        zm['QTRESVEL_INT'].attrs['long name'] = 'integrated QTRESVEL'
        zm['QTRESVEL_INT'].attrs['units'] = 'kg/kg'
        zm['QTETFD_INT'] = zm['AOA'] * 0
        zm['QTETFD_INT'].attrs['long name'] = 'integrated qtendetfd'
        zm['QTETFD_INT'].attrs['units'] = 'kg/kg'
        zm['QTDIFF_INT'] = zm['AOA'] * 0
        zm['QTDIFF_INT'].attrs['long name'] = 'integrated QTDIFF'
        zm['QTDIFF_INT'].attrs['units'] = 'kg/kg'
        
        if(qi == 1):
            # compute AOA tendency by parameterized tracer source
            # this won't be exactly right below 700 hPa since we've interpolated 
            # to pure pressure levels...
            zm['QTSOURCE'] = zm['AOA']/zm['AOA'] # clock tracer above 700 hPa, 1 s/s
            zm['QTSOURCE'] = zm['QTSOURCE'] / 86400 # scale clock tracer to units of day/s
            zm['QTSOURCE'].loc[{'plev':slice(700, 10000)}] = 0
            zm['QTSOURCE'].attrs['units'] = 'day/s'
            # AOA has no sinks
            zm['QTSINK'] = zm['AOA'] * 0
            zm['QTSINK'].attrs['units'] = 'day/s'
            # sum sources, sinks
            zm['QTSRCSNK'] = zm['QTSOURCE'] + zm['QTSINK']
            zm['QTSRCSNK'].attrs['units'] = 'day/s'
            # add this to the total tendency computed above
            zm['QTTOTAL'] = zm['QTTOTAL'] + zm['QTSRCSNK']
            zm['QTTOTAL'].attrs['units'] = 'day/s'
            # compute difference of TEM AOA tendency and measured AOA tendency
            zm['QTDIFF'] = zm['QTTOTAL'] - zm['AOATEND']
            zm['QTDIFF'].attrs['units'] = 'day/s'
       
            if(histnum == 1):
                # get integrated AOA tendency components
                # first sanity check on total AOATEND
                tmp = zm['AOA'] + zm['AOATEND'] * spd 
                tmp2 = zm['AOA']
                tmp2[{'time':slice(1, None)}] = tmp.isel(time=slice(None, -1)).values
                if(np.sum(np.abs(tmp2 - zm['AOA'])) > tol):
                    raise RuntimeError('AOATEND integration failed sanity check')
                
                tmp = zm['AOA'] + zm['QTRESVEL'] * spd 
                zm['QTRESVEL_INT'] = zm['AOA'].copy()
                zm['QTRESVEL_INT'][{'time':slice(1, None)}] = tmp.isel(time=slice(None, -1)).values
                zm['QTRESVEL_INT'].attrs['units'] = 'day'
                
                tmp = zm['AOA'] + tem['qtendetfd'] * spd 
                zm['UTETFD_INT'] = zm['AOA'].copy()
                zm['UTETFD_INT'][{'time':slice(1, None)}] = tmp.isel(time=slice(None, -1)).values
                zm['UTETFD_INT'].attrs['units'] = 'day'
                
                tmp = zm['AOA'] + zm['QTDIFF'] * spd 
                zm['QTDIFF_INT'] = zm['AOA'].copy()
                zm['QTDIFF_INT'][{'time':slice(1, None)}] = tmp.isel(time=slice(None, -1)).values
                zm['QTDIFF_INT'].attrs['units'] = 'day'

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
            # compute difference of TEM E90 tendency and measured E90 tendency
            zm['QTDIFF'] = zm['QTTOTAL'] - zm['E90TEND']
            
            if(histnum == 1):
                # get integrated E90 tendency components
                # first sanity check on total E90TEND
                tmp = zm['E90j'] + zm['E90TEND'] * spd 
                tmp2 = zm['E90j']
                tmp2[{'time':slice(1, None)}] = tmp.isel(time=slice(None, -1)).values
                if(np.sum(np.abs(tmp2 - zm['E90j'])) > tol):
                    raise RuntimeError('E90TEND integration failed sanity check')
                
                tmp = zm['E90j'] + zm['QTRESVEL'] * spd 
                zm['QTRESVEL_INT'] = zm['E90j'].copy()
                zm['QTRESVEL_INT'][{'time':slice(1, None)}] = tmp.isel(time=slice(None, -1)).values
                zm['QTRESVEL_INT'].attrs['units'] = 'kg/kg'
                
                tmp = zm['E90j'] + tem['qtendetfd'] * spd 
                zm['UTETFD_INT'] = zm['E90j'].copy()
                zm['UTETFD_INT'][{'time':slice(1, None)}] = tmp.isel(time=slice(None, -1)).values
                zm['UTETFD_INT'].attrs['units'] = 'kg/kg'
                
                tmp = zm['E90j'] + zm['QTDIFF'] * spd 
                zm['QTDIFF_INT'] = zm['E90j'].copy()
                zm['QTDIFF_INT'][{'time':slice(1, None)}] = tmp.isel(time=slice(None, -1)).values
                zm['QTDIFF_INT'].attrs['units'] = 'kg/kg'
            

    if(not dry):
        print('writing out...')
        zm = zm.drop_vars(orig_vars)
        zm.to_netcdf(outfile)
