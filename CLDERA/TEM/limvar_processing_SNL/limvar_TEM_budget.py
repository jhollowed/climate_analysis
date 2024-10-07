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
    dataloc = '/gpfs/cldera/data/e3sm/e3sm-cldera/CLDERA-historical/limvar_src_tag/tem_processed/timeConcat_zonalMeans'
    outloc = '/gpfs/cldera/data/e3sm/e3sm-cldera/CLDERA-historical/limvar_src_tag/tem_processed/timeConcat_zonalMeans'
    sfx = ''
else:
    raise RuntimeError('histnum must be 0 or 1')

spd = 24*60*60  # seconds per day
tol = 1e-4      # integrated wind sanity check tolerance in m/s

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
            zm['utendgw'] = zm['BUTGWSPEC'] + zm['UTGWSPEC'] + zm['UTGWORO']
            zm['utendgw'].attrs['long_name'] = 'sum of zonal mean BUTGWSPEC, UTGWSPEC, UTGWORO'
            zm['utendgw'].attrs['units'] = 'm/s2'

        # compute total residual velocity U tendency
        zm['utendresvel'] = tem['utendvtem'] + tem['utendwtem']
        zm['utendresvel'].attrs['long_name'] = 'sum of zonal mean utendvtem, utendwtem'
        zm['utendresvel'].attrs['units'] = 'm/s'
        # also compute the gradient-normal vectors of the streamfunction (parallel to contours)
        
        def sph_grad(x, phi, p):
            # phi in degrees, p in hPa
            latdim = x.dims.index('lat')
            pdim   = x.dims.index('plev')
            args   = {'dims':x.dims, 'coords':x.coords}
            # radius of earth, scale height in meters
            a, H = 6.37123e6, 7e3
            # compute gradient components
            dxdphi = 1/a * xr.DataArray(np.gradient(x, np.deg2rad(phi), axis=latdim), **args)
            dxdp   = -(p/H) * xr.DataArray(np.gradient(x, p, axis=pdim), **args)
            dxdp   = dxdp.transpose(*args['dims'])
            return dxdphi, dxdp

        lat, plev, time = tem['psitem'].lat, tem['psitem'].plev, tem['psitem'].time
        gphi, gp   = sph_grad(tem['psitem'], lat, plev) # gradient vector components 
        gnphi, gnp = gp, gphi # gradient-normal vector; sign is clockwise around local maxima
        zm['psitem_gnlat'] = gnphi  
        zm['psitem_gnp']   = gnp
        zm['psitem_gnlat'].attrs['units']     = 'kg/s/m'
        zm['psitem_gnlat'].attrs['long_name'] = 'meridional component of the TEM streamfunction '\
                                                'gradient-normal vector field'
        zm['psitem_gnp'].attrs['units']       = 'kg/s/m'
        zm['psitem_gnp'].attrs['long_name']   = 'vertical component of the TEM streamfunction'\
                                                'gradient-normal vector field'

        # compute total U tendency
        if(histnum == 0):
            # monthly data; total includes gravity waves
            zm['utendtotal'] = zm['utendgw'] + tem['utendepfd'] + tem['utendvtem'] + tem['utendwtem']
            zm['utendtotal'].attrs['long_name'] = 'sum of zonal mean utendgw, utendresvel, utendepfd'
            zm['utendtotal'].attrs['units'] = 'm/s2'
            # compute difference of TEM U tendency and measured U tendency
            zm['utenddiff'] = zm['utendtotal'] - zm['UTEND']
            zm['utenddiff'].attrs['long_name'] = 'difference of utendtotal and UTEND'
            zm['utenddiff'].attrs['units'] = 'm/s2'
        if(histnum == 1):

            # first compute the total U tendency (this overwrites the tendencies computed
            # in limvar_zonalmeans.py, which erroneously place a nan at the end of each month...)
            # Those are fine to use for the monthly averaged data, since the nans will have 
            # been effectively removed in the monthly averaging process, but not for the daily data
            # If we recompute the tendencies here, then since the data is concatenated in time 
            # at this stage, there will only be one nan at the very end.
            print('---- computing tendencies...')
            tmp = zm['U'].diff('time', label='lower') # in m/s/day
            zm['UTEND'] = tmp/spd  # in m/s/s
            zm['UTEND'].name = 'UTEND'
            zm['UTEND'].attrs['long_name'] = 'U tendency'
            zm['UTEND'].attrs['units']     = 'm/s/s'
            orig_vars.remove('UTEND')
           
            # daily data; total does not include gravity waves
            zm['utendtotal_nogw'] = tem['utendepfd'] + tem['utendvtem'] + tem['utendwtem']
            zm['utendtotal_nogw'].attrs['long_name'] = 'sum of zonal mean utendresvel, utendepfd'
            zm['utendtotal_nogw'].attrs['units'] = 'm/s2'
            # compute difference of TEM U tendency and measured U tendency
            zm['utenddiff'] = zm['UTEND'] - zm['utendtotal_nogw']
            zm['utenddiff'].attrs['long_name'] = 'difference of utendtotal_nogw and UTEND'
            zm['utenddiff'].attrs['units'] = 'm/s2'
            
            # ================== integrate U tendencies ==================
            # The method below deposits the cumulative sum of all tendency data up to time N 
            # at time N+1 (in other words, UTEND_int at time N is the cumulative sum of UTEND 
            # at time N-1).
            # This ensures that the tendency at the time n does not contribute to the integrated 
            # value at time n.
            # This involves a forward shift of the data by one index, where the integrated 
            # quantitiy at time 0 is set to 0.
            # After applying this recipe, the initial condition is added uniformly to the 
            # integrated tendency across time, which is just U at time 0

            def integrate(x0, dxdt):
                out = spd * dxdt.cumsum('time')
                out = x0 + out.shift(time=1, fill_value=0)
                return out
            
            # --- initial condition
            u0 = zm['U'].isel(time=0, drop=True) # initial condition
            
            zm['UTEND_int'] = integrate(u0, zm['UTEND'])
            zm['UTEND_int'].attrs['long_name'] = 'integrated UTEND'
            zm['UTEND_int'].attrs['units'] = 'm/s'
            # for UTEND, this method should exactly recover U...
            if(np.max(np.abs(zm['UTEND_int'] - zm['U'])) > tol):
                raise RuntimeError('UTEND integration failed sanity check')
            
            zm['utendresvel_int'] = integrate(u0, zm['utendresvel'])
            zm['utendresvel_int'].attrs['long_name'] = 'integrated utendresvel'
            zm['utendresvel_int'].attrs['units'] = 'm/s'
            
            zm['utendvtem_int'] = integrate(u0, tem['utendvtem'])
            zm['utendvtem_int'].attrs['long_name'] = 'integrated utendvtem'
            zm['utendvtem_int'].attrs['units'] = 'm/s'
            
            zm['utendwtem_int'] = integrate(u0, tem['utendwtem'])
            zm['utendwtem_int'].attrs['long_name'] = 'integrated utendwtem'
            zm['utendwtem_int'].attrs['units'] = 'm/s'
            
            zm['utendepfd_int'] = integrate(u0, tem['utendepfd'])
            zm['utendepfd_int'].attrs['long_name'] = 'integrated utendepfd'
            zm['utendepfd_int'].attrs['units'] = 'm/s'
            
            zm['utenddiff_int'] = integrate(u0, zm['utenddiff'])
            zm['utenddiff_int'].attrs['long_name'] = 'integrated utenddiff'
            zm['utenddiff_int'].attrs['units'] = 'm/s'

            # sanity checks...
            dd = zm['UTEND_int'].isel(time=0)
            assert zm['UTEND_int'].isel(time=0).equals(dd), "UTEND_int broken IC!"
            assert zm['utendresvel_int'].isel(time=0).equals(dd), "utendresvel_int broken IC!"
            assert zm['utendvtem_int'].isel(time=0).equals(dd), "utendvtem_int broken IC!"
            assert zm['utendwtem_int'].isel(time=0).equals(dd), "utendwtem_int broken IC!"
            assert zm['utendepfd_int'].isel(time=0).equals(dd), "utendepfd_int broken IC!"
            assert zm['utenddiff_int'].isel(time=0).equals(dd), "utenddiff_int broken IC!"

    else:
        
        # as above for U, replace the total tendency with new ones computed here
        tmp = zm_aoa.diff('time', label='lower') # in day/day
        zm['AOATEND'] = tmp/spd  # in day/s
        zm['AOATEND'].name = 'AOATEND'
        zm['AOATEND'].attrs['long_name'] = 'AOA tendency'
        zm['AOATEND'].attrs['units']     = 'day/s'
       
        tmp = zm_e90.diff('time', label='lower') # in kg/kg/days
        zm['E90TEND'] = tmp/spd  # in kg/kg/s
        zm['E90TEND'].name = 'E90TEND'
        zm['E90TEND'].attrs['long_name'] = 'E90 tendency'
        zm['E90TEND'].attrs['units']     = 'kg/kg/s'
      
        # compute total residual velocity tracer tendency
        zm['qtendresvel'] = tem['qtendvtem'] + tem['qtendwtem']
        zm['qtendresvel'].attrs['long_name'] = 'sum of zonal mean qtendvtem, qtendwtem'
        zm['qtendresvel'].attrs['units'] = '1/s'
        
        # compute total tracer tendency
        zm['qtendtotal'] = tem['qtendetfd'] + tem['qtendvtem'] + tem['qtendwtem']
        zm['qtendtotal'].attrs['long_name'] = 'sum of zonal mean qtendetfd, qtendresvel'
        zm['qtendtotal'].attrs['units'] = '1/s'
        
        # template for total sources/sinks
        zm['qtendsource'] = zm['AOA'] * 0 
        zm['qtendsource'].attrs['long_name'] = 'tracer source'
        zm['qtendsource'].attrs['units'] = '1/s'
        zm['qtendsink'] = zm['AOA'] * 0 
        zm['qtendsink'].attrs['long_name'] = 'tracer sink'
        zm['qtendsink'].attrs['units'] = '1/s'
        # template for tendency difference
        zm['qtenddiff'] = zm['AOA'] * 0 
        zm['qtenddiff'].attrs['units'] = '1/s'
        # templates for integrated tendencies
        zm['qtendresvel_int'] = zm['AOA'] * 0
        zm['qtendresvel_int'].attrs['long_name'] = 'integrated qtendresvel'
        zm['qtendresvel_int'].attrs['units'] = 'kg/kg'
        zm['qtendetfd_int'] = zm['AOA'] * 0
        zm['qtendetfd_int'].attrs['long_name'] = 'integrated qtendetfd'
        zm['qtendetfd_int'].attrs['units'] = 'kg/kg'
        zm['qtenddiff_int'] = zm['AOA'] * 0
        zm['qtenddiff_int'].attrs['long_name'] = 'integrated qtenddiff'
        zm['qtenddiff_int'].attrs['units'] = 'kg/kg'
        
        if(qi == 1):
            # compute AOA tendency by parameterized tracer source
            # this won't be exactly right below 700 hPa since we've interpolated 
            # to pure pressure levels...
            zm['qtendsource'] = zm['AOA']/zm['AOA'] # clock tracer above 700 hPa, 1 s/s
            zm['qtendsource'] = zm['qtendsource'] / 86400 # scale clock tracer to units of day/s
            zm['qtendsource'].loc[{'plev':slice(700, 10000)}] = 0
            zm['qtendsource'].attrs['units'] = 'day/s'
            # AOA has no sinks
            zm['qtendsink'] = zm['AOA'] * 0
            zm['qtendsink'].attrs['units'] = 'day/s'
            # add this to the total tendency computed above
            zm['qtendtotal'] = zm['qtendtotal'] + zm['qtendsource'] + zm['qtendsink']
            zm['qtendtotal'].attrs['units'] = 'day/s'
            # compute difference of TEM AOA tendency and measured AOA tendency
            zm['qtenddiff'] = zm['qtendtotal'] - zm['AOATEND']
            zm['qtenddiff'].attrs['long_name'] = 'difference of qtendtotal and AOATEND'
            zm['qtenddiff'].attrs['units'] = 'day/s'
       
            if(histnum == 1):
                # get integrated AOA tendency components 
                q0 = zm['AOA'].isel(time=0, drop=True) # initial condition
                
                zm['AOATEND_int'] = integrate(q0, zm['AOATEND'])
                # for AOATEND, this method should exactly recover AOA...
                if(np.max(np.abs(zm['AOATEND_int'] - zm['AOA'])) > tol*np.max(zm['E90'])):
                    raise RuntimeError('AOATEND integration failed sanity check')
                zm['AOATEND_int'].attrs['units'] = 'day'
                
                zm['qtendetfd_int'] = integrate(q0, zm['qtendetfd'])
                zm['qtenddiff_int'].attrs['units'] = 'day'
                
                zm['qtendresvel_int'] = integrate(q0, zm['qtendresvel'])
                zm['qtendresvel_int'].attrs['units'] = 'day'
                 
                zm['qtenddiff_int'] = integrate(q0, zm['qtenddiff'])
                zm['qtenddiff_int'].attrs['units'] = 'day'

        if(qi == 2):
            # compute E90 tendency by parameterized tracer sources and sinks
            # E90 has no source
            zm['qtendsource'] = zm['E90j'] * 0
            # e-folding decay with timescale of 90 days
            tau = 90 * 24 * 60 * 60
            zm['qtendsink'] = -1/tau * zm['E90j']
            # add this to the total tendency computed above
            zm['qtendtotal'] = zm['qtendtotal'] + zm['qtendsource'] + zm['qtendsink']
            # compute difference of TEM E90 tendency and measured E90 tendency
            zm['qtenddiff'] = zm['qtendtotal'] - zm['E90TEND']
            zm['qtenddiff'].attrs['long_name'] = 'difference of qtendtotal and E90TEND'
            
            if(histnum == 1):
                # get integrated E90 tendency components
                q0 = zm['E90'].isel(time=0, drop=True) # initial condition
                
                zm['E90TEND_int'] = integrate(q0, zm['E90TEND'])
                # for E90TEND, this method should exactly recover E90...
                if(np.max(np.abs(zm['E90TEND_int'] - zm['E90'])) > tol*np.max(zm['E90'])):
                    raise RuntimeError('E90TEND integration failed sanity check')
                zm['E90TEND_int'].attrs['units'] = 'day'
                
                zm['qtendetfd_int'] = integrate(q0, zm['qtendetfd'])
                zm['qtenddiff_int'].attrs['units'] = 'day'
                
                zm['qtendresvel_int'] = integrate(q0, zm['qtendresvel'])
                zm['qtendresvel_int'].attrs['units'] = 'day'
                 
                zm['qtenddiff_int'] = integrate(q0, zm['qtenddiff'])
                zm['qtenddiff_int'].attrs['units'] = 'day'
            

    if(not dry):
        print('writing out...')
        zm = zm.drop_vars(orig_vars)
        zm.to_netcdf(outfile)
