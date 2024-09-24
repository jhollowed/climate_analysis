'''
Joe Hollowed 
University of Michigan 2024

This script performs vertical interpolaton on a single ensemble member of the 
cldera limvar dataset. Command line arguments are:

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
    index of the first history file to include in the interpolation
tmaxi : int
    index of the last history file to include in the interpolation
overwrite : int
    whether or not to overwrite the interpolated data if it already exists
    at the specified outpit file. 0 = False, 1 = True, 
    If False, then the script skips the interpolation for any file that 
    already exists
dry : int
    whether or not to do a "dry run" of the function, printing information
    about the execution but skipping the actual interpolation. 0 = False, 1 = True

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
from geocat.comp import interpolation as gcinterp

outdir = '/gpfs/cldera/data/e3sm/e3sm-cldera/CLDERA-historical/limvar_src_tag/tem_processed/interp_to_pressure_3D'
loc    = '/gpfs/cldera/data/e3sm/e3sm-cldera/CLDERA-historical/limvar_src_tag'

Nens      = int(sys.argv[1])
histnum   = int(sys.argv[2])
massMag   = int(sys.argv[3])
tmini     = int(sys.argv[4])
tmaxi     = int(sys.argv[5])
overwrite = bool(int(sys.argv[6]))
dry       = bool(int(sys.argv[7]))
print('args: {}'.format(sys.argv))

cfb = massMag == 0   # counterfactual flag

# ----------------------------------------------------------------

# get ensemble data
print('locating data...')
ensdirs = np.array(sorted(glob.glob('{}/*{}Tg.ens{}'.format(loc, massMag, Nens))))
mask    = np.array([not 'hud' in s for s in ensdirs])
ensdirs = ensdirs[mask]
ensdir  = ensdirs[0]
enshist = sorted(glob.glob('{}/archive/atm/hist/*h{}*'.format(ensdir, histnum)))[tmini:tmaxi+1]

cfdirs  = np.array(sorted(glob.glob('{}/*ens{}.cf'.format(loc, Nens))))
mask    = np.array([not 'hud' in s for s in cfdirs])
cfdirs  = cfdirs[mask]
cfdir   = cfdirs[0]
cfhist  = sorted(glob.glob('{}/archive/atm/hist/*h{}*'.format(cfdir, histnum)))[tmini:tmaxi+1]
if(cfb):
    ensdirs=cfdirs
    ensdir=cfdir
    enshist=cfhist

# do interpolation

print('-------- working on ENS {}, mass {}Tg, times {}-{}, --------'.format(
       Nens, massMag, tmini, tmaxi))

print('first file is: {}'.format(enshist[0].split('/')[-1]))
print('last file is: {}'.format(enshist[-1].split('/')[-1]))

for i in range(len(enshist)):
    
    fname = enshist[i].split('/')[-1]
    print('\n ----- working on {}'.format(fname))
    
    outfile = '{}/{}_pinterp.nc'.format(outdir, fname.split('.nc')[0])
    if(os.path.isfile(outfile) and not overwrite):
        print('output exists; skipping...')
        continue

    data = xr.open_dataset(enshist[i])
    
    # first select only the data in the region where hybm > 0
    shape =  data['U'].shape
    hybm_mask = data.hybm > 0
    # flip the bool of the lower-most level with hybm=0 to avoid edge nans 
    # in the temperature interpolation
    hybm_mask[len(hybm_mask) - list(hybm_mask.values)[::-1].index(0) - 1] = True
    data = data.isel(lev=(hybm_mask))
    print('shape after cut on hybm: {} --> {}'.format(shape, data['U'].shape))
    
    # get coords, etc
    p0        = data['P0']
    PS        = data['PS']
    lat       = data['lat']
    lev       = data['lev']
    hybm      = data['hybm']
    time      = data['time']
    time_bnds = data['time_bnds']
    attrs     = data.attrs
    g0 = 9.80665
    
    if(dry): exit(0)

    # get native data precision
    ppua, ppva, ppwap, ppaoa, ppe90, ppta, ppp0, pplat, pplev = \
    data['U'].dtype, data['V'].dtype, data['OMEGA'].dtype, data['AOA'].dtype, data['E90j'].dtype, \
    data['T'].dtype, data['P0'].dtype, data['lat'].dtype, data['lev'].dtype
    if(histnum == 0):
        ppgw1, ppgw2, ppgw3 = data['BUTGWSPEC'].dtype, data['UTGWORO'].dtype, data['UTGWSPEC'].dtype

    # ------ interpolate to isobars
    print('Generating interpolation...')
    interp_args = {'ps':data['PS'], 'hyam':data['hyam'], 'hybm':data['hybm'], 'p0':float(p0.values),
                   'new_levels':lev.values*100, 'method':'log', 'extrapolate':True, 
                   'variable':'other'}
    ua  = gcinterp.interp_hybrid_to_pressure(data['U'], **interp_args)
    va  = gcinterp.interp_hybrid_to_pressure(data['V'], **interp_args)
    wap = gcinterp.interp_hybrid_to_pressure(data['OMEGA'], **interp_args)
    aoa = gcinterp.interp_hybrid_to_pressure(data['AOA'], **interp_args)
    e90 = gcinterp.interp_hybrid_to_pressure(data['E90j'], **interp_args)
    ua['plev']  = ua.plev/100
    va['plev']  = va.plev/100
    wap['plev'] = wap.plev/100
    aoa['plev'] = aoa.plev/100
    e90['plev'] = e90.plev/100

    t_interp_args = copy.deepcopy(interp_args)
    t_interp_args['variable'] = 'temperature'
    t_interp_args['t_bot']    = data['T'].isel(lev=-1)
    t_interp_args['phi_sfc']  = data['Z3'].isel(lev=-1) * g0
    ta = gcinterp.interp_hybrid_to_pressure(data['T'], **t_interp_args)
    ta['plev'] = ta.plev/100
    ta = ta.drop_vars('lev')
        
    print('Computing interpolation...')
    ua, va, ta, wap, aoa, e90 = dask.compute(ua, va, ta, wap, aoa, e90)
    
    # get isobars
    plev = ua['plev']
    
    # now remove single datapoint with hybm=0
    ta   = ta.isel(plev=slice(1, len(plev)))
    ua   = ua.isel(plev=slice(1, len(plev)))
    va   = va.isel(plev=slice(1, len(plev)))
    wap  = wap.isel(plev=slice(1, len(plev)))
    aoa  = aoa.isel(plev=slice(1, len(plev)))
    e90  = e90.isel(plev=slice(1, len(plev)))
    
    # if monthly data, then also interpolate the gravity wave forcing 
    if(histnum == 0):
        gw1  = gcinterp.interp_hybrid_to_pressure(data['BUTGWSPEC'], **interp_args)
        gw2  = gcinterp.interp_hybrid_to_pressure(data['UTGWORO'], **interp_args)
        gw3  = gcinterp.interp_hybrid_to_pressure(data['UTGWSPEC'], **interp_args)
        gw1['plev'] = gw1.plev/100
        gw2['plev'] = gw2.plev/100
        gw3['plev'] = gw3.plev/100

        print('Computing GW interpolation for monthly data...')
        gw1, gw2, gw3 = dask.compute(gw1, gw2, gw3) 
    
        gw1 = gw1.isel(plev=slice(1, len(plev)))
        gw2 = gw2.isel(plev=slice(1, len(plev)))
        gw3 = gw3.isel(plev=slice(1, len(plev)))
    plev = plev.isel(plev=slice(1, len(plev)))
     
    # cast all variables to their original type (precision was likely promoted by dask.compute)
    print('casting variables to original data types...')
    ua   = ua.astype(ppua)
    va   = va.astype(ppva)
    wap  = wap.astype(ppwap)
    aoa  = aoa.astype(ppaoa)
    e90  = e90.astype(ppe90)
    ta   = ta.astype(ppta)
    p0   = p0.astype(ppp0)
    lat  = lat.astype(pplat)
    plev = plev.astype(pplev)
    varout = [ua, va, wap, aoa, e90, ta]
    if(histnum == 0):
        gw1 = gw1.astype(ppgw1)
        gw2 = gw2.astype(ppgw2)
        gw3 = gw3.astype(ppgw3)
        varout = varout + [gw1, gw2, gw3]
    
    # merge all interpolated data variables into Dataset
    print('merging variables...')
    data = xr.merge(varout)

    # sanity check for nans
    numnan = sum([int(np.isnan(data[dv]).sum()) for dv in data.data_vars])
    if(numnan > 0):
        raise RuntimeError('nans found in interpolated variables! Debug')

    # add other variables
    data['plev']      = plev
    data['lat']       = lat
    data['time']      = time
    data['time_bnds'] = time_bnds
    data['P0']        = p0
    data.attrs        = attrs
    data.attrs['history'] = '{}, interpolated to pure pressure levels on {}'.format(
                            attrs['history'], datetime.today().strftime('%Y-%m-%d'))

    # add attribute to the data which shall point to the original dataset, 
    # which will need to be accessed for all data with hybm = 0
    data.attrs['parent_dataset'] = enshist[i]
    
    # done, write out
    print('writing to {}'.format(outfile.split('/')[-1]))
    data.to_netcdf(outfile)
