'''
Joe Hollowed 2024

This script takes the output of a previously generated set of vertically interpolated
data (from limvar_zinterp.py) and casts all of the variables as their data types from 
the native dataset. This is because in an old version of the limvar_zinterp.py script, 
I negelected to handle the fact that dask.compute() was implicity promoting many variables
from float32 to float64. That script has been updated to correctly handle this, but rather 
than re-running loimvar_zinterp.py, we can instead run the present script, which skips
the interpolation and does only the casting.
Command line arguments are:

Usage
-----
python ./limvar_zinterp.py [Nens] [histnum] [massMag] [tmini] [tmaxi] [dry]

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
dry       = bool(int(sys.argv[6]))
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

cfdirs  = np.array(sorted(glob.glob('{}/*ens{}*cf'.format(loc, Nens))))
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
    if(not os.path.isfile(outfile)):
        print('output does not exist; skipping...')
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
    p0        = float(data['P0'].values)
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

    # done with native data, open interpolation output
    data = xr.open_dataset(outfile)

    # cast all variables to their original type (precision was likely promoted by dask.compute)
    print('casting variables to original data types...')
    data['U']   = data['U'].astype(ppua)
    data['V']   = data['V'].astype(ppva)
    data['OMEGA']  = data['OMEGA'].astype(ppwap)
    data['AOA']  = data['AOA'].astype(ppaoa)
    data['E90j']  = data['E90j'].astype(ppe90)
    data['T']   = data['T'].astype(ppta)
    data['P0']   = data['P0'].astype(ppp0)
    data['lat']  = data['lat'].astype(pplat)
    data['plev'] = data['plev'].astype(pplev)
    if(histnum == 0):
        data['BUTGWSPEC'] = data['BUTGWSPEC'].astype(ppgw1)
        data['UTGWORO'] = data['UTGWORO'].astype(ppgw2)
        data['UTGWSPEC'] = data['UTGWSPEC'].astype(ppgw3)

    # done, write out
    print('writing to {}'.format(outfile.split('/')[-1]))
    data.to_netcdf(outfile)
