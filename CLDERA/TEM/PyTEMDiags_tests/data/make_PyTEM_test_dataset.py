
import pdb
import numpy as np
import xarray as xr

f = '/pscratch/sd/j/jhollo/E3SM/historical_data/TEM_test_data/v2.LR.WCYCL20TR.pmcpu.limvar.ens1.eam.h1.1998-11-23_TEM_VARIABLES.nc'
dat = xr.open_dataset(f)

var = dat['U'].mean('time')
var.attrs = {}
var.name = 'X'
var.values = np.ones(var.shape)

lat = dat['lat']
lon = dat['lon']

testdat = xr.Dataset({'X':var, 'lat':lat, 'lon':lon})
testdat.to_netcdf('testdat.nc')
