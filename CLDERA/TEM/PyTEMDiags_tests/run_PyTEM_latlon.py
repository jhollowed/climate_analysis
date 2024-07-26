# Joe Hollowed
# University of Michigan 2023
#
# Code to perform test runs of the PyTEM package, as well as compare results
# to implementations by Christiane Jablonowski
#
# This script is the same as run_PyTEM.py, except that this version sends a raveled set of
# structured latlon data to the PyTEM utility, for testing purposes.

import pdb
import PyTEMDiags

import numpy as np
import xarray as xr
import climate_toolbox as ctb
import matplotlib.pyplot as plt

mapsdir = '/pscratch/sd/j/jhollo/PyTEMDiags/maps'
outdir = '/pscratch/sd/j/jhollo/PyTEMDiags/output'
nc = '/pscratch/sd/j/jhollo/E3SM/historical_data/TEM_test_data/v2.LR.WCYCL20TR.pmcpu.limvar.ens1.eam.h1.1998-11-23_TEM_VARIABLES_remap_180x360_aave.nc'

data = xr.open_dataset(nc)
data = data.thin({'lat':1, 'lon':3}) #sub-sample data...
data = data.drop_vars(('lon_bnds'))
lev  = data.lev
p0   = float(data['P0'].values)
L    = 45
 
data = data.sel({'lev':slice(0, 170)})
data, lat_weights = PyTEMDiags.tem_util.format_latlon_data(data)
lat_weights = None #tmp

overwrite_map = False

print('creating aveaging object...')
ua, va, ta, wap, ps, lat = data['U'], data['V'], data['T'],\
                           data['OMEGA'], data['PS'], data['lat']
p = ctb.compute_hybrid_pressure(ps, data['hyam'], data['hybm'], dims_like=ta, p0=p0)
tem = PyTEMDiags.TEMDiagnostics(ua, va, ta, wap, p, lat, p0=p0, L=L, lat_weights=lat_weights,
                                overwrite_map=overwrite_map, debug=True, 
                                grid_name='180x360', map_save_dest=mapsdir)
tem.to_netcdf(loc=outdir, include_attrs=True)

# ------------------------------------------------------------
# ---------------- analyze variables -------------------------

temdat = xr.open_dataset(tem.out_file).isel(time=0)
lat = temdat.lat
lev = temdat.lev
figdir = '/global/homes/j/jhollo/repos/climate_analysis/CLDERA/TEM/figs/individual_variable_sanity_check'

exit(0)

i=0
for var in list(temdat.keys()):
    if('ncol' in temdat[var].dims):
        continue
    i+=1
    vv = temdat[var].transpose('lev', 'lat')
    print('plotting {}...'.format(vv.name))
    plt.contourf(lat, lev, vv, levels=30)
    plt.gca().set_yscale('log')
    if(i==1): plt.gca().invert_yaxis()
    plt.xlabel('lat')
    plt.ylabel('lev')
    plt.title(vv.name)
    plt.savefig('{}/{}_{}.png'.format(figdir, tem.out_file.split('/')[-1].split('.nc')[0], vv.name))
    plt.show()
