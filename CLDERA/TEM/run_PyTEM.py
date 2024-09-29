# Joe Hollowed
# University of Michigan 2023
#
# Code to perform test runs of the PyTEM package, as well as compare results
# to implementations by Tom Ehrmann and Christiane Jablonowski

import pdb
import PyTEMDiags
import xarray as xr
import climate_toolbox as ctb
import matplotlib.pyplot as plt

mapsdir = '/pscratch/sd/j/jhollo/PyTEMDiags/maps'
outdir = '/pscratch/sd/j/jhollo/PyTEMDiags/output'
nc = '/pscratch/sd/j/jhollo/E3SM/historical_data/TEM_test_data/v2.LR.WCYCL20TR.pmcpu.limvar.ens1.eam.h1.1998-11-23_TEM_VARIABLES.nc'
data = xr.open_dataset(nc)
p0 = float(data['P0'].values)
do_p_var = True
do_p_coord = False
overwrite_map = False
L = 150 # this value seems optimal for ne30pg2, according to sensitivty test
    
data=data.sel({'lev':slice(0, 170)})

# -----------------------------------------

if(do_p_var):
    ua, va, ta, wap, ps, lat = data['U'], data['V'], data['T'],\
                               data['OMEGA'], data['PS'], data['lat']
    p = ctb.compute_hybrid_pressure(ps, data['hyam'], data['hybm'], dims_like=ta, p0=p0)
    tem = PyTEMDiags.TEMDiagnostics(ua, va, ta, wap, p, lat, p0=p0, L=L, 
                                    overwrite_map=overwrite_map, debug=True, 
                                    grid_name='ne30pg2', map_save_dest=mapsdir)
    tem.to_netcdf(loc=outdir, include_attrs=True)

# -----------------------------------------

if(do_p_coord):
    ua, va, ta, wap, ps, lat = data['U'], data['V'], data['T'],\
                               data['OMEGA'], data['PS'], data['lat']
    p = ctb.compute_hybrid_pressure(ps, data['hyam'], data['hybm'], dims_like=ta, p0=p0)
    tem = PyTEMDiags.TEMDiagnostics(ua, va, ta, wap, p, lat, p0=p0, L=L, 
                                    overwrite_map=overwrite_map, debug=True, 
                                    grid_name='ne30pg2', map_save_dest=mapsdir)
    tem.to_netcdf(loc=outdir, include_attrs=True)

# ------------------------------------------------------------
# ---------------- analyze variables -------------------------

temdat = xr.open_dataset(tem.out_file).isel(time=0)
lat = temdat.lat
lev = temdat.lev
figdir = '/global/homes/j/jhollo/repos/climate_analysis/CLDERA/TEM/figs/individual_variable_sanity_check'

pdb.set_trace()

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
