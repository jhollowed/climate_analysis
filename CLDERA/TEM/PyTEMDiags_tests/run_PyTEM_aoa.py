# Joe Hollowed
# University of Michigan 2023
#
# Code to perform test runs of the PyTEM package, as well as compare results
# to implementations by Tom Ehrmann and Christiane Jablonowski

import sys
import pdb
import dask
import PyTEMDiags
import numpy as np
import xarray as xr
import climate_toolbox as ctb
import matplotlib.pyplot as plt
from geocat.comp import interpolation as gcinterp

mapsdir = '/pscratch/sd/j/jhollo/PyTEMDiags/maps'
outdir = '/pscratch/sd/j/jhollo/PyTEMDiags/output'
# ensemble member from limvar_src_tag with daily data
nc = '/pscratch/sd/j/jhollo/E3SM/historical_data/cldera_limvar/limvar_src_tag/v2.LR.WCYCL20TR.pmcpu.ctools.lv.3tagso4.10Tg.ens1.eam.h1.1991-06-01-00000.nc'

data = xr.open_dataset(nc)
p0 = float(data['P0'].values)
lat = data['lat']
lev = data['lev']
PS = data['PS']
L = 90
g0 = 9.80665

# ------ interpolate to isobars
print('Generating interpolation...')
interp_args = {'ps':data['PS'], 'hyam':data['hyam'], 'hybm':data['hybm'], 'p0':p0, 
               'new_levels':lev.values*100, 'method':'log', 'extrapolate':True, 'variable':'other'}

ua  = gcinterp.interp_hybrid_to_pressure(data['U'], **interp_args)
va  = gcinterp.interp_hybrid_to_pressure(data['V'], **interp_args)
wap = gcinterp.interp_hybrid_to_pressure(data['OMEGA'], **interp_args)
aoa = gcinterp.interp_hybrid_to_pressure(data['AOA'], **interp_args)
e90 = gcinterp.interp_hybrid_to_pressure(data['E90j'], **interp_args)
ua['plev'] = ua.plev/100
va['plev'] = va.plev/100
wap['plev'] = wap.plev/100
aoa['plev'] = aoa.plev/100
e90['plev'] = e90.plev/100

interp_args['variable'] = 'temperature'
interp_args['t_bot']    = data['T'].isel(lev=-1)
interp_args['phi_sfc']  = data['Z3'].isel(lev=-1) * g0
ta = gcinterp.interp_hybrid_to_pressure(data['T'], **interp_args)
ta['plev'] = ta.plev/100

print('Computing interpolation...')
ua, va, ta, wap, aoa, e90 = dask.compute(ua, va, ta, wap, aoa, e90)
plev = ua['plev']

# ------ compute TEM
print('Computing TEM...')
tem = PyTEMDiags.TEMDiagnostics(ua, va, ta, wap, plev, lat, q=[aoa, e90], p0=p0, L=L, 
                                overwrite_map=False, debug_level=1, 
                                grid_name='ne30pg2', map_save_dest=mapsdir)

print('writing out tem...')
tem.to_netcdf(loc=outdir, include_attrs=False)
tem.q_to_netcdf(loc=outdir, include_attrs=True)

# test 
#print('generating test figure...')
#ps = tem.ZM.sph_zonal_mean(PS.transpose())
#ax = fig.add_subplot()
#pdb.set_trace()
#ax.contourf(tem.lat, tem.plev, tem.ub.mean('time').T, levels=20, cmap='rainbow')
#ax.plot(tem.lat, ps.mean('time')/100, '-k', lw=2)
#ax.set_yscale('log')
#ax.invert_yaxis()
#plt.show()


# ------------------------------------------------------------
# ---------------- analyze variables -------------------------

temdat = xr.open_dataset(tem.out_file).mean('time')
lat = temdat.lat
lev = temdat.plev
figdir = '/global/homes/j/jhollo/repos/climate_analysis/CLDERA/TEM/figs/individual_variable_sanity_check'

i=0
for var in list(temdat.keys()):
    if('ncol' in temdat[var].dims):
        continue
    i+=1
    vv = temdat[var].transpose('plev', 'lat')
    print('plotting {}...'.format(vv.name))
    plt.contourf(lat, lev, vv, levels=30)
    plt.gca().set_yscale('log')
    plt.gca().invert_yaxis()
    plt.xlabel('lat')
    plt.ylabel('lev')
    plt.title(vv.name)
    plt.savefig('{}/{}_{}.png'.format(figdir, tem.out_file.split('/')[-1].split('.nc')[0], vv.name))

for qi in range(len(tem.q_out_file)):
    temdat = xr.open_dataset(tem.q_out_file[qi]).mean('time')
    tracer = tem.q_out_file[qi].split('.nc')[0].split('-')[-1]
    i=0
    for var in list(temdat.keys()):
        if('ncol' in temdat[var].dims):
            continue
        i+=1
        vv = temdat[var].transpose('plev', 'lat')
        print('plotting {}...'.format(vv.name))
        plt.contourf(lat, lev, vv, levels=30)
        plt.gca().set_yscale('log')
        plt.gca().invert_yaxis()
        plt.xlabel('lat')
        plt.ylabel('lev')
        plt.title(vv.name)
        plt.savefig('{}/{}_{}_{}.png'.format(figdir, tem.out_file.split('/')[-1].split('.nc')[0], 
                                             vv.name, tracer))
