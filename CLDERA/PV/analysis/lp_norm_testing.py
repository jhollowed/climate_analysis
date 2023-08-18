import pdb
import numpy as np
import xarray as xr
import climate_toolbox as ctb
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable

# ---------- read latlon
loc = '/pscratch/sd/j/jhollo/E3SM/E3SMv2_cases/pv_cases/highvar_PVPTtest_ens/injection__ens01/run'
f = '{}/injection__ens01.eam.h0.0001-01-01-00000.regrid.91x180_aave.nc'.format(loc)
dat = xr.open_dataset(f)

print('reading latlon...')
lat, lev, lon = dat['lat'], dat['lev'], dat['lon']
PS, hyam, hybm, hyai, hybi = dat['PS'], dat['hyam'], dat['hybm'], dat['hyai'], dat['hybi']
T, PT = dat['T'], dat['PT']
LAT, LEV = np.meshgrid(lat, lev)

# ---------- read native
loc_n = '/pscratch/sd/j/jhollo/E3SM/E3SMv2_cases/pv_cases/highvar_PVPTtest_ens/injection__ens01/run'
f_n = '{}/injection__ens01.eam.h0.0001-01-01-00000.nc'.format(loc_n)
dat_n = xr.open_dataset(f_n)

print('reading native...')
ncol = dat_n['ncol']
PS_n = dat_n['PS']
T_n, PT_n = dat_n['T'], dat_n['PT']
hw = dat_n['area']

# ---------- test various use cases of lp_norm 
# (which effectively tests use cases of the ctb weighted-mean functions)
# lp norm returned should be zero for all cases

print('\n---- 2D latlon, single time')
res = ctb.lp_norm(T[0].isel(lev=0), T[0].isel(lev=0), p=4)
print('NORM: {}'.format(res.values))

print('\n---- 2D native, single time')
res = ctb.lp_norm(T_n[0].isel(lev=0), T_n[0].isel(lev=0), p=4, horz_weights=hw)
print('NORM: {}'.format(res.values))

print('\n---- 2D latlon, all time')
res = ctb.lp_norm(T.isel(lev=0), T.isel(lev=0), p=4)
print('NORM: {}'.format(res.sum('time').values))

print('\n---- 2D native, all time')
res = ctb.lp_norm(T_n.isel(lev=0), T_n.isel(lev=0), p=4, horz_weights=hw)
print('NORM: {}'.format(res.sum('time').values))


print('\n---- 3D latlon, single time')
res = ctb.lp_norm(T[0], T[0], p=4, ps=PS[0], hyai=hyai, hybi=hybi)
print('NORM: {}'.format(res.values))

print('\n---- 3D native, single time')
res = ctb.lp_norm(T_n[0], T_n[0], p=4, horz_weights=hw, ps=PS_n[0], hyai=hyai, hybi=hybi)
print('NORM: {}'.format(res.values))

print('\n---- 3D latlon, all time')
res = ctb.lp_norm(T, T, p=4, ps=PS, hyai=hyai, hybi=hybi)
print('NORM: {}'.format(res.sum('time').values))

print('\n---- 3D native, all time')
res = ctb.lp_norm(T_n, T_n, p=4, horz_weights=hw, ps=PS_n, hyai=hyai, hybi=hybi)
print('NORM: {}'.format(res.sum('time').values))