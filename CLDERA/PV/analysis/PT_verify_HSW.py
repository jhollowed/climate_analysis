import pdb
import numpy as np
import xarray as xr
import climate_toolbox as ctb
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable

loc = '/pscratch/sd/j/jhollo/E3SM/E3SMv2_cases/pv_cases/highvar_PVPTtest_ens/injection__ens01/run'
tmp = '{}/PROCESSED'.format(loc)
f = '{}/injection__ens01.eam.h0.0001-01-01-00000.regrid.91x180_aave.nc'.format(loc)
dat = xr.open_dataset(f)

print('reading latlon...')
lat, lev, lon = dat['lat'], dat['lev'], dat['lon']
PS, hyam, hybm, hyai, hybi = dat['PS'], dat['hyam'], dat['hybm'], dat['hyai'], dat['hybi']
T, PT = dat['T'], dat['PT']
LAT, LEV = np.meshgrid(lat, lev)

print('computing theta...')
theta_tmp = '{}/theta.nc'.format(tmp)
force_write = False
try:
    if(force_write): raise FileNotFoundError
    theta = xr.open_dataset(theta_tmp)['THETA']
    print('    read')
except FileNotFoundError:
    theta = ctb.compute_theta(T, PS, hyam, hybm)
    theta.to_netcdf(theta_tmp)
    print('    writing')

print('taking means...')
Tzm_tmp = '{}/Tzm.nc'.format(tmp)
PTzm_tmp = '{}/PTzm.nc'.format(tmp)
thetazm_tmp = '{}/thetazm.nc'.format(tmp)
force_write = False
try:
    if(force_write): raise FileNotFoundError
    Tzm = xr.open_dataset(Tzm_tmp)['T']
    PTzm = xr.open_dataset(PTzm_tmp)['PT']
    thetazm = xr.open_dataset(thetazm_tmp)['THETA']
    print('    read')
except FileNotFoundError:
    Tzm = T.mean(['lon', 'time'])
    Tzm.to_netcdf(Tzm_tmp)
    PTzm = PT.mean(['lon', 'time'])
    PTzm.to_netcdf(PTzm_tmp)
    thetazm = theta.mean(['lon', 'time'])
    thetazm.to_netcdf(thetazm_tmp)
    print('    writing')

print('plotting...')
fig = plt.figure(figsize=(10, 5))
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)
pt_lev = np.arange(250, 800, 25)

im = ax1.contourf(LAT, LEV, Tzm, levels=20, cmap=plt.cm.plasma)
cn = ax1.contour(LAT, LEV, thetazm, levels=pt_lev, colors='w')
ax1.clabel(cn)
ax1.set_ylim([10, 1000])
ax1.invert_yaxis()
ax1.set_yscale('log')
ax1.set_title('computed PT [K] (white contours) \n30-day zonal-mean')
ax1.set_xlabel('lat  [deg]')
ax1.set_ylabel('lev  [hPa]')
ax1.yaxis.set_major_formatter(ScalarFormatter())
divider = make_axes_locatable(ax1)
cax = divider.append_axes('right', size='5%', pad=0.05)
cb = fig.colorbar(im, cax=cax, orientation='vertical')
cb.set_label('T  [K]')

im = ax2.contourf(LAT, LEV, Tzm, levels=20, cmap=plt.cm.plasma)
cn = ax2.contour(LAT, LEV, PTzm, levels=pt_lev, colors='w')
ax2.clabel(cn)
ax2.set_ylim([10, 1000])
ax2.invert_yaxis()
ax2.set_yscale('log')
ax2.set_title('diagnosed PT [K] (white contours) \n30-day zonal-mean')
ax2.set_xlabel('lat  [hPa]')
ax2.yaxis.set_major_formatter(ScalarFormatter())
divider = make_axes_locatable(ax2)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')
cb.set_label('T  [K]')

plt.tight_layout()
plt.savefig('PT_verify.png', dpi=200)
plt.show()
