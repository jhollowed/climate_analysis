import sys
import pdb
import glob
import numpy as np
import xarray as xr
from matplotlib import colors
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

TICK_FS = 11
LABEL_FS = 14
plt.rc('font', size=TICK_FS)          # controls default text sizes
plt.rc('axes', titlesize=LABEL_FS)     # fontsize of the axes title
plt.rc('axes', labelsize=LABEL_FS)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=TICK_FS)    # fontsize of the tick labels
plt.rc('ytick', labelsize=TICK_FS)    # fontsize of the tick labels
plt.rc('legend', fontsize=TICK_FS)    # legend fontsize
plt.rc('figure', titlesize=LABEL_FS)  # fontsize of the figure title

# --------------------------------

# --- year-round monthly data
pdir    = '/pscratch/sd/j/jhollo/E3SM/historical/CLDERA_2ndHistorical_1950-2014_monthly'
f_E90   = '{}/histoircal_h0_E90j_1850-2014.regrid.fv0.9x1.25_bilinear.nc'.format(pdir)
f_ST80  = '{}/histoircal_h0_ST80_25j_1850-2014.regrid.fv0.9x1.25_bilinear.nc'.format(pdir)
f_AOA   = '{}/histoircal_h0_AOA_1850-2014.regrid.fv0.9x1.25_bilinear.nc'.format(pdir)
f_PS    = '{}/histoircal_h0_PS_1850-2014.regrid.fv0.9x1.25_bilinear.nc'.format(pdir)
f_T     = '{}/histoircal_h0_T_1850-2014.regrid.fv0.9x1.25_bilinear.nc'.format(pdir)
zm_E90  = '{}/histoircal_h0_E90j_1850-2014.regrid.fv0.9x1.25_bilinear.zonalmean.nc'.format(pdir)
zm_ST80 = '{}/histoircal_h0_ST80_25j_1850-2014.regrid.fv0.9x1.25_bilinear.zonalmean.nc'.format(pdir)
zm_AOA  = '{}/histoircal_h0_AOA_1850-2014.regrid.fv0.9x1.25_bilinear.zonalmean.nc'.format(pdir)
zm_PS   = '{}/histoircal_h0_PS_1850-2014.regrid.fv0.9x1.25_bilinear.zonalmean.nc'.format(pdir)
zm_T    = '{}/histoircal_h0_T_1850-2014.regrid.fv0.9x1.25_bilinear.zonalmean.nc'.format(pdir)

try:
    print('reading zonal means...')
    dat  = xr.open_dataset(zm_E90)
    e90  = dat['E90j']
    st80 = xr.open_dataset(zm_ST80)['ST80_25j']
    aoa  = xr.open_dataset(zm_AOA)['AOA']
    ps   = xr.open_dataset(zm_PS)['PS']
    t   = xr.open_dataset(zm_T)['T']
except FileNotFoundError:
    print('zonal means not found; reading data, averaging...')
    
    runall=True
    if(sys.argv[1] is not None): runall=False

    if(sys.argv[1] == 'e90' or runall):
        print('working on E90...')
        dat_E90  = xr.open_dataset(f_E90)
        print('averaging...')
        e90  = dat_E90['E90j'].mean('lon') * 1e9     # E90 in ppb
        print('writing...')
        e90.to_netcdf(zm_E90)
        dat_E90.close()

    if(sys.argv[1] == 'st80' or runall):
        print('working on ST80...')
        dat_ST80 = xr.open_dataset(f_ST80)
        print('averaging...')
        st80 = dat_ST80['ST80_25j'].mean('lon') * 1e9 # ST80 in ppb 
        print('writing...')
        st80.to_netcdf(zm_ST80)
        dat_ST80.close()

    if(sys.argv[1] == 'aoa' or runall):
        print('working on AOA...')
        dat_AOA  = xr.open_dataset(f_AOA)
        print('averaging...')
        aoa  = dat_AOA['AOA'].mean('lon') * 1/365    # AOA in years
        print('writing...')
        aoa.to_netcdf(zm_AOA)
        dat_AOA.close()
    
    if(sys.argv[1] == 'ps' or runall):
        print('working on PS...')
        dat_PS  = xr.open_dataset(f_PS)
        print('averaging...')
        aoa  = dat_PS['PS'].mean('lon')
        print('writing...')
        aoa.to_netcdf(zm_PS)
        dat_PS.close()
        
    if(sys.argv[1] == 't' or runall):
        print('working on T...')
        dat_T  = xr.open_dataset(f_T)
        print('averaging...')
        aoa  = dat_T['T'].mean('lon')
        print('writing...')
        aoa.to_netcdf(zm_T)
        dat_T.close()

# --- take time mean
e90  = e90.mean('time')
st80 = st80.mean('time')
aoa  = aoa.mean('time')


lat = dat['lat']
lev = dat['lev']
LAT, LEV = np.meshgrid(lat, lev)

fig = plt.figure(figsize=(14, 5))
ax1 = fig.add_subplot(131)
ax2 = fig.add_subplot(132)
ax3 = fig.add_subplot(133)


# ---------- E90
print('plotting E90...')
levels = [0.1, 1, 10, 20, 30, 40, 50, 60, 70, 
          80, 90, 100, 110, 120, 130, 140, 150]
lablevels = [1, 10, 50, 100, 110, 130, 150]
divnorm=colors.TwoSlopeNorm(vmin=0, vcenter=50, vmax=150)
cf = ax1.contourf(LAT, LEV, e90, levels=levels, cmap=plt.cm.RdYlBu_r, norm=divnorm, extend='both')
ax1.contour(LAT, LEV, e90, levels=levels, colors=['k'], linewidths=[0.4])          # draw contours
ax1.contour(LAT, LEV, e90, levels=[90], colors=['k'], linewidths=[2.5])            # bold 90 ppb
cfl = ax1.contour(LAT, LEV, e90, levels=lablevels, colors=['k'], linewidths=[0.4]) # contours to label
cflm01 = ax1.contour(LAT, LEV, e90, levels=[0.1], colors=['k'], linewidths=[0.4])  # label 0.1 ppb

ax1.clabel(cfl, inline=1, fontsize=TICK_FS, fmt='%1.0f')
ax1.clabel(cflm01, inline=1, fontsize=TICK_FS, fmt='%1.1f')
cb  = fig.colorbar(cf, ax=ax1, location='top')
cbt = [0, 50, 100, 150]
cb.set_ticks(cbt)
cb.set_ticklabels([str(t) for t in cbt])
cb.set_label('E90 [ppb]')

ax1.set_ylim([30, 900])
ax1.set_yscale('log')
ax1.invert_yaxis()
ax1.set_ylabel('lev [hPa]')
ax1.set_xlabel('lat [deg]')
ax1.yaxis.set_major_formatter(ScalarFormatter())
ax1.yaxis.set_minor_formatter(ScalarFormatter())
for label in ax1.get_yticklabels(minor=True)[::2]:
    label.set_visible(False)


# ---------- ST80
print('plotting ST80...')
levels=[20, 40, 60, 80, 100, 120, 140, 160, 180, 199]
lablevels=[20, 60, 100, 140, 180, 199]
cf = ax2.contourf(LAT, LEV, st80, levels=levels, cmap=plt.cm.YlOrRd, extend='both')
ax2.contour(LAT, LEV, st80, levels=levels, colors=['k'], linewidths=[0.4])          # draw contours
cfl = ax2.contour(LAT, LEV, st80, levels=lablevels, colors=['k'], linewidths=[0.4]) # contours to label

ax2.clabel(cfl, inline=1, fontsize=TICK_FS, fmt='%1.0f')
cb  = fig.colorbar(cf, ax=ax2, location='top')
cb.set_label('ST80_25 [ppb]')

ax2.set_ylim([50, 200])
ax2.set_yscale('log')
ax2.invert_yaxis()
ax2.set_ylabel('lev [hPa]')
ax2.set_xlabel('lat [deg]')
ax1.yaxis.set_major_formatter(ScalarFormatter())
ax1.yaxis.set_minor_formatter(ScalarFormatter())


# ---------- AOA
print('plotting AOA...')
levels=[0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5]
lablevels=[1, 2, 3, 4, 5]
cf = ax3.contourf(LAT, LEV, aoa, levels=levels, cmap=plt.cm.YlGnBu_r, extend='both')
ax3.contour(LAT, LEV, aoa, levels=levels, colors=['k'], linewidths=[0.4])          # draw contours
cfl = ax3.contour(LAT, LEV, aoa, levels=lablevels, colors=['k'], linewidths=[0.4]) # contours to label

ax3.clabel(cfl, inline=1, fontsize=TICK_FS, fmt='%1.0f')
cb  = fig.colorbar(cf, ax=ax3, location='top')
cb.set_label('AOA [years]')

ax3.set_ylim([1, 100])
ax3.set_yscale('log')
ax3.yaxis.set_label_position("right")
ax3.yaxis.tick_right()
ax3.invert_yaxis()
ax3.set_ylabel('lev [hPa]')
ax3.set_xlabel('lat [deg]')
ax3.yaxis.set_major_formatter(ScalarFormatter())
ax3.yaxis.set_minor_formatter(ScalarFormatter())
for label in ax3.get_yticklabels(minor=True)[1::2]:
    label.set_visible(False)

plt.show()
