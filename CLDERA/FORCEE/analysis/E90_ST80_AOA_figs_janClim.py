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

# --- jan instant. data
#pdir = '/global/homes/j/jhollo/repos/climate_analysis/CLDERA/FORCEE'
#f = '{}/fullvar_2ndHistorical_1850-2014/yearly_ic_avg.regrid._bilinear.nc'.format(pdir)
# --- year-round monthly data
pdir   = '/pscratch/sd/j/jhollo/E3SM/CLDERA_2ndHistorical_1950-2014_monthly'
f_E90  = '{}/histoircal_h0_E90j_1850-2014.regrid._bilinear.nc'.format(pdir)
f_ST80 = '{}/histoircal_h0_ST80_25j_1850-2014.regrid._bilinear.nc'.format(pdir)
f_AOA  = '{}/histoircal_h0_AOA_1850-2014.regrid._bilinear.nc'.format(pdir)
dat = xr.open_dataset(f)
dat = dat.sel(time=dat['time'][0]) # remove (length-0) time dimension

print('reading data, taking zonal means...')
e90  = dat['E90j'].mean('lon') * 1e9     # E90 in ppb
st80 = dat['ST80_25j'].mean('lon') * 1e9 # ST80 in ppb
aoa  = dat['AOA'].mean('lon') * 1/365    # AOA in years

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
ax2.set_xlabel('lat [deg]')
for label in ax2.get_yticklabels():
    label.set_visible(False)
for label in ax2.get_yticklabels(minor=True):
    label.set_visible(False)


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
