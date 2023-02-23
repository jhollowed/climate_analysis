import numpy as np
import matplotlib.pyplot as plt
import pdb
import glob
import xarray as xr
import cartopy.crs as ccrs
from climate_artist import horizontal_slice as plthor
from climate_artist import horizontal_slice as pltvert

# ==================== read data ==========

loc = '/global/cfs/cdirs/m4014/data/HSW/post/release_011423/ens_members_latlon'
out = '/global/homes/j/jhollo/repos/climate_analysis/CLDERA/SAI/analysis/tmp'

ens = [xr.open_dataset(glob.glob('{}/ens0{}/AOD.h2*'.format(loc, i+1))[0])['AOD'] for i in range(5)]

ensmean = np.zeros(ens[0].shape)
ensstd = np.zeros(ens[0].shape)

print('computing ensemble mean')
for i in range(5):
    ensmean = ensmean + ens[i]
ensmean = ensmean / 5

print('computing ensemble std')
for i in range(5):
    ensstd = ensstd + (ens[i] - ensmean)**2
ensstd = ensstd / 4
ensstd = np.sqrt(ensstd)

nlat = len(ensmean['lat'])
nlon = len(ensmean['lon'])
nt = len(ensmean['time'])
lat = ensmean['lat']
lon = ensmean['lon']
time = ensmean['time']

print('writing out')
ensmean.to_netcdf('{}/ensmean.nc'.format(out))
ensstd.to_netcdf('{}/ensstd.nc'.format(out))

# ==================== plot ===========

tstart = 91
clev = np.linspace(0, 0.5, 6)
cm = plt.cm.jet
pltargs = {'levels':clev, 'cmap':cm, 'extend':'both'}
LON, LAT = np.meshgrid(lon, lat)
for j in range(nt-tstart-1):
    
    print('creating figure...')
    fig = plt.figure(figsize=(10, 14))
    ax1 = fig.add_subplot(321, projection=ccrs.PlateCarree())
    ax2 = fig.add_subplot(322, projection=ccrs.PlateCarree())
    ax3 = fig.add_subplot(323, projection=ccrs.PlateCarree())
    ax4 = fig.add_subplot(324, projection=ccrs.PlateCarree())
    ax5 = fig.add_subplot(325, projection=ccrs.PlateCarree())
    ax6 = fig.add_subplot(326, projection=ccrs.PlateCarree())
    ax = [ax1, ax2, ax3, ax4, ax5, ax6]

    print('plotting time {}...'.format(j+1))
    for i in range(5):
        print('    plotting ens {}...'.format(i+1))
        var_dict = [{'var':ens[i][j], 'plotType':'contourf', 'plotArgs':pltargs, 'colorFormatter':None}]
        plthor(lon, lat, var_dict, ax=ax[i], annotation='ens0{}'.format(i+1))
        #ax[i].contourf(LON, LAT, ens[i][tstart+j], levels=clev, cmap=cm, extend='both')
    plt.show()
