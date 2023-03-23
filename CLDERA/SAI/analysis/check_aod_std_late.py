import numpy as np
import matplotlib.pyplot as plt
import pdb
import glob
import xarray as xr
import cartopy.crs as ccrs
from climate_artist import horizontal_slice as plthor
from climate_artist import horizontal_slice as pltvert

# ==================== read data ==========

# --- pm
#loc = '/global/cfs/cdirs/m4014/data/HSW/post/release_011423/ens_members_latlon'
#out = '/global/homes/j/jhollo/repos/climate_analysis/CLDERA/SAI/analysis/tmp'
# --- miniroomba
loc = '/Users/joe/tmp'
out = '/Users/joe/tmp'
figs = '/Users/joe/repos/climate_analysis/CLDERA/SAI/analysis/figs/aod_check_late'

ens = [xr.open_dataset(glob.glob('{}/ens0{}/AOD.h2*'.format(loc, i+1))[0])['AOD'] for i in range(5)]

# zero out almost-zeros
#for i in range(5):
#    pdb.set_trace()

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

tstart = 135
cm = plt.cm.rainbow
clev = np.linspace(0, 0.07, 6) + 1e-4
pltargs = {'levels':clev, 'cmap':cm, 'extend':'both'}

clev_c = np.linspace(0.01, 0.06, 3) + 1e-4
pltargs_c = {'levels':clev_c, 'colors':None, 'linewidths':1.25}

clev_mean = np.linspace(0, 0.07, 7) + 1e-4
clev_mean = clev_mean*1000000 #tmp
cm_mean = plt.cm.Greys
pltargs_mean = {'levels':clev_mean, 'cmap':cm_mean, 'extend':'both', 'alpha':1}

clev_std = np.linspace(0, 0.01, 6) + 1e-4
clev_std = clev_std * 10 #tmp
cm_std = plt.cm.rainbow
pltargs_std = {'levels':clev_std, 'cmap':cm_std, 'extend':'both', 'alpha':1}

LON, LAT = np.meshgrid(lon, lat)
colors = ['r', 'b', 'g', 'm', 'c']

fac = 5
clev = clev * fac
clev_c = clev_c * fac

for j in range(nt-tstart):
    
    print('creating figure...')
    fig = plt.figure(figsize=(10, 5))
    ax1 = fig.add_subplot(221, projection=ccrs.PlateCarree())
    ax2 = fig.add_subplot(222, projection=ccrs.PlateCarree())
    ax3 = fig.add_subplot(223, projection=ccrs.PlateCarree())
    ax4 = fig.add_subplot(224, projection=ccrs.PlateCarree())
    #ax5 = fig.add_subplot(325, projection=ccrs.PlateCarree())
    #ax6 = fig.add_subplot(326, projection=ccrs.PlateCarree())
    ax = [ax1, ax2, ax3, ax4]
    fig2 = plt.figure(figsize=(12, 3))
    ax6 = fig2.add_subplot(121, projection=ccrs.PlateCarree())
    ax7 = fig2.add_subplot(122, projection=ccrs.PlateCarree())
    fig.suptitle('Day {:.0f}'.format((tstart+j) * 2), fontsize=14)
    fig2.suptitle('Day {:.0f}'.format((tstart+j) * 2), fontsize=14)
        
    var_dict = [{'var':ensmean[j+tstart], 'plotType':'contourf', 
                 'plotArgs':pltargs_mean, 'colorFormatter':None}]
    plthor(lon, lat, var_dict, ax=ax6, annotation=None, include_contours=False, 
           coastlinesArgs={'lw':0.5, 'alpha':0.75})

    var_dict = [{'var':ensstd[j+tstart], 'plotType':'contourf', 
        'plotArgs':pltargs_std}]
    plthor(lon, lat, var_dict, ax=ax7, annotation=None, include_contours=True, 
           coastlinesArgs={'lw':0.5, 'alpha':0.75})

    print('plotting time {}...'.format(j+1))
    for i in range(5):
        if(i<4):
            print('    plotting ens {}...'.format(i+1))
            var_dict = [{'var':ens[i][j+tstart], 
                         'plotType':'contourf', 'plotArgs':pltargs, 'colorFormatter':None}]
            plthor(lon, lat, var_dict, ax=ax[i], annotation='ens0{}'.format(i+1), 
                   include_contours=False)
        for k in range(len(clev_c)):
            pltargs_c['colors']=colors[i]
            pltargs_c['alpha']=(k+1)/len(clev_c)
            pltargs_c['levels']=[clev_c[k]]
            var_dict = [{'var':ens[i][j+tstart], 'plotType':'contour', 
                         'plotArgs':pltargs_c, 'colorFormatter':None}]
            plthor(lon, lat, var_dict, ax=ax6, annotation=None, include_contours=False, 
                   coastlines=False)
        ax6.plot([0,0], [0,0], color=colors[i], label='ens0{}'.format(i+1))
        #ax6.set_extent((-180, 180, -40, 60))
        #ax7.set_extent((-180, 180, -40, 60))
        #ax6.set_aspect(1.5)
        #ax7.set_aspect(1.65)
    ax6.legend(loc='lower left', ncol=2)
    fig.tight_layout()
    fig2.tight_layout()
    plt.show()
    fig.savefig('{}/ens{:03d}.png'.format(figs, j))
    fig2.savefig('{}/std{:03d}.png'.format(figs, j))
