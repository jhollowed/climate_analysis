'''
Joe Hollowed
University of Michigan, March 2023

Script to plot spread between ensemble members as 1D time series of zonally-averaged variables
at specified (lat, lev) positions
'''

import pdb
import glob
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# point to ensemble location
loc = '/project/projectdirs/m4014/data/HSW/outputs/release_030123/pert_ens'
out = './out'  # output directory for figs and any temporary data
var = 'U'     # variable to plot
N = 8         # number of ensemble members (doesn't need to be hardcoded, could be read from 'loc')

# som hardcoded presets for T and U; other variables could be added as well
if(var == 'T'):
    unit = 'K'
    color = 'r'
if(var == 'U'):
    unit = 'm/s'
    color = 'b'

# get variable data for each ensemble member
# note that we assume the presence in the data directory of remapped lat-lon data files,
# which contain the substring "regrid". We also specify "h1" history files here
ensd = [xr.open_dataset(glob.glob('{}/pert0{}/*h1*regrid*'.format(loc, i+1))[0]) for i in range(N)]
ens = [ensd[i][var] for i in range(N)]

# Next, compute ensemble mean, ensemble std
# This is wrapped in a try;except which first tries to read these stats from file
# On the first run of this script, the try will fail with FileNotFound, and we compute
#  the mean, std from scratch, and then save the result
try:
    ensmean = xr.open_dataset('{}/pertmean_{}.nc'.format(out, var))[var]
    ensstd = xr.open_dataset('{}/pertstd_{}.nc'.format(out, var))[var]
    print('read')
except FileNotFoundError:
    ensmean = ens[0].copy(deep=True, data=np.zeros(ens[0].shape))
    ensstd = ens[0].copy(deep=True, data=np.zeros(ens[0].shape))
    print('computing ensemble mean')
    # sum each ensemble member, divide by N
    for i in range(N):
        ensmean = ensmean + ens[i]
    ensmean = ensmean / N

    print('computing ensemble std')
    # sum the squared deviation of each ensemble member, normalize by (N-1), sqrt
    for i in range(N):
        ensstd = ensstd + (ens[i] - ensmean)**2
    ensstd = ensstd / (N-1)
    ensstd = np.sqrt(ensstd)

    # write to netcdf to avoid computing on next execution of this script
    print('writing out')
    ensmean.to_netcdf('{}/pertmean_{}.nc'.format(out, var))
    ensstd.to_netcdf('{}/pertstd_{}.nc'.format(out, var))

# get data size and coords
nlat = len(ens[0]['lat'])
nlon = len(ens[0]['lon'])
nt = len(ens[0]['time'])
lat = ens[0]['lat']
lon = ens[0]['lon']
time = ensd[0]['ndcur'] + ensd[0]['nscur'] / 60 / 60 / 24

# ==================== plot ===========

print('plotting')
fig = plt.figure()
gs = gridspec.GridSpec(3,1)
ax = fig.add_subplot(gs[0:2,:]) # plot var vs. time on ax1
ax2 = fig.add_subplot(gs[2,:])  # plot (std / mean) vs. time on ax2
lev = [1000, 100, 50] # these are the levels that will be plotted, in hPa
lat = [0, 35, 70]     # these are the latitudes that will be plotted, in deg
mcolor = 'k'          # color for the ensemble mean

# loop through all ensemble members and compute zonal mean. Again save this result to file, 
# which we check for in the try;except block
print('taking zonal mean of members')
pperti = [0]*N # the zonal mean of each ensemble member will be loaded into this list
for i in range(N):
    print(i)
    try:
        pperti[i] = xr.open_dataset('{}/0{}zmean_{}.nc'.format(out, i+1, var))
        print('read')
    except FileNotFoundError:
        print('computing')
        pperti[i] = ens[i].mean('lon')
        pperti[i].to_netcdf('{}/0{}zmean_{}.nc'.format(out, i+1, var))

# same thing for zonal mean of ensemble mean
print('zonal mean of mean')
try:
    ppmean = xr.open_dataset('{}/ppmean_{}.nc'.format(out, var))[var]
    print('read')
except FileNotFoundError:
    print('computing')
    ppmean = ensmean.mean('lon')
    ppmean.to_netcdf('{}/ppmean_{}.nc'.format(out, var))
print('std of zonal mean')
pstd = np.zeros(ppmean[:,0,0].shape)

# now we will loop through all lat,lev pairs specified for plotting above. For each pair, 
# we will take that lat,lev position in the "3d" zonal-mean time series to reduce the data to 
# a 1d time series per ensemble member (e.g. "zonal mean T at 0deg, 50hPa, for all time")
perti = [0]*N # this list will hold N 1d time series at a time, being re-used for each lat-lev pairing 
for j in range(len(lat)):
    print('---------------- lat {}'.format(lat[j]))
    for k in range(len(lev)):
        print('--- lev {}'.format(lev[k]))
        pmean = ppmean.sel({'lev':lev[k], 'lat':lat[j]}, method='nearest')
        pstd[:] = 0
        for i in range(N):
            perti[i] = pperti[i].sel({'lev':lev[k], 'lat':lat[j]}, method='nearest')[var]
            # plot individual ensemble members with low alpha
            ax.plot(time, perti[i], color='k', alpha=0.1)
            # now compute the std of this 1d time series... we don't actually end up using the 
            # gridpoint-wise std that we computed earlier in this script
            pstd = pstd + ((perti[i] - pmean)**2).values
        pstd = pstd / (N-1)
        pstd = np.sqrt(pstd)
        
        # plot var vs. time on ax1
        ax.plot(time, pmean, color=mcolor, lw=1.5) # plot the ensemble mean
        ax.fill_between(time, pmean-pstd, pmean+pstd, color=color, alpha=0.4) # plot 1std band
        ax.fill_between(time, pmean-2*pstd, pmean+2*pstd, color=color, alpha=0.2) # plot 2std band
        ax.set_ylabel('zonal-mean {}  [{}]'.format(var, unit))
        ax.set_title('lev={} hPa, lat={} degN'.format(lev[k], lat[j]))
        ax.set_xlim([0, 300])
        ax.grid(alpha=0.5) 
        
        #plot (std / mean) vs. time on ax2
        ax2.plot(time, np.abs(pstd/pmean) * 100, color=color, alpha=0.4, label='1std')
        ax2.plot(time, np.abs(2*pstd/pmean) * 100, color=color, alpha=0.2, label='2std')
        ax2.legend(loc='upper left')
        ax2.set_ylabel('{} std/mean  [%]'.format(var))
        ax2.set_xlabel('time  [days]')
        ax2.set_xlim([0, 300])
        ax2.grid(alpha=0.5)
        
        # hardcoded axis limits for different cases, may need to play around and add more
        if(var == 'U' and lev[k] == 1000):
            ax2.set_ylim([-5, 100]) #  for U
        if(var == 'U' and lev[k] < 1000):
            ax2.set_ylim([-5, 50]) #  for U
        if(var == 'T'):
            ax2.set_ylim([-0.5, 2.5]) # for T

        # save figure
        plt.savefig('{}/{}_lat{}_lev{}.png'.format(out, var, lat[j], lev[k]), dpi=300)

        # clear axes for next lat,lev pair
        ax.clear()
        ax2.clear()
