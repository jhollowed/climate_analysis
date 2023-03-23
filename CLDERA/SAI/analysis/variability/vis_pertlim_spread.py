import numpy as np
import glob
import pdb
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import xarray as xr
from metpy.units import units as u
import metpy.constants as const

loc = '/project/projectdirs/m4014/data/HSW/outputs/release_030123/pert_ens'
out = 'figs'
var = 'U'
if(var == 'T'):
    unit = 'K'
    color = 'r'
if(var == 'U'):
    unit = 'm/s'
    color = 'b'
N = 8

ensd = [xr.open_dataset(glob.glob('{}/pert0{}/*h1*regrid*'.format(loc, i+1))[0]) for i in range(N)]
ens = [ensd[i][var] for i in range(N)]

try:
    ensmean = xr.open_dataset('{}/pertmean_{}.nc'.format(out, var))[var]
    ensstd = xr.open_dataset('{}/pertstd_{}.nc'.format(out, var))[var]
    print('read')
except FileNotFoundError:
    ensmean = ens[0].copy(deep=True, data=np.zeros(ens[0].shape))
    ensstd = ens[0].copy(deep=True, data=np.zeros(ens[0].shape))
    print('computing ensemble mean')
    for i in range(N):
        ensmean = ensmean + ens[i]
    ensmean = ensmean / N

    print('computing ensemble std')
    for i in range(N):
        ensstd = ensstd + (ens[i] - ensmean)**2
    ensstd = ensstd / (N-1)
    ensstd = np.sqrt(ensstd)

    print('writing out')
    ensmean.to_netcdf('{}/pertmean_{}.nc'.format(out, var))
    ensstd.to_netcdf('{}/pertstd_{}.nc'.format(out, var))
    
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
ax = fig.add_subplot(gs[0:2,:])
ax2 = fig.add_subplot(gs[2,:])
lev = [1000, 100, 50]
lat = [0, 35, 70]
mcolor = 'k'
ls = ['-', '--', ':']

print('zonal mean of mems')
pperti = [0]*N
for i in range(N):
    print(i)
    try:
        pperti[i] = xr.open_dataset('{}/0{}zmean_{}.nc'.format(out, i+1, var))
        print('read')
    except FileNotFoundError:
        print('computing')
        pperti[i] = ens[i].mean('lon')
        pperti[i].to_netcdf('{}/0{}zmean_{}.nc'.format(out, i+1, var))

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

perti = [0]*N
for j in range(len(lat)):
    print('---------------- lat {}'.format(lat[j]))
    for k in range(len(lev)):
        print('--- lev {}'.format(lev[k]))
        pmean = ppmean.sel({'lev':lev[k], 'lat':lat[j]}, method='nearest')
        pstd[:] = 0
        for i in range(N):
            perti[i] = pperti[i].sel({'lev':lev[k], 'lat':lat[j]}, method='nearest')[var]
            ax.plot(time, perti[i], color='k', alpha=0.1)
            pstd = pstd + ((perti[i] - pmean)**2).values
        pstd = pstd / (N-1)
        pstd = np.sqrt(pstd)

        ax.plot(time, pmean, color=mcolor, lw=1.5)
        ax.fill_between(time, pmean-pstd, pmean+pstd, color=color, alpha=0.4)
        ax.fill_between(time, pmean-2*pstd, pmean+2*pstd, color=color, alpha=0.2)
        ax.grid(alpha=0.5)
        ax.set_ylabel('zonal-mean {}  [{}]'.format(var, unit))
        ax.set_title('lev={} hPa, lat={} degN'.format(lev[k], lat[j]))

        ax2.plot(time, np.abs(pstd/pmean) * 100, color=color, alpha=0.4, label='1std')
        ax2.plot(time, np.abs(2*pstd/pmean) * 100, color=color, alpha=0.2, label='2std')
        ax2.grid(alpha=0.5)
        if(var == 'U' and lev[k] == 1000):
            ax2.set_ylim([-5, 100]) #  for U
        if(var == 'U' and lev[k] < 1000):
            ax2.set_ylim([-5, 50]) #  for U
        if(var == 'T'):
            ax2.set_ylim([-0.5, 2.5]) # for T
        ax2.legend(loc='upper left')
        ax2.set_ylabel('{} std/mean  [%]'.format(var))
        ax2.set_xlabel('time  [days]')

        ax.set_xlim([0, 300])
        ax2.set_xlim([0, 300])
        plt.savefig('{}/{}_t300_lat{}_lev{}.png'.format(out, var, lat[j], lev[k]), dpi=300)

        ax.set_xlim([0, 120])
        ax2.set_xlim([0, 120])
        plt.savefig('{}/{}_t120_lat{}_lev{}.png'.format(out, var, lat[j], lev[k]), dpi=300)

        ax.clear()
        ax2.clear()
