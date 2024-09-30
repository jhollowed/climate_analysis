import os
import pdb
import sys
import glob
import math
import numpy as np
import xarray as xr
import climate_toolbox as ctb
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import ScalarFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
import wrappers
import climate_toolbox as ctb

# -------------------------------------------------------------------------------------

# ---------- locate data ----------
fig_dir = '/global/homes/j/jhollo/repos/climate_analysis/CLDERA/PV/analysis/historical_figs'
ens_dir = '/pscratch/sd/j/jhollo/E3SM/E3SMv2_cases/historical_cases/'
tmp_dir = '{}/post-processed'.format(ens_dir)

if(len(sys.argv) == 1): raise RuntimeError('ensemble member number or mean not specified')
ensNum = sys.argv[1] #either an int, or the string 'mean'
if(ensNum == 'mean'): ensMean = True

timeSamples = 30
LEVSEL=50
clev = 10

do_isentropes = False
lev_spec = ['_plev', '_thetalev']
sfx = lev_spec[do_isentropes]

refdat = xr.open_dataset(glob.glob('{}/*ens1/run/*eam.h0*'.format(ens_dir))[0])
dat    = xr.open_dataset('{}/dat_native_ens{}{}.nc'.format(tmp_dir, ensNum, sfx))
cfdat  = xr.open_dataset('{}/cfdat_native_ens{}{}.nc'.format(tmp_dir, ensNum, sfx))
cfdiff = xr.open_dataset('{}/dat_native_ens{}_cfdiff{}.nc'.format(tmp_dir, ensNum, sfx))

time = dat['time']
day = ctb.time2day(time)
lat  = refdat['lat']
lon  = refdat['lon']
timeIdx = np.round(np.linspace(0, len(time) - 1, timeSamples)).astype(int)

for j in timeIdx:
    
    fig = plt.figure(figsize=(10,5))
    axpt1 = fig.add_subplot(231)
    axpt2 = fig.add_subplot(232)
    axpt3 = fig.add_subplot(233)
    axpv1 = fig.add_subplot(234)
    axpv2 = fig.add_subplot(235)
    axpv3 = fig.add_subplot(236)

    axpt = [axpt1, axpt2, axpt3]
    ptvars = [dat['PT'], dat['PT_INCONSISTENCY'], cfdiff['PT_INCONSISTENCY']]
    ptcmap = [plt.cm.plasma, plt.cm.YlOrRd, plt.cm.YlOrRd]
    pttitles = ['PT', 'PT inconsistency', 'PT inconsistency diff']
    ptclev = [np.linspace(450, 550, clev), np.linspace(-50, 60, clev), np.linspace(-25, 25, clev)]

    axpv = [axpv1, axpv2, axpv3]
    pvvars = [dat['PV']/1e-5, dat['PV_INCONSISTENCY']/1e-5, cfdiff['PV_INCONSISTENCY']/1e-5]
    pvcmap = [plt.cm.viridis, plt.cm.YlGnBu, plt.cm.YlGnBu]
    pvtitles = ['PV', 'PV inconsistency', 'PV inconsistency diff']
    pvclev = [np.linspace(-7, 7, clev), np.linspace(-30, 30, clev), np.linspace(-20, 20, clev)]


    print('--- working on time {}, lev {}'.format(j, LEVSEL))
    fig.suptitle("lev = {} hPa, day {:.2f}".format(LEVSEL, day[j]), fontsize=14)
    selargs = {'indexers':{'time':time[j], 'lev':LEVSEL}, 'method':'nearest'}
    
    for ax in axpv: ax.set_xlabel('lon')
    for i in range(len(axpt)):
        axpt[i].set_title(pttitles[i])
        axpv[i].set_title(pvtitles[i])
        if i == 1: 
            axpt[i].yaxis.set_ticklabels([])
            axpv[i].yaxis.set_ticklabels([])
        if i == 2: 
            axpt[i].yaxis.tick_right()
            axpv[i].yaxis.tick_right()
        if i != 1: 
            axpt[i].set_ylabel('lat')
            axpv[i].set_ylabel('lat')

    for i in range(3):
        print('var {}...'.format(i))
        cpt = axpt[i].tricontourf(lon, lat, ptvars[i].sel(**selargs), cmap=ptcmap[i], levels=ptclev[i], extend='both')
        #axpt[i].tricontour(lon, lat, ptvars[i].sel(**selargs), colors='k', levels=clev, linewidths=0.5)
        cpv = axpv[i].tricontourf(lon, lat, pvvars[i].sel(**selargs), cmap=pvcmap[i], levels=pvclev[i], extend='both')
        #axpv[i].tricontour(lon, lat, pvvars[i].sel(**selargs), colors='k', levels=clev, linewidths=0.5)

        cbpt = fig.colorbar(cpt, ax=axpt[i], location='bottom')
        cbpv = fig.colorbar(cpv, ax=axpv[i], location='bottom')

    plt.tight_layout()
    plt.savefig('{}/day{}.png'.format(fig_dir, day[j]))
    #if(j == timeIdx[-1]): plt.show()
    
