import pdb
import glob
import numpy as np
import xarray as xr
import climate_toolbox as ctb
import matplotlib.pyplot as plt
from scipy.interpolate import griddata as gridd

params = {'text.usetex': False, 'mathtext.fontset': 'stixsans'}
plt.rcParams.update(params)

SCRATCH='/glade/scratch/jhollowed/CAM/cases/sai_runs'

# first let's concatenate year1 outputs (since I had nhtfrq too low)
f1N = 'tmp/dat1N.nc'
f1S = 'tmp/dat1S.nc'
try:
    dat1N = xr.open_dataset(f1N)['U']
    dat1S = xr.open_dataset(f1S)['U']
    print('read year1')
except FileNotFoundError:
    print('computing year1 jet max')
    year1_run = '{}/SE_ne16L72_whs_sai_spinup_FIRSTYEAR/run'.format(SCRATCH)
    year1 = '{}/outputs.concat.zonalMean.nc'.format(year1_run)
    dat1 = ctb.concat_run_outputs(year1_run, histnum=0, mean=['lon'], outFile=year1)
    dat1 = dat1['U'].sel({'lev':slice(0.2, 50)})
    dat1N = dat1.sel({'lat':slice(45, 85)}).max(dim=['lat', 'lev'])
    dat1S = dat1.sel({'lat':slice(-85, -45)}).max(dim=['lat', 'lev'])
    dat1N.to_netcdf(f1N) 
    dat1S.to_netcdf(f1S) 
t1 = ctb.time2day(dat1N['time'])


# get year 2
f2N = 'tmp/dat2N.nc'
f2S = 'tmp/dat2S.nc'
try:
    dat2N = xr.open_dataset(f2N)['U']
    dat2S = xr.open_dataset(f2S)['U']
    print('read year2')
except FileNotFoundError:
    print('computing year2 jet max')
    year2 = '{}/SE_ne16L72_whs_sai_spinup/run/SE_ne16L72_whs_sai_spinup.cam.h0.0001-01-01-00000.nc'\
            .format(SCRATCH)
    dat2 = xr.open_dataset(year2).mean('lon')
    dat2 = dat2['U'].sel({'lev':slice(0.2, 50)})
    dat2N = dat2.sel({'lat':slice(45, 85)}).max(dim=['lat', 'lev'])
    dat2S = dat2.sel({'lat':slice(-85, -45)}).max(dim=['lat', 'lev'])
    dat2N.to_netcdf(f2N) 
    dat2S.to_netcdf(f2S) 
t2 = ctb.time2day(dat2N['time']) + 360

# concat year1, year2
t = np.hstack([t1, t2])
datN = np.hstack([dat1N.values, dat2N.values])
datS = np.hstack([dat1S.values, dat2S.values])

# difference of jet maxima
datdiff = datN - datS

# get datdiff at inithist times
inithist = glob.glob('{}/SE_ne16L72_whs_sai_spinup/remap_inithist/*.i.*'.format(SCRATCH))
tinit_symmetric_mask = np.zeros(len(inithist), dtype=bool)
for i in range(len(inithist)):
    init = xr.open_dataset(inithist[i]).mean('lon').sel({'lev':slice(0.2, 50)})
    initN = init['U'].sel({'lat':slice(45, 85)}).max(dim=['lat', 'lev'])
    initS = init['U'].sel({'lat':slice(-85, -45)}).max(dim=['lat', 'lev'])
    initdiff = abs(initN-initS)
    if initdiff <10:
        tinit_symmetric_mask[i] = True
# get inithist times (first day of each month in year 2)
td2 = np.array([tt.day for tt in dat2N['time'].values])
tm2 = np.array([tt.month for tt in dat2N['time'].values])
td2mask1 = np.logical_and(td2 == 1, tm2 != 1)
td2mask2 = np.logical_and(td2 == 2, tm2 != 1)
tinit = np.hstack([t2[td2mask1], t2[td2mask2]-1])
tinit_symmetric = tinit[tinit_symmetric_mask]


# plot jet symmetry at +-74deg, 3hPa

# -- ax1
fig = plt.figure(figsize=(8.5, 6))
ax = fig.add_subplot(211)
ax.plot(t, datN, color='r', label='northern jet')
ax.plot(t, datS, color='k', label='southern jet')

ylim = ax.get_ylim()
[ax.plot([ti,ti], ylim, '--k', lw=0.85, alpha=0.5) for ti in tinit]
ax.plot([0,0], [0,0], '--k', lw=0.85, alpha=0.5)
[ax.plot([ti,ti], ylim, '-c', lw=1) for ti in tinit_symmetric]
ax.plot([0,0], [0,0], '-c', lw=1)
ax.set_ylim(ylim)
ax.set_xlim([0, 720])

ax.set_ylabel('max($U$) [m/s]\npolar jet zonal wind maxima', fontsize=11)
ax.tick_params(axis="x", direction='in')
ax.xaxis.set_ticklabels([])
ax.legend(frameon=False, fontsize=12, loc='lower left')
ax.set_title('CAM-SE HSW ne16L72 spinup\n(12-hour means in year 1, 48 hours instantaneous in year 2)', 
             fontsize=12)

# -- ax2
ax = fig.add_subplot(212)
ax.plot(t, datN-datS, color='b')
ax.plot(t, np.zeros(len(t)), '-k', label='symmetric jet maxima')

ylim = ax.get_ylim()
[ax.plot([ti,ti], ylim, '--k', lw=0.85, alpha=0.5) for ti in tinit]
ax.plot([0,0], [0,0], '--k', lw=0.85, alpha=0.5, label='inithist output times')
[ax.plot([ti,ti], ylim, '-c', lw=1) for ti in tinit_symmetric]
ax.plot([0,0], [0,0], '-c', lw=1, label='inithist w/ $\Delta U < 10$ m/s')
ax.set_ylim(ylim)
ax.set_xlim([0, 720])

ax.tick_params()
ax.set_ylabel('$(\Delta U)$ [m/s]\ndiff of NH,SH jet maxima', fontsize=11)
ax.set_xlabel('time [days]')
ax2 = ax.secondary_xaxis("top")
ax2.tick_params(axis="x", direction='in')
ax2.xaxis.set_ticklabels([])
ax.legend(frameon=False, fontsize=12, loc='upper left')

plt.tight_layout()
plt.subplots_adjust(hspace=0)
plt.savefig('./figs/check_hsw_jet_symmetry_spinup.png', dpi=300)
