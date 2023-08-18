# Joe Hollowed
# Sandia National Labs 2023
#
# This script analyzes 3-day runs of the HSW++ idealized volcanic injection model
# which integrated the evolution of a diagnostic and tracer potentical vorticity (PV)
# and potential temperature (PT) fields.
#
# This includes steps of preprocessing which does reductions of the data to the required
# spatial means, ensemble mean, and remapping. Preprocessed data is written out to 
# intermediate netcdf files.

import pdb
import glob
import numpy as np
import xarray as xr
import climate_toolbox as ctb
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import ScalarFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable

# -------------------------------------------------------------------------------------

# ---------- locate data ----------
ens_dir = '/pscratch/sd/j/jhollo/E3SM/E3SMv2_cases/pv_cases/highvar_PVPTtest_ens/'
tmp_dir = '{}/post-processed'.format(ens_dir)
# get counterfactual and injection run native and remapped files
ens_native_files  = sorted(glob.glob('{}/injection*/run/*h0*00.nc'.format(ens_dir)))
ens_latlon_files  = sorted(glob.glob('{}/injection*/run/*regrid*'.format(ens_dir)))
cf_native_files   = sorted(glob.glob('{}/counterfactual*/run/*h0*00.nc'.format(ens_dir)))
cf_latlon_files   = sorted(glob.glob('{}/counterfactual*/run/*regrid*'.format(ens_dir)))
# lists for holding dataset objects
N = len(ens_native_files)
ens_native_dat = []
ens_latlon_dat = []
cf_native_dat  = []
cf_latlon_dat  = []

LATSLICE = slice(-10, 30)
LEVSEL = 30

# -----------------------------------------------------------------------------------
# ------------------------- Global inconsistency line plots -------------------------

# ---------- prepare figure ----------
#fig = plt.figure(figsize=(10,12))
#gs  = gridspec.GridSpec(2, 2)
#ax1 = fig.add_subplot(gs[0, 0])
#ax2 = fig.add_subplot(gs[0, 1])
#ax3 = fig.add_subplot(gs[1, :])

fig = plt.figure()
ax1 = fig.add_subplot(111)
fig2 = plt.figure()
ax2 = fig2.add_subplot(111)

# loop through ensemble, take ensemble means, zonal means, and interpolate to isentropes
for i in range(N):

    dat = xr.open_dataset(ens_latlon_files[i]).sel(lat=LATSLICE)
    cfdat = xr.open_dataset(cf_latlon_files[i]).sel(lat=LATSLICE)
    
    if(i == 0):
        lat, lev, lon = dat['lat'], dat['lev'], dat['lon']
        time = ctb.time2day(dat.time)
        hyam, hybm = dat['hyam'], dat['hybm']
        LAT, LEV = np.meshgrid(lat, lev)
        LAT, LON = np.meshgrid(lat, lon)
        
    print('reading...')
    PS = dat['PS']
    PT, PT_TRCR = dat['PT'], dat['PT_TRCR']
    PV, PV_TRCR = dat['PV'], dat['PV_TRCR']
    cfPT, cfPT_TRCR = cfdat['PT'], cfdat['PT_TRCR']
    cfPV, cfPV_TRCR = cfdat['PV'], cfdat['PV_TRCR']
    
    # ----------------------------------
    PVdiff = PV[1::24] - PV_TRCR[1::24]
    PTdiff = PT[1::24] - PT_TRCR[1::24]
    cfPVdiff = cfPV[1::24] - cfPV_TRCR[1::24]
    cfPTdiff = cfPT[1::24] - cfPT_TRCR[1::24]
    PVdiff.to_netcdf('data/PVdiff.nc')
    PTdiff.to_netcdf('data/PTdiff.nc')
    cfPVdiff.to_netcdf('data/cfPVdiff.nc')
    cfPTdiff.to_netcdf('data/cfPTdiff.nc')
    exit()
    # ----------------------------------
    
    PV_lev        = PV.sel(lev=LEVSEL, method='nearest')
    PV_TRCR_lev   = PV_TRCR.sel(lev=LEVSEL, method='nearest')
    cfPV_lev      = cfPV.sel(lev=LEVSEL, method='nearest')
    cfPV_TRCR_lev = cfPV_TRCR.sel(lev=LEVSEL, method='nearest')
    
    PT_lev        = PT.sel(lev=LEVSEL, method='nearest')
    PT_TRCR_lev   = PT_TRCR.sel(lev=LEVSEL, method='nearest')
    cfPT_lev      = cfPT.sel(lev=LEVSEL, method='nearest')
    cfPT_TRCR_lev = cfPT_TRCR.sel(lev=LEVSEL, method='nearest')
    
    print('computing l4 norm...')
    PV_lev_norm   = ctb.lp_norm(PV_lev, PV_TRCR_lev, p=4)
    cfPV_lev_norm = ctb.lp_norm(cfPV_lev, cfPV_TRCR_lev, p=4)
    PT_lev_norm   = ctb.lp_norm(PT_lev, PT_TRCR_lev, p=4)
    cfPT_lev_norm = ctb.lp_norm(cfPT_lev, cfPT_TRCR_lev, p=4)
    #normval = np.max([np.max(PV_lev_norm), np.max(cfPV_lev_norm)])
    #PV_lev_norm   = PV_lev_norm / normval
    #cfPV_lev_norm = cfPV_lev_norm / normval
    #print('{}, {}'.format(np.max(PV_lev_norm).values, np.max(cfPV_lev_norm).values))
    
    print('computing ens mean norm...')
    if(i == 0):
        PV_lev_norm_ensmean   = xr.zeros_like(PV_lev_norm)
        cfPV_lev_norm_ensmean = xr.zeros_like(PV_lev_norm)
        PT_lev_norm_ensmean   = xr.zeros_like(PT_lev_norm)
        cfPT_lev_norm_ensmean = xr.zeros_like(PT_lev_norm)
    PV_lev_norm_ensmean   = PV_lev_norm_ensmean + PV_lev_norm
    cfPV_lev_norm_ensmean = cfPV_lev_norm_ensmean + cfPV_lev_norm
    PT_lev_norm_ensmean   = PT_lev_norm_ensmean + PT_lev_norm
    cfPT_lev_norm_ensmean = cfPT_lev_norm_ensmean + cfPT_lev_norm
    
    print('plotting norm, distribution...')
    ax1.plot(time, PV_lev_norm, '-b', alpha=0.1)
    ax1.plot(time, cfPV_lev_norm, '-k', alpha=0.1)
    ax2.plot(time, PT_lev_norm, '-r', alpha=0.1)
    ax2.plot(time, cfPT_lev_norm, '-k', alpha=0.1)
    
PV_lev_norm_ensmean  = PV_lev_norm_ensmean / N
cfPV_lev_norm_ensmean = cfPV_lev_norm_ensmean / N
ax1.plot(time, PV_lev_norm_ensmean, '-b', lw=2, label='injection ensemble mean')
ax1.plot(time, cfPV_lev_norm_ensmean, '-k', lw=2, label='counterfactual ensemble mean')
ax1.set_xlabel('time [days]')
ax1.set_title('{} hPa'.format(LEVSEL))

PT_lev_norm_ensmean  = PT_lev_norm_ensmean / N
cfPT_lev_norm_ensmean = cfPT_lev_norm_ensmean / N
ax2.plot(time, PT_lev_norm_ensmean, '-r', lw=2, label='injection ensemble mean')
ax2.plot(time, cfPT_lev_norm_ensmean, '-k', lw=2, label='counterfactual ensemble mean')
ax2.set_xlabel('time [days]')
ax2.set_title('{} hPa'.format(LEVSEL))

#ax1.set_ylabel('normalized\ndiagnostic-tracer PV inconsistency')
#fig.savefig('PVinconsistency_normalized.png', dpi=200)

ax1.set_ylabel('diagnostic-tracer PV inconsistency [%]')
fig.savefig('PVinconsistency_{}hPa.png'.format(LEVSEL), dpi=200)

ax2.set_ylabel('diagnostic-tracer PT inconsistency')
fig2.savefig('PTinconsistency_{}hPa.png'.format(LEVSEL), dpi=200)


# -----------------------------------------------------------------------------------
# ------------------------- PV/PT histograms at day 0,15 ----------------------------

fig = plt.figure()
ax = fig.add_subplot()
fig2 = plt.figure()
ax2 = fig2.add_subplot()

default=False
if(LEVSEL == 30):
    nbins = 30
    PV_bins = np.linspace(-2, 6, nbins+1)
    PT_bins = np.linspace(540, 580, nbins+1)
    PV_ylim = [0, 0.55]
    PT_ylim = [0, 0.25]
elif(LEVSEL == 50):
    nbins = 30
    PV_bins = np.linspace(-1.25, 3, nbins+1)
    PT_bins = np.linspace(450, 485, nbins+1)
    PV_ylim = [0, 0.95]
    PT_ylim = [0, 0.33]
elif(LEVSEL == 80):
    nbins = 30
    PV_bins = np.linspace(-0.75, 1.8, nbins+1)
    PT_bins = np.linspace(385, 420, nbins+1)
    PV_ylim = [0, 1.8]
    PT_ylim = [0, 0.52]
elif(LEVSEL == 10):
    nbins = 45
    PV_bins = np.linspace(-7, 23, nbins+1)
    PT_bins = np.linspace(750, 816, nbins+1)
    PV_ylim = [0, 0.15]
    PT_ylim = [0, 0.2]
else:
    PV_bins = 30
    PT_bins = 30
    PV_ylim = [0, 1]
    PT_ylim = [0, 1]
    default=True

PV_hist_args = {'histtype':'step', 'alpha':0.2, 'lw':1, 
                'bins':PV_bins, 'density':True, 'zorder':0}
PV_hist_args_ensmean = {'lw':2, 'bins':PV_bins, 'density':True}
PT_hist_args = {'histtype':'step', 'alpha':0.2, 'lw':1, 
                'bins':PT_bins, 'density':True, 'zorder':0}
PT_hist_args_ensmean = {'lw':2, 'bins':PT_bins, 'density':True}

#final_day = 30
final_day = np.arange(30)

for fday in final_day:
    print('------ day: {}'.format(fday))
    fidx = fday*24
    ax.clear()
    ax2.clear()

    for i in range(N):

        dat = xr.open_dataset(ens_latlon_files[i]).sel(lat=LATSLICE)
        cfdat = xr.open_dataset(cf_latlon_files[i]).sel(lat=LATSLICE)
    
        print('reading...')
        PT, PT_TRCR = dat['PT'], dat['PT_TRCR']
        PV, PV_TRCR = dat['PV'], dat['PV_TRCR']
        cfPT, cfPT_TRCR = cfdat['PT'], cfdat['PT_TRCR']
        cfPV, cfPV_TRCR = cfdat['PV'], cfdat['PV_TRCR']
    
        PV_lev        = PV.sel(lev=LEVSEL, method='nearest')
        PV_TRCR_lev   = PV_TRCR.sel(lev=LEVSEL, method='nearest')
        cfPV_lev      = cfPV.sel(lev=LEVSEL, method='nearest')
        cfPV_TRCR_lev = cfPV_TRCR.sel(lev=LEVSEL, method='nearest')
    
        PT_lev        = PT.sel(lev=LEVSEL, method='nearest')
        PT_TRCR_lev   = PT_TRCR.sel(lev=LEVSEL, method='nearest')
        cfPT_lev      = cfPT.sel(lev=LEVSEL, method='nearest')
        cfPT_TRCR_lev = cfPT_TRCR.sel(lev=LEVSEL, method='nearest')
    
        print('computing ens mean norm...')
        if(i == 0):
            PV_lev_ensmean   = xr.zeros_like(PV_lev)
            cfPV_lev_ensmean = xr.zeros_like(PV_lev)
            PT_lev_ensmean   = xr.zeros_like(PT_lev)
            cfPT_lev_ensmean = xr.zeros_like(PT_lev)
        PV_lev_ensmean   = PV_lev_ensmean + PV_lev
        cfPV_lev_ensmean = cfPV_lev_ensmean + cfPV_lev
        PT_lev_ensmean   = PT_lev_ensmean + PT_lev
        cfPT_lev_ensmean = cfPT_lev_ensmean + cfPT_lev

        ax.hist(np.ravel(PV_lev[fidx].values)/1e-5, color='b', 
                ls=(0, (1, 1)), **PV_hist_args)
        ax.hist(np.ravel(cfPV_lev[0].values)/1e-5, color='k', 
                ls='-', **PV_hist_args)
        ax.hist(np.ravel(cfPV_lev[fidx].values)/1e-5, color='k', 
                ls=(0, (1, 1)), **PV_hist_args)
    
        ax2.hist(np.ravel(PT_lev[fidx].values), color='r', 
                 ls=(0, (1, 1)), **PT_hist_args)
        ax2.hist(np.ravel(cfPT_lev[0].values), color='k', 
                 ls='-', **PT_hist_args)
        ax2.hist(np.ravel(cfPT_lev[fidx].values), color='k', 
                 ls=(0, (1, 1)), **PT_hist_args)

    PV_lev_ensmean  = PV_lev_ensmean / N
    cfPV_lev_ensmean = cfPV_lev_ensmean / N
    
    ax.hist(np.ravel(PV_lev_ensmean[fidx].values)/1e-5, color='b', 
            ls=(0, (1, 1)), **PV_hist_args_ensmean, histtype='step', zorder=99)
    ax.hist(np.ravel(PV_lev_ensmean[fidx].values)/1e-5, color='b', 
            ls=(0, (1, 1)), **PV_hist_args_ensmean, alpha=0.3, zorder=99, 
            label='ensemble mean')
    
    ax.hist(np.ravel(cfPV_lev_ensmean[0].values)/1e-5, color='k', 
            ls='-', **PV_hist_args_ensmean, histtype='step', zorder=50, 
            label='initial distribution')
    
    ax.hist(np.ravel(cfPV_lev_ensmean[fidx].values)/1e-5, color='k', 
            ls=(0, (1, 1)), **PV_hist_args_ensmean, histtype='step', zorder=60)
    ax.hist(np.ravel(cfPV_lev_ensmean[fidx].values)/1e-5, color='k', 
            ls=(0, (1, 1)), **PV_hist_args_ensmean, alpha=0.3, zorder=60, 
            label='counterfactual \nensemble mean')
    
    ax.plot([],[],'-k', lw=PV_hist_args['lw'], alpha=PV_hist_args['alpha'],
            label='ensemble member initial \ndistributions')
    ax.plot([],[],'k', lw=PV_hist_args['lw'], alpha=PV_hist_args['alpha'],
            label='counterfactual ensemble \nmembers', ls=(0, (1, 1)))
    ax.plot([],[],'b', lw=PV_hist_args['lw'], alpha=PV_hist_args['alpha'],
            label='ensemble members', ls=(0, (1, 1)))

    PT_lev_ensmean  = PT_lev_ensmean / N
    cfPT_lev_ensmean = cfPT_lev_ensmean / N
    
    ax2.hist(np.ravel(PT_lev_ensmean[fidx].values), color='r', 
             ls=(0, (1, 1)), **PT_hist_args_ensmean, histtype='step', zorder=99)
    ax2.hist(np.ravel(PT_lev_ensmean[fidx].values), color='r', 
             ls=(0, (1, 1)), **PT_hist_args_ensmean, alpha=0.3, zorder=99, 
             label='ensemble mean')
    
    ax2.hist(np.ravel(cfPT_lev_ensmean[0].values), color='k', ls='-', 
             **PT_hist_args_ensmean, histtype='step', zorder=50, 
             label='initial distribution')
    
    ax2.hist(np.ravel(cfPT_lev_ensmean[fidx].values), color='k', 
             ls=(0, (1, 1)), **PT_hist_args_ensmean, histtype='step', zorder=60)
    ax2.hist(np.ravel(cfPT_lev_ensmean[fidx].values), color='k', 
             ls=(0, (1, 1)), **PT_hist_args_ensmean, alpha=0.3, zorder=60, 
             label='counterfactual \nensemble mean')
    
    ax2.plot([],[],'-k', lw=PT_hist_args['lw'], alpha=PT_hist_args['alpha'],
            label='ensemble member initial \ndistributions')
    ax2.plot([],[],'k', lw=PT_hist_args['lw'], alpha=PT_hist_args['alpha'],
            label='counterfactual ensemble \nmembers', ls=(0, (1, 1)))
    ax2.plot([],[],'r', lw=PT_hist_args['lw'], alpha=PT_hist_args['alpha'],
            label='ensemble members', ls=(0, (1, 1)))

    ax.set_xlabel('PV values [PVU]')
    ax.set_ylabel('pdf')
    ax.set_ylim(PV_ylim)
    if(not default): ax.set_xlim([min(PV_bins), max(PV_bins)])
    ax.set_title('day {}'.format(fday))
    ax.legend(loc='upper right', ncol=2)
    
    ax2.set_xlabel('PT values [K]')
    ax2.set_ylabel('pdf')
    ax2.set_ylim(PT_ylim)
    if(not default): ax2.set_xlim([min(PT_bins), max(PT_bins)])
    ax2.set_title('{} hPa, day {}'.format(LEVSEL, fday))
    ax2.legend(loc='upper right', ncol=2, fontsize=10)
    
    fig.savefig('PV_hists/PV_hist_{}hPa_day{:02d}.png'.format(LEVSEL, fday), dpi=200)
    fig2.savefig('PT_hists/PT_hist_{}hPa_day{:02d}.png'.format(LEVSEL, fday), dpi=200)