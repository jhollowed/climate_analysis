import os
import sys

print(sys.executable)

import pdb
import copy
import glob
import warnings
import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import cartopy.util as cutil
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import colors
import matplotlib.ticker as mticker
from matplotlib.gridspec import GridSpec
import matplotlib.transforms as mtransforms
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

sys.path.append(os.path.abspath("/global/homes/j/jhollo/repos/climate_util"))
from PyTEMDiags import sph_zonal_averager
import climate_toolbox as ctb
import artist_utils as au

# ---------- matplotlib settings

mpl.rcParams.update(mpl.rcParamsDefault)

SMALL_SIZE  = 9
MEDIUM_SIZE = 11
BIG_SIZE    = 12
BIGGER_SIZE = 14

plt.rc('font', size=MEDIUM_SIZE)         # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)    # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

params = {'text.usetex': False, 'mathtext.fontset': 'stixsans'}
plt.rcParams.update(params)

print('done')

# ==================================
# ========= read data ==============

data_dest = '/pscratch/sd/j/jhollo/tmp'

print('reading ensemble...')
ens_members = sorted(glob.glob('/global/cfs/cdirs/m4014/data/HSW/outputs/release_011423/netcdf/ens_members/ens*/*h2*.nc'))
ens_members = [xr.open_dataset(mem) for mem in ens_members]
ens_members = [mem.assign_coords(time = ctb.time2day(mem['time'])) for mem in ens_members]

time   = ens_members[0]['time']
datlat = ens_members[0]['lat']
lat    = np.linspace(-90, 90, 181)

AOD_mems    = [mem['AOD'] for mem in ens_members]

print('taking AOD zonal means...')
zm = sph_zonal_averager(datlat, lat, L=70, debug=False)
zm.sph_compute_matrices()
AOD_mems = [zm.sph_zonal_mean(AOD_mem.T).T for AOD_mem in AOD_mems]

print('getting ensemble mean AOD...')
AOD_em = xr.zeros_like(AOD_mems[0])
for i in range(len(AOD_mems)):
    AOD_em = AOD_em + AOD_mems[i]
AOD_em = AOD_em / len(AOD_mems)

print('done')

# ===============================================
# ========= read supplemental data ==============

print('reading ensemble rerun...')
ens_members = sorted(glob.glob('/pscratch/sd/j/jhollo/E3SM/E3SMv2_cases/sai_cases/ic_ens/extra_outputs_for_hsw_paper/ens_members/*541day*ens*/run/*.h2.*.nc'))
ens_members = [xr.open_dataset(mem) for mem in ens_members]
ens_members = [mem.assign_coords(time = ctb.time2day(mem['time'])) for mem in ens_members]

time2 = ens_members[0]['time']
TIME2, LAT2 = np.meshgrid(time2, lat)

# get heating at lowest model level
heat_mems   = [mem['SAI_HEAT'].sel(lev=1000, method='nearest') for mem in ens_members]

print('getting zonal mean HEAT...')
heat_mems    = [zm.sph_zonal_mean(mem.T).T for mem in heat_mems]

print('getting ensemble mean HEAT...')
heat_em = xr.zeros_like(heat_mems[0])
for i in range(len(heat_mems)):
    heat_em = heat_em + heat_mems[i]
heat_em = heat_em / len(heat_mems)

print('done')

print('cleaning weird polar artifacts at injection time...')
for i in range(len(lat)):
    if np.abs(lat[i]-15.15) > 15:
        for k in range(len(time)):
            if(time[k] < 183): AOD_em[k,i] = 0
    if np.abs(lat[i]-15.15) > 25:
        for k in range(len(time2)):
            if time2[k] < 183: heat_em[k,i] = 0
print('done')

# =============================================
# ===============  plotting  ==================

TIME, LAT   = np.meshgrid(time, lat)

fig   = plt.figure(figsize=(6, 5))
gs    = fig.add_gridspec(nrows=5, ncols=1, hspace=0.4)
ax1   = fig.add_subplot(gs[:2])
ax2   = fig.add_subplot(gs[2:])

erupt_time = 180
erupt_pos  = 15.15
#AODcmap    = au.WhiteYellowOrangeRed
AODcmap    = plt.cm.OrRd
AODcmap    = plt.cm.GnBu

# --- for matching fraction of pre-eruption space taking in plots
minlin, maxlin = 170, 270
maxlog         = 1000
erupt_time     = 180
f              = (180-minlin)/maxlin
minlog         = (maxlog*f)/(f-1) - erupt_time/(f-1)

# -----------------------------------------------
# --------- UPPER LINEAR-SCALED PANEL -----------

timelim   = [minlin, maxlin]
timeticks = -np.arange(-maxlin, -minlin+15, 15)[::-1]
latlim    = [-30, 90]
latticks  = [-30, 0, 30, 60, 90]

AOD_var   = AOD_em.T
AODlevels = np.arange(0.05, 0.351, 0.05)
AODlw     = 1
AODccolor = 'w'
AODcalpha = 0.66

cbarlabels  = list(AODlevels)
for i in range(len(cbarlabels)):
    l = AODlevels[i]
    if(np.isclose((l*10)%1, 0)): 
        cbarlabels[i] = str(np.around(l, decimals=1))
    else: 
        cbarlabels[i] = ' '

H_var     = heat_em.T
Hcolor    = 'darkblue'
Hcolor='crimson'
Hlw       = 1
Halpha    = 1
#Hlevels   = -np.arange(0.1, 0.41, 0.1)[::-1]
#Hlabpos   = [[230, 20], [252, 30], [255, 55], [260, 70]]
#Hlabrot   = [0, 0, 0, 0]
#Hlabspace = [-10.5, -10.5, -10, -10]
Hlevels   = [-0.1, -0.25, -0.4][::-1]
Hlabpos   = [[230, 20], [247, 40], [260, 70]]
Hlabrot   = [0, 0, 0]
Hlabspace = [-10.5, -13.5, -10]
Hlabel    = r'cooling rate [K day$^{-1}$]'
Hls       = '-.'

cf  = ax1.contourf(TIME, LAT, AOD_var, cmap=AODcmap, levels=AODlevels, vmax=np.max(AODlevels)*1.25, extend='max')
cfl = ax1.contour(TIME, LAT, AOD_var, levels=AODlevels, colors=AODccolor, alpha=AODcalpha, linewidths=AODlw, linestyles='-')
cfl2 = ax1.contour(TIME2, LAT2, H_var, levels=Hlevels, colors=Hcolor, linewidths=Hlw, linestyles=Hls, alpha=Halpha)

ax1.plot([0], [-1000], ls=Hls, color=Hcolor, lw=Hlw, alpha=Halpha, label=Hlabel) # dummy

ax1.plot(erupt_time, erupt_pos, '^k', ms=13, mew=1.5, mec='w')

cbar = fig.colorbar(cf, ax=ax1, orientation='vertical', location='right', extendfrac=0.025,
                    extend='both', extendrect=False, aspect=8)
cbar.set_label('AOD')
cbar.ax.set_yticks(AODlevels, labels=cbarlabels)

#ax1.axhline(0, linestyle=':', color='k', linewidth=1.5, alpha=0.5)

for i in range(len(Hlabspace)):
    labs = ax1.clabel(cfl2, inline=True, fontsize=SMALL_SIZE, manual=[Hlabpos[i]], 
                      inline_spacing=Hlabspace[i], fmt=lambda x: '{:.3g}'.format(x))
for i in range(len(labs)):
    labs[i].set_rotation(Hlabrot[i])
     
ax1.set_yticks(latticks)
ax1.set_ylim(latlim)
ax1.set_ylabel('latitude')
ax1.set_xticks(timeticks)
ax1.set_xlim(timelim)
ax1.set_xlabel('time [days]')
ax1.tick_params(top=True, bottom=True, left=True, right=True, 
                labelleft=True, labelright=False, labeltop=True, labelbottom=False, which='both')
ax1.xaxis.set_label_position('top') 

ax1.legend(loc='upper left', ncol=2, fancybox=False, shadow=False, frameon=False)

# --------------------------------------------
# --------- LOWER LOG-SCALED PANEL -----------

AOD_var       = np.log10(AOD_em.T)
AODlevels     = np.arange(-2.5, -0.45, 0.25)
AODboldlevels = [-2, -1]
AODclevels    = np.setdiff1d(AODlevels, AODboldlevels)
    
cbarlabels  = list(AODlevels)
for i in range(len(cbarlabels)):
    l = AODlevels[i]
    if(np.isclose((l*10)%1, 0)): 
        cbarlabels[i] = str(np.around(l, decimals=1))
    else: 
        cbarlabels[i] = ' '

AODlw      = 0.75
AODccolor  = 'k'
AODcalpha  = 0.5

timelim   = [minlog, maxlog]
timeticks = -np.arange(-max(timelim), -minlog+100, 100)[::-1]
latlim    = [-90, 90]
latticks  = [-90, -60, -30, 0, 30, 60, 90]

cf  = ax2.contourf(TIME, LAT, AOD_var, levels=AODlevels, cmap=AODcmap, vmax=np.max(AODlevels)*0.5)
cfl = ax2.contour(TIME, LAT, AOD_var, levels=AODclevels, colors=AODccolor, alpha=AODcalpha, linewidths=AODlw, linestyles='-')
cfl_bold = ax2.contour(TIME, LAT, AOD_var, levels=AODboldlevels, colors=AODccolor, alpha=0.75, linewidths=AODlw*2, linestyles='-')

ax2.plot(erupt_time, erupt_pos, '^k', ms=13, mew=1.5, mec='w')

cbar = fig.colorbar(cf, ax=ax2, orientation='vertical', location='right', extendfrac=0.025,
                    aspect=12, extend='both', extendrect=False)
cbar.set_label(r'$\mathit{log}_{10}$(AOD)')
cbar.ax.set_yticks(AODlevels, labels=cbarlabels)

labs = ax2.clabel(cfl_bold, inline=True, fontsize=SMALL_SIZE, inline_spacing=4, 
                  manual=[[650, -30], [525, 35]])

ax2.axhline(0, linestyle=':', color='k', linewidth=1.5, alpha=0.5)
        
ax2.set_ylim(latlim)
ax2.set_ylabel('latitude')
ax2.set_yticks(latticks)
ax2.set_xlim(timelim)
ax2.set_xlabel('time [days]')
ax2.set_xticks(timeticks)
ax2.tick_params(top=True, bottom=True, left=True, right=True, 
                labelleft=True, labelright=False, labeltop=False, labelbottom=True, which='both')

# ------------------------------------------

axs = [ax1, ax2]
axlabs = ['(a)', '(b)']
for i in range(len(axs)):
    trans = mtransforms.ScaledTranslation(-20/72, 7/72, fig.dpi_scale_trans)
    pos = [-0.07, 0.87]
    axs[i].text(pos[0], pos[1], axlabs[i], transform=axs[i].transAxes + trans,
               fontsize=BIGGER_SIZE, va='bottom', fontfamily='sans-serif')
    
#ax1.grid()

plt.tight_layout()
fig.savefig('figs/AODnew.pdf', dpi=300, bbox_inches='tight')
plt.show()