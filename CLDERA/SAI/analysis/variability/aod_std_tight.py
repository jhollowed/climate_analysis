'''
Joe Hollowed
University of Michigan, March 2023

Script to plot spread in AOD between members of an input ensemble; AOD distributions are plotted 
in the horizontal and saved as an image for each time sample in the data.
'''

import pdb
import glob
import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from climate_artist import horizontal_slice as plthor
from climate_artist import vertical_slice as pltvert

# ==================== read data ==========
# point to location of ensemble
loc = '/global/cscratch1/sd/jhollo/E3SM/E3SMv2_cases/sai_cases/pertlim_ic_ens/'
# point to location for figure and temporary file outputs
out = '/global/homes/j/jhollo/repos/climate_analysis/CLDERA/SAI/analysis/variability/figs/aod'
N = 5   # number of ensemble members (doesn't need to be hardcoded, could be read from 'loc')

# get AOD data for each ensemble member
# note that we assume the presence in the data directory of remapped lat-lon data files,
# which contain the substring "regrid". We also specify "h1" history files here
ens = [xr.open_dataset(glob.glob('{}/*ens0{}/run/*eam.h1*regrid*'.format(loc, i+1))[0])['AOD'] for i in range(N)]


# compute ensemble mean; sum each ensemble member, divide by N
print('computing ensemble mean')
ensmean = np.zeros(ens[0].shape)
for i in range(N):
    ensmean = ensmean + ens[i]
ensmean = ensmean / N

# compute the ensemble std; 
# sum the squared deviation of each ensemble member, normalize by (N-1), sqrt
print('computing ensemble std')
ensstd = np.zeros(ens[0].shape)
for i in range(N):
    ensstd = ensstd + (ens[i] - ensmean)**2
ensstd = ensstd / (N-1)
ensstd = np.sqrt(ensstd)

# Here we write out teh mean, std data, but these files aren't actually used;
# We could add try;except blocks above to avoid re-computing the mean and std each time this 
# script is run, as in vis_pertlim_spread.py
# If the script is running slow, this would be recommended
print('writing out')
ensmean.to_netcdf('{}/ensmean.nc'.format(out))
ensstd.to_netcdf('{}/ensstd.nc'.format(out))

# get data size and coords
nlat = len(ensmean['lat'])
nlon = len(ensmean['lon'])
nt = len(ensmean['time'])
lat = ensmean['lat']
lon = ensmean['lon']
time = ensmean['time']

# ==================== plot ===========

# time to start plotting at; 
# this is in units of (time samples), and not number of days, so may need to change if different 
# history files are being used, and also if ensemble members with different injection delays are
# used
tstart = 0

# define a "late" flag, which allows us to specify different hard-coded choices of the AOD contour 
# levels for plotting at "late" times
# This was done since I created these plots during the first weeks of the injection, but also for 
# the well-mixed tracer population months after the injection. At these two times, the max AOD values
# are very different
if(tstart > 40): late=True
else: late=False

# --- set contour values for different cases
# these will probably need to be experimented with

# AOD contour levels for each ensemble member
if(not late): clev_c = np.linspace(0.25, 0.75, 3) + 1e-3
else:         clev_c = (np.linspace(0.01, 0.06, 3) + 1e-4) * 5
pltargs_c = {'levels':clev_c, 'colors':None, 'linewidths':1.25}

# AOD contour levels for ensemble mean
if(not late): clev_mean = np.linspace(0, 0.4, 7) + 1e-3
else:         clev_mean = (np.linspace(0, 0.07, 7) + 1e-4 )
cm_mean = plt.cm.Greys
pltargs_mean = {'levels':clev_mean, 'cmap':cm_mean, 'extend':'both', 'alpha':1}

# contour levels for ensemble std
if(not late): clev_std = np.linspace(0, 1, 6) + 1e-3
else:         clev_std = (np.linspace(0, 0.01, 6) + 1e-4) * 10
cm_std = plt.cm.rainbow
pltargs_std = {'levels':clev_std, 'cmap':cm_std, 'extend':'both', 'alpha':1}

LON, LAT = np.meshgrid(lon, lat)
colors = ['r', 'b', 'g', 'm', 'c'] # contour colors for each ensemble member; length must be N


# --- start plotting; outer loop = time, inner loop = ensemble members
for j in range(nt-tstart):
    
    print('plotting time {}...'.format(j+1))
    print('creating figure...')
    fig = plt.figure(figsize=(12, 3))
    ax1 = fig2.add_subplot(121, projection=ccrs.Mercator()) # AOD contour plots
    ax2 = fig2.add_subplot(122, projection=ccrs.Mercator()) # std plot
    fig.suptitle('Day {:.0f}'.format((tstart+j) * 2), fontsize=14) 
       
    # the plotting is done with my own custom plotting routines which wrap 
    # pyplot and cartopy functions; arguments are passed via dictionaries
    # these functions are decently-well documented in inline comments, see
    # the 'horizontal_slice' class of 'climate_artist'
    # alternatively, these calls could be removed and replaced with standard
    # matplotlib utilities...

    # plot ens mean AOD contours
    var_dict = [{'var':ensmean[j+tstart], 'plotType':'contourf', 
                 'plotArgs':pltargs_mean, 'colorFormatter':None}]
    plthor(lon, lat, var_dict, ax=ax1, annotation=None, include_contours=False, 
           coastlinesArgs={'lw':0.5, 'alpha':0.75})

    # plot std colormap 
    var_dict = [{'var':ensstd[j+tstart], 'plotType':'contourf', 
        'plotArgs':pltargs_std}]
    plthor(lon, lat, var_dict, ax=ax2, annotation=None, include_contours=True, 
           coastlinesArgs={'lw':0.5, 'alpha':0.75})

    # loop through ensemble members
    for i in range(N):
        # loop through desired contour levels; each one will be plotted individually
        # so it can be assigned a decreasing alpha value
        for k in range(len(clev_c)):
            pltargs_c['colors']=colors[i]
            pltargs_c['alpha']=(k+1)/len(clev_c)
            pltargs_c['levels']=[clev_c[k]]
            var_dict = [{'var':ens[i][j+tstart], 'plotType':'contour', 
                         'plotArgs':pltargs_c, 'colorFormatter':None}]
            plthor(lon, lat, var_dict, ax=ax1, annotation=None, include_contours=False, 
                   coastlines=False)
        # dummy plot for ensemble member legend label
        ax1.plot([0,0], [0,0], color=colors[i], label='ens0{}'.format(i+1))
    
    # set the cartopy plot extent to focus on the tropical to sub-tropical regions, with 
    # preference toward the NH. Values in tuple arg are latitudes
    ax1.set_extent((-179, 180, -40, 60))
    ax2.set_extent((-179, 180, -40, 60))
    # this is a bit of an abuse, but cartopy by default enforces a 1:1 scale for longitude 
    # and latitude, which makes the figures awkwardly long and narrow; adjust the aspect
    # the be something more reasonable here (though continent coastlines will now be distorted)
    ax1.set_aspect(1.5)
    ax2.set_aspect(1.65)
    
    # all done
    ax1.legend(loc='lower left', ncol=2)
    fig.tight_layout()
    plt.show()
    #fig.savefig('{}/aod_std{:03d}.png'.format(out, j))

