import sys
import pdb
import glob
import scipy
import numpy as np
import xarray as xr
from matplotlib import colors
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import ScalarFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
import wrappers

TICK_FS = 11
LABEL_FS = 14
plt.rc('font', size=TICK_FS)          # controls default text sizes
plt.rc('axes', titlesize=LABEL_FS)    # fontsize of the axes title
plt.rc('axes', labelsize=LABEL_FS)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=TICK_FS)    # fontsize of the tick labels
plt.rc('ytick', labelsize=TICK_FS)    # fontsize of the tick labels
plt.rc('legend', fontsize=TICK_FS)    # legend fontsize
plt.rc('figure', titlesize=LABEL_FS)  # fontsize of the figure title

# --------------------------------------------------------------


# --------------- ANALYSIS SETTINGS ---------------
WEIGHT_LAT = True   # whether or not to weight horizontal averages in by cos(lat)
WEIGHT_LEV = True   # whether or not to weight veritcal averages by pdel
DO_ANNUAL  = False  # whether or not to use annually-averaged data fror this analysis
                    # if true, climatological anomalies will not be taken in E90, ST80, 
                    # or the tropopause; averaged raw ppb values will be plotting instead

        
# --------------- read data ---------------

print('---- reading data...')
pdir         = '/pscratch/sd/j/jhollo/E3SM/historical/CLDERA_2ndHistorical_1950-2014_monthly'
# full data (we'll pull meta data, coordiante info, etc from here)
E90     = '{}/histoircal_h0_E90j_1850-2014.regrid.fv0.9x1.25_bilinear.nc'.format(pdir)
# zonal mean data
zm_E90  = '{}/histoircal_h0_E90j_1850-2014.regrid.fv0.9x1.25_bilinear.zonalmean.nc'.format(pdir)
zm_ST80 = '{}/histoircal_h0_ST80_25j_1850-2014.regrid.fv0.9x1.25_bilinear.zonalmean.nc'.format(pdir)
zm_T    = '{}/histoircal_h0_T_1850-2014.regrid.fv0.9x1.25_bilinear.zonalmean.nc'.format(pdir)
zm_PS   = '{}/histoircal_h0_PS_1850-2014.regrid.fv0.9x1.25_bilinear.zonalmean.nc'.format(pdir)

# get coordinates etc.
dat            = xr.open_dataset(E90) 
lat, lev       = dat['lat'], dat['lev'] 
P0, hyai, hybi = dat['P0'], dat['hyai'], dat['hybi']

# get zonally averaged data
e90   = xr.open_dataset(zm_E90)['E90j']      # in ppb
st80  = xr.open_dataset(zm_ST80)['ST80_25j'] # in ppb
T     = xr.open_dataset(zm_T)['T']
PS    = xr.open_dataset(zm_PS)['PS']


# --------------- get monthly climatology ---------------
# these are climatological averages of each zonally-averaged variable
# result is a 3D dataset in (lat, lev, month) with 12 months

e90_climo  = e90.groupby('time.month').mean('time')
st80_climo = st80.groupby('time.month').mean('time')
T_climo    = T.groupby('time.month').mean('time')
PS_climo   = PS.groupby('time.month').mean('time')


# --------------- get annual means ---------------
# if DO_ANNUAL is true, here we take annual means of the data
# if not, we do nothing other than to add a "year" coordinate to
# the monthly-mean data for later use, which is decimal-valued and
# includes the current month. This was done to make the use of "year"
# in the following analysis and plotting steps more consistent

if(DO_ANNUAL):
    print('---- taking annnual means...')
    e90  = e90.groupby('time.year').mean('time')
    e90  = e90.rename({'year':'time'})
    e90  = e90.assign_coords(year = e90['time'])
    st80 = st80.groupby('time.year').mean('time')
    st80 = st80.rename({'year':'time'})
    st80 = st80.assign_coords(year = st80['time'])
    T    = T.groupby('time.year').mean('time')
    T    = T.rename({'year':'time'})
    T    = T.assign_coords(year = T['time'])
    PS   = PS.groupby('time.year').mean('time')
    PS   = PS.rename({'year':'time'})
    PS   = PS.assign_coords(year = PS['time'])
else:
    e90  = e90.assign_coords(year = (e90['time.year'] + e90['time.month']/12))
    st80 = st80.assign_coords(year = (st80['time.year'] + st80['time.month']/12))
    T    = T.assign_coords(year = (T['time.year'] + T['time.month']/12))
    PS   = PS.assign_coords(year = (PS['time.year'] + PS['time.month']/12))

    
# --------------- find tropopause as function of time, lat ---------------
# we will need to the tropopause position to analyze, as well as for masking
# the tracer quantities into "strat" and "trop" populatiions.
# the tropopause is found with NCL's "trop_wmo" builtin, which computes the 
# pressure at the thermal (static) tropopause for pressures and temperatures 
# following the 1992 WMO definition of the height of the tropopause:
# 
# The first tropopause is defined as the lowest level at which
# the lapse rate decreases to 2 deg K per kilometer or less,
# provided also the average lapse rate between this level and
# all higher levels within 2 kilometers does not exceed 2 deg K.

print('---- finding tropopause...') 
trop_t = wrappers.trop_wmo(T.lev.values, 
                           T.transpose('time', 'lat', 'lev')).values

# now create troposphere/stratosphere masks. These will evolve in time, 
# so that on each month, these data is masks according to the current 
# monthly-mean tropopause

# this is basically a "meshgrid" of lev...
LEV = (xr.broadcast(T.time, T.lat, T.lev)[2]) 
# transposes and expansion of the of the array by a new dimension of length 72 
# needed to make the return of trop_wmo agree in shape/size with LEV 
# (for later masking)...
trop_mask  = LEV >  np.array([trop_t.T] * 72).T
trop_mask  = trop_mask.transpose('time', 'lev', 'lat')
# stratosphere mask is the inverse of the troposphere mask...
strat_mask = ~trop_mask

# now do the same thing for the climatological tropopause (will allow us 
# to later find tropopuse height anomalies)
LEV_climo = (xr.broadcast(T_climo.month, T_climo.lat, T_climo.lev)[2])
trop_t_climo = wrappers.trop_wmo(T_climo.lev.values, 
                                 T_climo.transpose('month', 'lat', 'lev')).values
trop_mask_climo  = LEV_climo >  np.array([trop_t_climo.T] * 72).T
trop_mask_climo  = trop_mask_climo.transpose('month', 'lev', 'lat')
strat_mask_climo = ~trop_mask_climo

# finally, let's find the mean pressure position of the tropically-averaged 
# tropopuase (on +- 20 degrees). This will be a 2D time series in (lev-time)
eq_mask = np.logical_and(lat > -20, lat < 20)
trop_t_eq_mean = np.mean(trop_t.T[eq_mask], axis=0)
trop_t_eq_mean_time = T.year.values
trop_t_climo_eq_mean = np.mean(trop_t_climo.T[eq_mask], axis=0)


# --------------- get e90, st80 in troposphere, stratosphere ---------------
# do this by simply applying the masks computed above. Array positions removed 
# by the mask will be filled wiht zeros

print('---- masking tracer fields...')
e90_trop  = e90.where(trop_mask, 0)
e90_strat = e90.where(strat_mask, 0)
st80_trop = st80.where(trop_mask, 0)
st80_strat = st80.where(strat_mask, 0)
e90_climo_trop   = e90_climo.where(trop_mask_climo, 0)
e90_climo_strat  = e90_climo.where(strat_mask_climo, 0)
st80_climo_trop  = st80_climo.where(trop_mask_climo, 0)
st80_climo_strat = st80_climo.where(strat_mask_climo, 0)

if(0):
    # debugging plots which will show the binary masks, the 
    # temperature, tracer fields, and masked tracer fields
    debug_level = 10
    plt.imshow(trop_mask[debug_level]); plt.savefig('trop10.png')
    plt.imshow(T[debug_level]); plt.savefig('t10.png')
    plt.imshow(e90_trop[debug_level]); plt.savefig('et10.png')
    plt.imshow(e90_strat[debug_level]); plt.savefig('es10.png')
    plt.imshow(e90[debug_level]); plt.savefig('e10.png')
    plt.imshow(st80_trop[debug_level]); plt.savefig('stt10.png')
    plt.imshow(st80_strat[debug_level]); plt.savefig('sts10.png')
    plt.imshow(st80[debug_level]); plt.savefig('st10.png')
    
    
# --------------- take weighted tropical averages ---------------
# here we compute tropical-mean averages for the tracers and their 
# climatologies. This will result in 2D time series in (lev-height).
# Contouring these later will give e.g. mean tropical E90 90ppb position

print('---------- computing tropical means...')
# ---- lat weighting
e90_eq  = e90.sel({'lat':slice(-20, 20)})
st80_eq = st80.sel({'lat':slice(-20, 20)})
e90_climo_eq  = e90_climo.sel({'lat':slice(-20, 20)})
st80_climo_eq = st80_climo.sel({'lat':slice(-20, 20)})
if(WEIGHT_LAT):
    # weight the datasets in laitutde by cos(lat)...
    print('---- weighting lat...')
    weights = np.cos(np.deg2rad(lat))
    weights.name = "weights"
    e90_eq_weighted  = e90_eq.sel({'lat':slice(-20, 20)}).weighted(weights)
    st80_eq_weighted = st80_eq.sel({'lat':slice(-20, 20)}).weighted(weights)
    e90_climo_eq_weighted  = e90_climo_eq.sel({'lat':slice(-20, 20)}).weighted(weights)
    st80_climo_eq_weighted = st80_climo_eq.sel({'lat':slice(-20, 20)}).weighted(weights)
else:
    e90_eq_weighted  = e90_eq
    st80_eq_weighted = st80_eq
    e90_climo_eq_weighted  = e90_climo_eq
    st80_climo_eq_weighted = st80_climo_eq
    
# done computing horizontal weights, do horizontal averaging...
e90_eq_mean   = e90_eq_weighted.mean('lat')
st80_eq_mean  = st80_eq_weighted.mean('lat')
e90_climo_eq_mean   = e90_climo_eq_weighted.mean('lat')
st80_climo_eq_mean  = st80_climo_eq_weighted.mean('lat')

    
# --------------- take weighted global averages ---------------
# here we compute global averages of the strat/trop masked tracer fields
# These global averages are optionally weighted in both lat and lev, controlled
# by the flags WEIGHT_LEV and WEIGHT_LAT. These should always be on; the toggle
# was only implemented for sanity checking

print('---------- computing global means...')
# ---- lev weighting
if(WEIGHT_LEV):
    # weight the datasets in laitutde by pressure thickness...
    # the pressure thickness is computed properly for each gridpoint
    # usign NCL's dpres_hybrid_ccm. This function computes the actual
    # pressure at each gridpoint using PS, P0, and the coordinate 
    # hybrid coefficients
    print('---- weighting lev...')
    weights = wrappers.dpres_hybrid_ccm(PS, P0, 
                                        hyai, hybi).transpose('ncl2', 'ncl1', 'ncl3')
    weights.name = "weights"
    # the way that these weights are returned by my wrapper strips away their 
    # dimension names, and also with the data in teh wrong shape. Rather annoying. 
    # Next few lines here are to reintroudce the dimension names and transpose
    # to a shape that matches the data....
    weights = weights.rename({'ncl2':'time', 'ncl1':'lev', 'ncl3':'lat'})
    weights = weights.reindex(indexers={'time':e90_trop.time.values, 
                                        'lev':e90_trop.lev.values,
                                        'lat':e90_trop.lat.values})
    e90_trop_weighted   = e90_trop.weighted(weights)
    e90_strat_weighted  = e90_strat.weighted(weights)
    st80_trop_weighted  = st80_trop.weighted(weights)
    st80_strat_weighted = st80_strat.weighted(weights)
    
    # now do the same thing for the climatological means...
    weights = wrappers.dpres_hybrid_ccm(PS_climo, P0, 
                                        hyai, hybi).transpose('ncl2', 'ncl1', 'ncl3')
    weights.name = "weights"
    weights = weights.rename({'ncl2':'month', 'ncl1':'lev', 'ncl3':'lat'})
    weights = weights.reindex(indexers={'month': e90_climo_trop.month.values, 
                                        'lev':e90_climo_trop.lev.values,
                                        'lat':e90_climo_trop.lat.values})
    e90_climo_trop_weighted   = e90_climo_trop.weighted(weights)
    e90_climo_strat_weighted  = e90_climo_strat.weighted(weights)
    st80_climo_trop_weighted  = st80_climo_trop.weighted(weights)
    st80_climo_strat_weighted = st80_climo_strat.weighted(weights)
    
else:
    e90_trop_weighted, e90_strat_weighted = e90_trop, e90_strat
    st80_trop_weighted, st80_strat_weighted = st80_trop, st80_strat
    e90_climo_trop_weighted, e90_climo_strat_weighted = \
        e90_climo_trop, e90_climo_strat
    st80_climo_trop_weighted, st80_climo_strat_weighted = \
        st80_climo_trop, st80_climo_strat

# done computing vertical weights, do vertical averaging...
e90_trop_mean   = e90_trop_weighted.mean('lev')
e90_strat_mean  = e90_strat_weighted.mean('lev')
st80_trop_mean  = st80_trop_weighted.mean('lev')
st80_strat_mean = st80_strat_weighted.mean('lev')
e90_climo_trop_mean   = e90_climo_trop_weighted.mean('lev')
e90_climo_strat_mean  = e90_climo_strat_weighted.mean('lev')
st80_climo_trop_mean  = st80_climo_trop_weighted.mean('lev')
st80_climo_strat_mean = st80_climo_strat_weighted.mean('lev')

# ---- lat weighting
if(WEIGHT_LAT):
    # weight the datasets in laitutde by cos(lat)...
    # same procedure as for the weighted tropical averages above...
    print('---- weighting lat...')
    weights = np.cos(np.deg2rad(lat))
    weights.name = "weights"
    e90_trop_weighted   = e90_trop_mean.weighted(weights)
    e90_strat_weighted  = e90_strat_mean.weighted(weights)
    st80_trop_weighted  = st80_trop_mean.weighted(weights)
    st80_strat_weighted = st80_strat_mean.weighted(weights)
    e90_climo_trop_weighted   = e90_climo_trop_mean.weighted(weights)
    e90_climo_strat_weighted  = e90_climo_strat_mean.weighted(weights)
    st80_climo_trop_weighted  = st80_climo_trop_mean.weighted(weights)
    st80_climo_strat_weighted = st80_climo_strat_mean.weighted(weights)
else:
    e90_trop_weighted, e90_strat_weighted = e90_trop_mean, e90_strat_mean
    st80_trop_weighted, st80_strat_weighted = st80_trop_mean, st80_strat_mean
    e90_climo_trop_weighted, e90_climo_strat_weighted = \
        e90_climo_trop_mean, e90_climo_strat_mean
    st80_climo_trop_weighted, st80_climo_strat_weighted = \
        st80_climo_trop_mean, st80_climo_strat_mean

# done computing horizontal weights, do horizontal averaging...
e90_trop_mean   = e90_trop_weighted.mean('lat')
e90_strat_mean  = e90_strat_weighted.mean('lat')
st80_trop_mean  = st80_trop_weighted.mean('lat')
st80_strat_mean = st80_strat_weighted.mean('lat')
e90_climo_trop_mean   = e90_climo_trop_weighted.mean('lat')
e90_climo_strat_mean  = e90_climo_strat_weighted.mean('lat')
st80_climo_trop_mean  = st80_climo_trop_weighted.mean('lat')
st80_climo_strat_mean = st80_climo_strat_weighted.mean('lat')


# --------------- get anomalies ---------------
# finally, compute monthly anomalies from the climatology, if 
# we are not anayzing annual data. 
# This is done for the tropically-averaged tracer fields and tropopause
# position, as well as the globally-avergad stratosphere and troposphere
# tracers

if(not DO_ANNUAL):
    print('---- computing tropical mean anomalies...')
    e90_eq_mean_anom  = e90_eq_mean.groupby('time.month') - e90_climo_eq_mean
    st80_eq_mean_anom = st80_eq_mean.groupby('time.month') - st80_climo_eq_mean
    
    # because of the way we built these tropopause arrays above, they cannot 
    # be simply grouped by month with xarray... here's an ugly manual solution
    # which just repeats the monthly climatology over a new dataset that fills
    # and aligns with the time domain of the data
    trop_t_climo_eq_mean_full = trop_t_eq_mean.copy()
    for i in range(len(trop_t_eq_mean_time)):
        year = trop_t_eq_mean_time[i]
        month = np.round((year - np.floor(year)) * 12).astype(int)
        trop_t_climo_eq_mean_full[i] = trop_t_climo_eq_mean[month]
    trop_t_eq_mean_anom = trop_t_eq_mean - trop_t_climo_eq_mean_full
    
    print('---- computing global mean anomalies...')
    e90_trop_mean_anom = e90_trop_mean.groupby('time.month') - e90_climo_trop_mean
    e90_strat_mean_anom = e90_strat_mean.groupby('time.month') - e90_climo_strat_mean
    st80_trop_mean_anom = st80_trop_mean.groupby('time.month') - st80_climo_trop_mean
    st80_strat_mean_anom = st80_strat_mean.groupby('time.month') - st80_climo_strat_mean


# ----------------------------------------------------------------------
# ------------------------------ plotting ------------------------------

# --------------- organize historical eruption data ---------------
# These data are moslty used for visual purposes, and will mark when eruptions occur
# int he historical time domain, and visually represent their injection magnotudes
print('---- plotting...')
# from Toohey+ EVA
eruptions = {'Agung':          {'year':1963, 'month':3,  'TgS':5.22},
             'Fuego':          {'year':1974, 'month':10,  'TgS':1.18},
             'El Chichon':     {'year':1982, 'month':4, 'TgS':3.5},
             'Pinatubo':       {'year':1991, 'month':6,  'TgS':9},
             'Manam':          {'year':2005, 'month':1,  'TgS':0.08},
             'Soufriere Hills':{'year':2006, 'month':5,  'TgS':0.07},
             'Rabaul':         {'year':2006, 'month':10,  'TgS':0.08},
             'Kasatochi':      {'year':2008, 'month':8, 'TgS':0.19},
             'Sarychev':       {'year':2009, 'month':6,  'TgS':0.28},
             'Merapi':         {'year':2010, 'month':11,  'TgS':0.05},
             'Nabro':          {'year':2011, 'month':6, 'TgS':0.184}}
# if we are analyzing monthly rather than annual data, for convenice
# we will update the "year" field of this eruption log to be a decimal
# which includes the month

emissions_file = '/pscratch/sd/j/jhollo/E3SM/historical/CLDERA_2ndHistorical_1950-2014_monthly/emissions'\
                 '/VolcanEESMv3.11_SO2_850-2016_Mscale_Zreduc_2deg_c180812.nc'
emissions_data = xr.open_dataset(emissions_file)
pdb.set_trace()



if(not DO_ANNUAL):
    for eruption in eruptions.keys():
        e = eruptions[eruption]
        eruptions[eruption]['year'] = e['year'] + e['month']/12

# now we build a dictionary of eruption magnitudes per-year. Most of these
# years will hold a zero, and years with multiple eruptions will hold
# the summed injection magnitude
eruption_mags = {year:0 for year in e90_trop_mean.year.values}
for eruption in eruptions.values():
    eruption_mags[eruption['year']] = eruption_mags[eruption['year']] + eruption['TgS']
    
    
# --------------- make anual/monthly plotting choices ---------------
# depending if we are analyzing annual or monthly data, we will have different 
# settings for colormaps, axis labels, etc.
# The variables named like "st80_trop_line_data", etc. will set which data to
# use for the "global tracer evolution" plot. If DO_ANNUAL, then this will be
# average tracer concentrations in ppb. If not, then this will be average tracer
# anomalies from the monthly climatology, in ppb.

if(DO_ANNUAL):
    colors = e90_trop_mean.year.values
    cmap = plt.cm.YlOrRd
    vmin, vmax = 1850, 2014
    mew = 0
    cmap_label = 'year'
    st80_trop_line_data, st80_strat_line_dat = st80_trop_mean, st80_strat_mean
    e90_trop_line_data, e90_strat_line_data  = e90_trop_mean, e90_strat_mean
    e90_label  = 'E90 global strat. mean [ppb]'
    st80_label = 'ST80 global trop. mean [ppb]'
else:
    colors = ((e90_strat_mean.year.values - 
               np.floor(e90_strat_mean.year.values)) * 12).astype(int)
    cmap = plt.cm.twilight_shifted
    cmap_label = 'month'
    vmin, vmax = 0, 11
    mew = 1.5
    st80_trop_line_data, st80_strat_line_dat = st80_trop_mean_anom, st80_strat_mean_anom
    e90_trop_line_data, e90_strat_line_data  = e90_trop_mean_anom, e90_strat_mean_anom 
    e90_label  = 'E90 global strat. mean anomaly [ppb]'
    st80_label = 'ST80 global trop. mean anomaly [ppb]'
    
    
# --------------- create figure ---------------
# Figure will have three rows. Top row will contin two tracer-tracer correlation plots
# side-by-side for E90, ST80 in the stratosphere, troposphere.
# Second row will havfe the "Global tracer evolution" plot, for which we selected data 
# above. 
# Third row will have a "tropopause evolution" plot, which shows the tropical E90 90ppb 
# contour, and mean tropical tropopause height for the time domain. If not DO_ANNUAL, then
# this will be shown as their anomalies from the monthly climatology.

fig = plt.figure(figsize=(10,12))
gs  = gridspec.GridSpec(3, 2)
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
ax3 = fig.add_subplot(gs[1, :])
ax4 = fig.add_subplot(gs[2, :])


# ---------------------------------------------------------------
# --------------- tracer-tracer correlation plots ---------------

# scatter marker area for years with no eruption
s_area = 1
# scale scatter marker areas with eruption magnitude
s_erupt_areas = (np.array(list(eruption_mags.values()))+s_area)**2.5*np.pi
# mask for years with/without an eruption
emask = s_erupt_areas > np.pi

# e90-st80 stratopshere tracer-tracer correlation. Non-eruption years marked with dots
ax1.scatter(e90_strat_mean, st80_strat_mean, marker='o', s=s_area**2*np.pi, 
            c=colors, cmap=cmap)
# Eruption years marked with stars, scaled for size
ax1.scatter(e90_strat_mean[emask], st80_strat_mean[emask], marker='*', edgecolors='lime',
            s=s_erupt_areas[emask], c=colors[emask], cmap=cmap, linewidth=mew, 
            vmin=vmin, vmax=vmax)
ax1.set_title('STRATOSPHERE')
ax1.set_xlabel('e90 global mean [ppb]')
ax1.set_ylabel('st80 global mean [ppb]')
# colorbar axis. This won't be used, and is just a padding to make sure both the strat 
# and trop correlation axes are the same size
divider = make_axes_locatable(ax1)
cax1 = divider.append_axes("right", size="5%", pad=0.05)
cax1.axis('off')

# e90-st80 tropopshere tracer-tracer correlation, point styling same as above
im = ax2.scatter(e90_trop_mean, st80_trop_mean, marker='o', s=s_area**2*np.pi, 
            c=colors, cmap=cmap)
ax2.scatter(e90_trop_mean[emask], st80_trop_mean[emask], marker='*', edgecolors='lime',
            s=s_erupt_areas[emask], c=colors[emask], cmap=cmap, linewidth=mew, 
            vmin=vmin, vmax=vmax)
ax2.set_title('TROPOSPHERE')
ax2.set_xlabel('e90 global mean [ppb]')
# colorbar axis
divider = make_axes_locatable(ax2)
cax2 = divider.append_axes("right", size="5%", pad=0.05)
# colorbar
cb = plt.colorbar(im, cax=cax2)
cb.set_label(cmap_label)


# -------------------------------------------------------
# --------------- global tracer evolution ---------------

# E90 strat mean on the left axis (black), ST80 trop mean on the right axis (grey)
ax33 = ax3.twinx()
ax33.plot(st80_trop_mean.year.values, st80_trop_line_data, '-k', alpha=0.44, lw=0.75)
ax3.plot(e90_strat_mean.year.values, e90_strat_line_data, '-k', lw=0.75)
ax3.set_title('Global tracer evolution')
ax3.set_xlabel('time')
ax3.set_ylabel(e90_label)
ax33.set_ylabel(st80_label)
ax33.yaxis.label.set_color('grey')

# draw vertical lines on eruption years, with the thickness corresponding to the 
# injection magnitude
for year in eruption_mags.keys():
    if(eruption_mags[year] > 0):
        ax3.axvline(year, color='r', lw=(eruption_mags[year]/9)*2.5)
ypos = np.mean(e90_strat_line_data)
# dummy plot for legend
ax3.plot([1980,1980], [ypos,ypos], color='r', label='eruption date\n(line width = TgS)')
ax3.legend(fancybox=False, loc='lower left')

# these were eyeballed and work well for this historical data
if(not DO_ANNUAL):
    ax3.set_xlim([1960, 2000])
    ax3.set_ylim([-0.4, 0.4])
    ax33.set_ylim([-0.3, 0.3])
    
# ------------------------
# ---- e90 90ppb evolution

# if we're not analyzing the annual data, then we will plot anomalies 
# here, and show them over a righter time window (near Pinatubo)
if(not DO_ANNUAL):
    
    # this is going to be a bit awkard...
    # we want to plot the anomalous 90ppb line, but this line is only
    # obtained by contouring the 2D lev-time distribution. Our
    # strategy will be do then do this twice for the tropically avergaed data, 
    # and the tropically averaged climatology. This will again require us to 
    # manually repeat the climatological data on a monthly basis over the full
    # time series to later take a difference (since we will not be able to use
    # xarray's by-month grouping abilities on the return object from ax.contour
    # First, let's build this repeated-climatology array
    e90_climo_eq_mean_full = e90_eq_mean.copy(deep=True)
    for i in range(len(e90_eq_mean.year.values)):
        year = e90_eq_mean.year.values[i]
        month = np.round((year - np.floor(year)) * 12).astype(int)
        e90_climo_eq_mean_full[i, :] = e90_climo_eq_mean[month, :]
    
    # now contour the tropical mean e90, and tropical mean climatology
    # clear the result from the axis since we only want the return (yes I know 
    # I could probably use numpy instead of pyplot for this)
    YEAR, LEV = np.meshgrid(e90_eq_mean.year.values, lev)
    cf  = ax4.contour(YEAR, LEV, e90_eq_mean.T, levels=[90])
    cfc = ax4.contour(YEAR, LEV, e90_climo_eq_mean_full.T, levels=[90], linestyles=['--'])
    ax4.clear()

    # extract the data for the single contour that we plotted above. This is acessible through
    # the "allsegs" attribute of the returned contour object from the plotting routing. Used
    # the debugger to see how to index and transpose this object. Do this for both the tropical
    # mean data and the climatology...
    ppb_mean  = cf.allsegs[0][0].T[1]
    ppb_mean_time = cf.allsegs[0][0].T[0]
    ppb_climo = cfc.allsegs[0][0].T[1]
    ppb_climo_time = cfc.allsegs[0][0].T[0]
    # Now, it seems that the resulting time positions for the contour data are not equal between
    # the data and the climatology. I don't know why, but plt.contour must be choosing these
    # positions by some method. The positions do share endpoints, and they are close enough that
    # we will just interpolate the climatological contour to the posiiton of the data...
    # (use a debugger to inspect the objects created above and verify that this claim is indeed
    # true). First obtain the interpolating function:
    ff = scipy.interpolate.interp1d(ppb_climo_time, ppb_climo)
    
    # now evaluate the interpolating function...
    #need the [:-1] here since we cant interp to the endpoint of the original dataset
    ppb_climo = ff(ppb_mean_time[:-1])
    ppb_climo_time = ppb_mean_time[:-1]
    ppb_mean_time = ppb_mean_time[:-1]
    ppb_mean = ppb_mean[:-1]
    
    # finally, take the anomaly
    ppb_anom = ppb_mean - ppb_climo

    # plot the result
    # We show botht the e90 90ppb position anomaly (blue), and the tropopause anoamly (grey)
    ax4.plot(ppb_mean_time, ppb_anom, '-b', lw=0.75, 
             label='e90 tropical mean (+- 20 deg) 90 ppb position anomaly')
    ax4.plot(trop_t_eq_mean_time, trop_t_eq_mean_anom, '-k', lw=0.75, alpha=0.5, 
            label='tropical mean (+- 20 deg) tropopause position anomaly')
    ax4.set_title('Tropopause evolution')
    ax4.set_xlabel('time [years]')
    ax4.set_ylabel('lev [hPa]')
    ax4.legend(fancybox=False, loc='lower left')
else:
    # if we are analyzing annual data, then this is easier. Don't compute anomalies, just plot
    # the e90 90ppb contour (blue) and tropopause position (grey)
    YEAR, LEV = np.meshgrid(e90_eq_mean.year.values, lev)
    cf  = ax4.contour(YEAR, LEV, e90_eq_mean.T, levels=[90], colors='cornflowerblue')
    ax4.plot([1900,1900], [100,100], color='cornflowerblue', 
             label='e90 tropical mean (+- 20 deg) 90 ppb position')
    ax4.invert_yaxis()
    ax4.set_title('Tropopause evolution')
    ax4.set_xlabel('time [years]')
    ax4.set_ylabel('lev [hPa]')

    ax4.plot(trop_t_eq_mean_time, trop_t_eq_mean, '-k', alpha=0.5, 
             label='tropical mean (+- 20 deg) tropopause position')
    ax4.legend(fancybox=False, loc='upper left')
    ax4.set_ylim([94, 112])

# draw vertical lines on eruption years, with the thickness corresponding to the 
# injection magnitude
for year in eruption_mags.keys():
    if(eruption_mags[year] > 0):
        ax4.axvline(year, color='r', lw=(eruption_mags[year]/9)*2.5)
# dummy plot for legend
ax4.plot([1980,1980], [0,0], color='r', label='eruption date\n(line width = TgS)')

# if we aren't analyzing the annual data, then zoom these fgures into the volcanically active 
# period contining Pinatubio and El Chichon
if(not DO_ANNUAL):
    ax4.set_xlim([1960, 2000])


# -------------------------------------------
# --------------- save figure ---------------
    
plt.tight_layout()
plt.savefig('SCATTER_{}_wlev{}_wlat{}.png'.format(['monthly', 'annual'][int(DO_ANNUAL)], int(WEIGHT_LEV), int(WEIGHT_LAT)))

