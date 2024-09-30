import sys
import pdb
import glob
import scipy
import numpy as np
import pandas as pd
import xarray as xr
import seaborn as sns
from scipy.stats import t
from matplotlib import colors
import matplotlib.pyplot as plt
from scipy.stats import linregress
from scipy.interpolate import interp1d
import matplotlib.gridspec as gridspec
from matplotlib.ticker import ScalarFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
import wrappers
import climate_toolbox as ctb

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
WEIGHT_LAT = True       # whether or not to weight horizontal averages in by cos(lat)
WEIGHT_LEV = True       # whether or not to weight veritcal averages by pdel
DO_ANNUAL  = True       # whether or not to use annually-averaged data fror this analysis
                        # if true, climatological anomalies will not be taken in E90, ST80, 
                        # or the tropopause; averaged raw ppb values will be plotting
                        # instead
USE_EVA    = False      # whether or not to use the minimal eruption dataset manually
                        # collected
                        # from Toohey+ (Easy Volcanic Aerosols, EVAv1.0). If False, will use 
                        # the much more expansive table collected from the CLDERA E3SM 
                        # emissions file metadata. Should be False.
PLOT_ERUPTIONS = True   # whether or not to display vertical lines on figures giving 
                        # times and strengths of eruptions
PLOT_TRENDS    = True   # whether or not to plot best-fit lines
ERUPT_FIT_MONTHS = 4    # number of months to use for post-eruption regression
TROPICAL_MEANS = True   # whether or not to restrict measures of the mean strat/trop tracer content
                        # to the tropics
MIN_LEV = 50            # "global" tracer distributions only considered up to MIN_LEV hPa
LEV_SLICE = slice(MIN_LEV,1200)
TROP_BOUND = 20         # the "tropical region" is defined on +- (TROP_BOUND) degrees in latitude
TROP_SLICE = slice(-TROP_BOUND,TROP_BOUND)

REG_CI = 95            # confidence interval to estiamte via bootstrap for regression (%)
REG_NBOOT = 10000      # number of bootstrap samples to draw for regressions ci
REG_MIN_TGS = 0.5      # minimum eruption mag (TgS) to include in regression

# --------------- organize historical eruption data ---------------
# These data are moslty used for visual purposes, and will mark when eruptions occur
# in the historical time domain, and visually represent their injection magnotudes

print('---- reading emissions data...')
if(USE_EVA):
    # This table manually entered from Toohey+ EVA. E3SM data read in in the section below
    # willbe parsed to filla dictionary with the same form as this one
    eruptions = {'Agung':          {'year':1963, 'month':3,  'TgS':5.22,  'lat':-8.342},
                 'Fuego':          {'year':1974, 'month':10, 'TgS':1.18,  'lat':14.473},
                 'El Chichon':     {'year':1982, 'month':4,  'TgS':3.5,   'lat':17.360},
                 'Pinatubo':       {'year':1991, 'month':6,  'TgS':9,     'lat':15.130},
                 'Manam':          {'year':2005, 'month':1,  'TgS':0.08,  'lat':-4.080},
                 'Soufriere Hills':{'year':2006, 'month':5,  'TgS':0.07,  'lat':16.720},
                 'Rabaul':         {'year':2006, 'month':10, 'TgS':0.08,  'lat':-4.270},
                 'Kasatochi':      {'year':2008, 'month':8,  'TgS':0.19,  'lat':52.180},
                 'Sarychev':       {'year':2009, 'month':6,  'TgS':0.28,  'lat':48.090},
                 'Merapi':         {'year':2010, 'month':11, 'TgS':0.05,  'lat':-7.540},
                 'Nabro':          {'year':2011, 'month':6,  'TgS':0.184, 'lat':13.370}}
else:
    # This is the volcanic SO2 emissions file used for CLDERA E3SM coupled runs. 
    # Parsing code here puts the data into the same format as seen for the EVA 
    # 'eruptions' dictionary above
    emissions_file = '/pscratch/sd/j/jhollo/E3SM/historical/CLDERA_2ndHistorical_1950-2014_monthly/'\
                     'emissions/VolcanEESMv3.11_SO2_850-2016_Mscale_Zreduc_2deg_c180812.nc'
    emissions = xr.open_dataset(emissions_file, decode_times=False)
    emissions_data = emissions.data_summary.split('\n')[53:]
    eruptions = {}
    for entry in emissions_data:
        eruption_dict = {}
        entry_data = entry.split()

        entry_name   = ' '.join(entry_data[9:])
        entry_year   = int(entry_data[0][:4])
        entry_altmin = float(entry_data[3])
        if(entry_year < 1850 or entry_year >= 2015): continue
        if(entry_altmin < 10): continue

        eruption_dict['year']   = entry_year
        eruption_dict['month']  = int(entry_data[0][4:6])
        eruption_dict['TgS']    = float(entry_data[5])
        eruption_dict['lat']    = float(entry_data[1])
        eruption_dict['AltMin'] = entry_altmin
        eruptions[entry_name]   = eruption_dict

# get the maximum emission size for later normalization
max_TgS = max([eruptions[e]['TgS'] for e in eruptions.keys()])

# if we are analyzing monthly rather than annual data, for convenice
# we will update the "year" field of this eruption log to be a decimal
# which includes the month
if(not DO_ANNUAL):
    for eruption in eruptions.keys():
        e = eruptions[eruption]
        eruptions[eruption]['year'] = e['year'] + e['month']/12

        
# --------------- read data ---------------

print('---- reading data...')
pdir         = '/pscratch/sd/j/jhollo/E3SM/historical/CLDERA_2ndHistorical_1950-2014_monthly'
# full data (we'll pull meta data, coordiante info, etc from here)
datf    = '{}/histoircal_h0_PS_1850-2014.regrid.fv0.9x1.25_bilinear.nc'.format(pdir)
# zonal mean data
zm_E90  = '{}/histoircal_h0_E90j_1850-2014.regrid.fv0.9x1.25_bilinear.zonalmean.nc'.format(pdir)
zm_ST80 = '{}/histoircal_h0_ST80_25j_1850-2014.regrid.fv0.9x1.25_bilinear.zonalmean.nc'.format(pdir)
zm_T    = '{}/histoircal_h0_T_1850-2014.regrid.fv0.9x1.25_bilinear.zonalmean.nc'.format(pdir)
zm_PS   = '{}/histoircal_h0_PS_1850-2014.regrid.fv0.9x1.25_bilinear.zonalmean.nc'.format(pdir)

# get coordinates etc.
dat                        = xr.open_dataset(datf)
lat, lev, ilev             = dat['lat'], dat['lev'], dat['ilev']
P0, hyam, hybm, hyai, hybi = dat['P0'], dat['hyai'], dat['hybi'], dat['hyai'], dat['hybi']

# get zonally averaged data
e90   = xr.open_dataset(zm_E90)['E90j']      # in ppb
st80  = xr.open_dataset(zm_ST80)['ST80_25j'] # in ppb
T     = xr.open_dataset(zm_T)['T']           # in K
PS    = xr.open_dataset(zm_PS)['PS']         # in Pa

# get gridpoint pressure
P = ctb.compute_hybrid_pressure(dat['PS'], hyam, hybm, dims_like=T)
pdb.set_trace()


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
time = e90.year.values


# --------------- yearly eruption history ----------
# first we build a dictionary of eruption magnitudes per-year. Most of these
# years will hold a zero, and years with multiple eruptions will hold
# the summed injection magnitude
eruption_mags = {year:0 for year in time}
for eruption in eruptions.values():
    eruption_mags[eruption['year']] = eruption_mags[eruption['year']] + eruption['TgS']
# now lets build a normalized version of the same list for later setting plot alpha values
eruption_mags_normed = dict(eruption_mags)
for year in time:
    eruption_mags_normed[year] = eruption_mags[year] / max_TgS
    
# let's also build a filtered version of these eruption magnitudes for tropical volcanoes
eruption_mags_eq = {year:0 for year in time}
for eruption in eruptions.values():
    if(abs(eruption['lat']) < 20):
        eruption_mags_eq[eruption['year']] = \
            eruption_mags_eq[eruption['year']] + eruption['TgS']
    else:
        continue
eruption_mags_eq_normed = dict(eruption_mags_eq)
for year in time:
    eruption_mags_eq_normed[year] = eruption_mags_eq[year] / max_TgS
    
    
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
eq_mask = np.logical_and(lat > -TROP_BOUND, lat < TROP_BOUND)
trop_t_eq_mean = np.mean(trop_t.T[eq_mask], axis=0)
trop_t_climo_eq_mean = np.mean(trop_t_climo.T[eq_mask], axis=0)


# --------------- get e90, st80 in troposphere, stratosphere ---------------
# do this by simply applying the masks computed above. Array positions removed 
# by the mask will be filled wiht zeros

print('---- masking tracer fields...')
e90_trop         = e90.where(trop_mask, 0)
e90_strat        = e90.where(strat_mask, 0)
st80_trop        = st80.where(trop_mask, 0)
st80_strat       = st80.where(strat_mask, 0)
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

    
# --------------- spatially reduce data ---------------
# reduce datasets to below 10hPa, and to the tropical region if desired

if(TROPICAL_MEANS):
    reduction = {'lev':LEV_SLICE, 'lat':TROP_SLICE}
    lat       = lat.sel({'lat':TROP_SLICE})
    PS        = PS.sel({'lat':TROP_SLICE})
    PS_climo  = PS_climo.sel({'lat':TROP_SLICE})
else:
    reduction = {'lev':LEV_SLICE}
lev  = lev.sel({'lev':LEV_SLICE})

# get ILEV_SLICE corresponding to LEV_SLICE
ILEV_ISEL = [i for i in range(len(ilev)) if ~np.isnan(ilev.where(ilev>=MIN_LEV))[i]]
ILEV_ISEL.insert(0, ILEV_ISEL[0]-1)
hyai = hyai.isel({'ilev':ILEV_ISEL})
hybi = hybi.isel({'ilev':ILEV_ISEL})

e90              = e90.sel(reduction)
st80             = st80.sel(reduction)
T                = T.sel(reduction)
e90_climo        = e90_climo.sel(reduction)
st80_climo       = st80_climo.sel(reduction)
e90_trop         = e90_trop.sel(reduction)
e90_strat        = e90_strat.sel(reduction)
st80_trop        = st80_trop.sel(reduction)
st80_strat       = st80_strat.sel(reduction)
e90_climo_trop   = e90_climo_trop.sel(reduction)
e90_climo_strat  = e90_climo_strat.sel(reduction)
st80_climo_trop  = st80_climo_trop.sel(reduction)
st80_climo_strat = st80_climo_strat.sel(reduction)
    
    
# --------------- take weighted tropical averages ---------------
# here we compute tropical-mean averages for the tracers and their 
# climatologies. This will result in 2D time series in (lev-height).
# Contouring these later will give e.g. mean tropical E90 90ppb position

print('---------- computing tropical means...')
# ---- lat weighting
e90_eq        = e90.sel({'lat':TROP_SLICE})
st80_eq       = st80.sel({'lat':TROP_SLICE})
e90_climo_eq  = e90_climo.sel({'lat':TROP_SLICE})
st80_climo_eq = st80_climo.sel({'lat':TROP_SLICE})
if(WEIGHT_LAT):
    # weight the datasets in laitutde by cos(lat)...
    print('---- weighting lat...')
    weights = np.cos(np.deg2rad(lat))
    weights.name = "weights"
    e90_eq_weighted        = e90_eq.weighted(weights)
    st80_eq_weighted       = st80_eq.weighted(weights)
    e90_climo_eq_weighted  = e90_climo_eq.weighted(weights)
    st80_climo_eq_weighted = st80_climo_eq.weighted(weights)
else:
    e90_eq_weighted        = e90_eq
    st80_eq_weighted       = st80_eq
    e90_climo_eq_weighted  = e90_climo_eq
    st80_climo_eq_weighted = st80_climo_eq
    
# done computing horizontal weights, do horizontal averaging...
e90_eq_mean         = e90_eq_weighted.mean('lat')
st80_eq_mean        = st80_eq_weighted.mean('lat')
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
    # dimension names, and also with the data in the wrong shape. Rather annoying. 
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


# --------------- e90 90ppb time series ---------------
METHOD = 'linear'
if(not DO_ANNUAL):
    # if we're not analyzing the annual data, then we will plot anomalies e90 90 ppb and
    # tropopause anomalies, rather than the raw time series, and show them over a tighter 
    # time window (near Pinatubo)
    
    # Interpolate along lev for time
    e90_90ppb_eq_lev = xr.zeros_like(e90_trop_mean)
    e90_90ppb_eq_lev.name = 'E90j_90ppb_lev'
    for i in range(len(time)):
        interp_func = interp1d(e90_eq_mean[i,:], lev, kind=METHOD) 
        e90_90ppb_eq_lev[i] = float(interp_func(90))

    e90_climo_90ppb_lev = xr.zeros_like(e90_climo_trop_mean)
    e90_climo_90ppb_lev.name = 'E90j_climo_90ppb_lev'
    for i in range(12):
        interp_func = interp1d(e90_climo_eq_mean[i,:], lev, kind=METHOD) 
        e90_climo_90ppb_lev[i] = float(interp_func(90))

else:
    # if we are analyzing annual data, don't compute anomalies, just plot the e90 90ppb
    # contour and tropopause position
    
    # Interpolate along lev for time
    e90_90ppb_eq_lev = xr.zeros_like(e90_trop_mean)
    e90_90ppb_eq_lev.name = 'E90j_90ppb_lev'
    for i in range(len(time)):
        interp_func = interp1d(e90_eq_mean[i,:], lev, kind=METHOD) 
        e90_90ppb_eq_lev[i] = float(interp_func(90))
        

# --------------- get anomalies ---------------
# finally, compute monthly anomalies from the climatology, if 
# we are not anayzing annual data. 
# This is done for the tropically-averaged tracer fields and tropopause
# position, as well as the globally-avergad stratosphere and troposphere
# tracers

if(not DO_ANNUAL):
    print('---- computing tropical mean anomalies...')
    
    # get 2D anomalies
    e90_eq_mean_anom  = e90_eq_mean.groupby('time.month') - e90_climo_eq_mean
    st80_eq_mean_anom = st80_eq_mean.groupby('time.month') - st80_climo_eq_mean
    
    # get 90ppb anomalies
    e90_90ppb_eq_lev_anom = e90_90ppb_eq_lev.groupby('time.month') - e90_climo_90ppb_lev
    
    # because of the way we built these tropopause arrays above, they cannot 
    # be simply grouped by month with xarray... here's an ugly manual solution
    # which just repeats the monthly climatology over a new dataset that fills
    # and aligns with the time domain of the data
    trop_t_climo_eq_mean_full = trop_t_eq_mean.copy()
    for i in range(len(time)):
        year = time[i]
        month = np.round((year - np.floor(year)) * 12).astype(int)
        trop_t_climo_eq_mean_full[i] = trop_t_climo_eq_mean[month]
    trop_t_eq_mean_anom = trop_t_eq_mean - trop_t_climo_eq_mean_full
    
    print('---- computing global mean anomalies...')
    e90_trop_mean_anom = e90_trop_mean.groupby('time.month') - e90_climo_trop_mean
    e90_strat_mean_anom = e90_strat_mean.groupby('time.month') - e90_climo_strat_mean
    st80_trop_mean_anom = st80_trop_mean.groupby('time.month') - st80_climo_trop_mean
    st80_strat_mean_anom = st80_strat_mean.groupby('time.month') - st80_climo_strat_mean


# --------------- trend analysis  ---------------

if(DO_ANNUAL):
    print('---- doing global linear regression...')
    # Here we compute linear fits for the e90 90 ppb line and tropopause position.
    # For the annual analysis, this is done on the raw time series data for each line.

    # get best fit, std error, r- and p-values of e90 global strat. mean and 
    # st80 global trop. mean
    e90_strat_fit = linregress(time, e90_strat_mean)
    st80_trop_fit = linregress(time, st80_trop_mean)
    # get best fit, std error, r- and p-values of e90 90ppb altitude and 
    # tropopause altitude
    e90_90ppb_eq_fit = linregress(time, e90_90ppb_eq_lev)
    trop_eq_fit = linregress(time, trop_t_eq_mean)
    
if(not DO_ANNUAL):
    print('---- doing per-eruption linear regression...')
    # Here we compute trends for N-months post-eruption for all tropical eruptions
    # This is only done for the monthly analysis
    
    e90_strat_anom_eruption_fits    = []
    st80_trop_anom_eruption_fits    = []
    e90_90ppb_eq_anom_eruption_fits = []
    trop_eq_anom_eruption_fits      = []
    eruption_fit_mags               = []
    eruption_eq_fit_mags            = []
    eruption_idx                    = []
    eruption_eq_idx                 = []
    ERM = ERUPT_FIT_MONTHS

    for i in range(len(eruption_mags)):
        year = list(eruption_mags.keys())[i]
        mag  = eruption_mags[year]
        
        if(mag > 0):
            eruption_idx.append(i)
            eruption_fit_mags.append(mag)
            # get best fit, std error, r- and p-values of e90 global strat. mean and 
            # st80 global trop.mean post-eruptions
            e90_strat_anom_eruption_fits.append(linregress(time[i:i+ERM], 
                                                           e90_strat_mean_anom[i:i+ERM]))
            st80_trop_anom_eruption_fits.append(linregress(time[i:i+ERM], 
                                                           st80_trop_mean_anom[i:i+ERM]))
            
            # get best fit, std error, r- and p-values of e90 90ppb altitude and 
            # tropopause altitude post-
            erupt_lat  = [e['lat'] for e in list(eruptions.values()) if e['year'] == year][0]
            if(abs(erupt_lat) <= 20):
                eruption_eq_idx.append(i)
                eruption_eq_fit_mags.append(mag)
                e90_90ppb_eq_anom_eruption_fits.append(linregress(time[i:i+ERM], 
                                                                  e90_90ppb_eq_lev_anom[i:i+ERM]))
                trop_eq_anom_eruption_fits.append(linregress(time[i:i+ERM], 
                                                             trop_t_eq_mean_anom[i:i+ERM]))
                
    eruption_fit_mags    = np.array(eruption_fit_mags)
    eruption_eq_fit_mags = np.array(eruption_eq_fit_mags)
        

# ----------------------------------------------------------------------
# ------------------------------ plotting ------------------------------

print('---- plotting...')
# --------------- make anual/monthly plotting choices ---------------
# depending if we are analyzing annual or monthly data, we will have different 
# settings for colormaps, axis labels, etc.
# The variables named like "st80_trop_line_data", etc. will set which data to
# use for the "global tracer evolution" plot. If DO_ANNUAL, then this will be
# average tracer concentrations in ppb. If not, then this will be average tracer
# anomalies from the monthly climatology, in ppb.

st80_color     = 'yellowgreen'
e90_color      = 'cornflowerblue'
trop_color     = 'grey'
eruption_color = 'r'

if(DO_ANNUAL):
    colors = time
    cmap = plt.cm.YlOrRd
    vmin, vmax = 1850, 2014
    mew = 0
    cmap_label = 'year'
    st80_trop_line_data, st80_strat_line_dat = st80_trop_mean, st80_strat_mean
    e90_trop_line_data, e90_strat_line_data  = e90_trop_mean, e90_strat_mean
    # scatter marker area for years with no eruption
    s_area = 2
    scatter_label = 'annual means'
    if(PLOT_TRENDS):
        linestyle = '-'
    else:
        linestyle = '-'
else:
    colors = ((time - np.floor(time)) * 12).astype(int)
    cmap = plt.cm.twilight_shifted
    cmap_label = 'month'
    vmin, vmax = 0, 11
    mew = 1.5
    st80_trop_line_data, st80_strat_line_dat = st80_trop_mean_anom, st80_strat_mean_anom
    e90_trop_line_data, e90_strat_line_data  = e90_trop_mean_anom, e90_strat_mean_anom 
    # scatter marker area for years with no eruption
    s_area = 1
    scatter_label = 'monthly means'
    linestyle = '-'
    
# make choices for plot titles, labels, etc. depending on whether we're restricitng total
# strat/trop tracers mixing ratios to the tropics
if(TROPICAL_MEANS):
    tracertracer_strat_title = 'TROPICAL STRATOSPHERE'
    tracertracer_trop_title  = 'TROPICAL TROPOSPHERE'
    tracertracer_xlabel      = 'E90 tropical mean [ppb]'
    tracertracer_ylabel      = 'ST89 tropical mean [ppb]'
    global_evol_title        = 'Global tropical tracer evolution'
    e90_eruption_fit_label   = 'E90 tropical strat. mean (+- stderr)'
    st80_eruption_fit_label  = 'ST80 tropical trop. mean (+- stderr)'
    if(DO_ANNUAL):
        e90_label  = 'E90 tropical strat. mean [ppb]'
        st80_label = 'ST80 tropical trop. mean [ppb]'
    else:
        e90_label  = 'E90 tropical strat. mean anomaly [ppb]'
        st80_label = 'ST80 tropical trop. mean anomaly [ppb]'
else:
    tracertracer_strat_title = 'STRATOSPHERE'
    tracertracer_trop_title  = 'TROPOSPHERE'
    tracertracer_xlabel      = 'E90 global mean [ppb]'
    tracertracer_ylabel      = 'ST89 global mean [ppb]'
    global_evol_title        = 'Global tracer evolution'
    e90_eruption_fit_label   = 'E90 global strat. mean (+- stderr)'
    st80_eruption_fit_label  = 'ST80 global trop. mean (+- stderr)'
    if(DO_ANNUAL):
        e90_label  = 'E90 global strat. mean [ppb]'
        st80_label = 'ST80 global trop. mean [ppb]'
    else:
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

# scale scatter marker areas with eruption magnitude
s_erupt_areas = (np.array(list(eruption_mags.values())))**2.5*np.pi
# mask for years with/without an eruption
emask = s_erupt_areas > np.pi

# e90-st80 stratopshere tracer-tracer correlation. Non-eruption years marked with dots
ax1.scatter(e90_strat_mean, st80_strat_mean, marker='o', s=s_area**2*np.pi, 
            c=colors, cmap=cmap)
# dummy plot for legend
ax1.scatter([], [], marker='o', color=cmap(0.8), linestyle='None', 
            s = s_area**2*np.pi, label=scatter_label)

if(PLOT_ERUPTIONS):
    # Eruption years marked with stars, scaled for size
    ax1.scatter(e90_strat_mean[emask], st80_strat_mean[emask], marker='*', 
                edgecolors='lime', s=s_erupt_areas[emask], c=colors[emask], 
                cmap=cmap, linewidth=mew, vmin=vmin, vmax=vmax)
    # dummy plot for legend
    ax1.scatter([], [], marker='*', color=cmap(0.8), linestyle='None', 
               s = s_area**6*np.pi, label='eruptions (size = TgS)')

ax1.set_title(tracertracer_strat_title)
ax1.set_xlabel(tracertracer_xlabel)
ax1.set_ylabel(tracertracer_ylabel)
legend = ax1.legend(loc='upper left', fancybox=False)


# colorbar axis. This won't be used, and is just a padding to make sure both the strat 
# and trop correlation axes are the same size
divider = make_axes_locatable(ax1)
cax1 = divider.append_axes("right", size="5%", pad=0.05)
cax1.axis('off')

# e90-st80 tropopshere tracer-tracer correlation, point styling same as above
im = ax2.scatter(e90_trop_mean, st80_trop_mean, marker='o', s=s_area**2*np.pi, 
            c=colors, cmap=cmap)

if(PLOT_ERUPTIONS):
    ax2.scatter(e90_trop_mean[emask], st80_trop_mean[emask], marker='*', edgecolors='lime',
                s=s_erupt_areas[emask], c=colors[emask], cmap=cmap, linewidth=mew, 
                vmin=vmin, vmax=vmax)
    
ax2.set_title(tracertracer_trop_title)
ax2.set_xlabel(tracertracer_xlabel)
# colorbar axis
divider = make_axes_locatable(ax2)
cax2 = divider.append_axes("right", size="5%", pad=0.05)
# colorbar
cb = plt.colorbar(im, cax=cax2)
cb.set_label(cmap_label)


# -------------------------------------------------------
# --------------- global tracer evolution ---------------

# E90 strat mean on the left axis (black), ST80 trop mean on the right axis (grey)
#ax33 = ax3.twinx()
ax3.set_zorder(99)
ax3.set_frame_on(False)

ax3.plot(time, e90_strat_line_data, linestyle, color=e90_color, lw=0.75)
#ax33.plot(time, st80_trop_line_data, linestyle, color=st80_color, lw=0.75, zorder=0)

# plot best fit lines and confidence intervals for the annual analysis
# we already manually performed the best fit to obtain r- and p-values
# we also want to plot confidence intervals... rather than doing the bootstrapping
# ourselves, we use seaborn's "regpot" function, and assume that the fit 
# is consistent with that found above
if(PLOT_TRENDS and DO_ANNUAL):
    #sns.regplot(data=st80_trop_line_data.to_dataframe(), color=st80_color, 
    #            x='year', y='ST80_25j', ci=REG_CI, n_boot=REG_NBOOT, scatter=False, ax=ax33)
    sns.regplot(data=e90_strat_line_data.to_dataframe(), color=e90_color,
                x='year', y='E90j', ci=REG_CI, n_boot=REG_NBOOT, scatter=False, ax=ax3, 
                label='E90 best fit w/ 95%ci (r,p) = ({:.2f}, {:.2e})'\
                       .format(e90_strat_fit.rvalue, e90_strat_fit.pvalue))
    # dummy plot for st80 legend
    ax3.plot([],[], '-', color=st80_color, label='ST80 best fit w/ 95%ci (r,p) = ({:.2f}, {:.2e})'\
                                .format(st80_trop_fit.rvalue, st80_trop_fit.pvalue))

if(PLOT_ERUPTIONS):
    # draw vertical lines on eruption years, with the thickness corresponding to the 
    # injection magnitude
    for year in eruption_mags.keys():
        if(eruption_mags[year] > 0):
            ax3.axvline(year, color=eruption_color, lw=(eruption_mags_normed[year])*2.5, zorder=0)
    ypos = np.mean(e90_strat_line_data)
    # dummy plot for legend
    ax3.plot([1980,1980], [ypos,ypos], color=eruption_color, label='eruption (line width = TgS)')

ax3.legend(fancybox=False, loc='lower left')

if(not DO_ANNUAL):
    # these were eyeballed and work well for this historical data
    ax3.set_xlim([1960, 2000])
    ax3.set_ylim([-0.4, 0.4])
    #ax33.set_ylim([-0.3, 0.3])
else:
    # buffer lower ylim by 20% to make room for legend
    ax3ylim = ax3.get_ylim()
    #ax33ylim = ax33.get_ylim()
    ax3.set_ylim([ax3ylim[0] - np.diff(ax3ylim)*0.2, ax3ylim[1]])
    #ax33.set_ylim([ax33ylim[0] - np.diff(ax33ylim)*0.2, ax33ylim[1]])
     # make st80 axis width match e90
    e90axwidth = np.diff(ax3.get_ylim())
    #st80axcenter = np.sum(ax33.get_ylim())/2
    #ax33.set_ylim([st80axcenter - e90axwidth/2, 
    #               st80axcenter + e90axwidth/2,])
    
ax3.set_title(global_evol_title)
ax3.set_xlabel('time')
ax3.set_ylabel(e90_label)
#ax33.set_ylabel(st80_label)
#ax33.yaxis.label.set_color(st80_color)
ax3.yaxis.label.set_color(e90_color)
    
    
    
# ------------------------------------------------------------
# --------------- e90 tropical 90ppb evolution ---------------

# if we're not analyzing the annual data, then we will plot anomalies 
# here, and show them over a righter time window (near Pinatubo)
if(not DO_ANNUAL):
    
    # We show both the e90 90ppb position anomaly (blue), and the tropopause anoamly (grey)
    ax4.plot(time, e90_90ppb_eq_lev_anom, linestyle, color=e90_color, lw=0.75, 
             label='e90 tropical mean (+- 20 deg) 90 ppb position anomaly')
    ax4.plot(time, trop_t_eq_mean_anom, linestyle, color=trop_color, lw=0.75, 
            label='tropical mean (+- 20 deg) tropopause position anomaly')
    ax4.set_title('Tropopause evolution')
    ax4.set_xlabel('time [years]')
    ax4.set_ylabel('lev [hPa]')
    ax4.legend(fancybox=False, loc='lower left')
    ax4.set_ylim([-11, 11])
else:

    # if we are analyzing annual data, don't show anomalies, just plot
    # the e90 90ppb contour (blue) and tropopause position (grey)
    cf  = ax4.plot(time, e90_90ppb_eq_lev, linestyle, color='cornflowerblue',
                   label='e90 tropical mean (+- 20 deg) 90 ppb position', lw=0.75)

    ax4.plot(time, trop_t_eq_mean, linestyle, color=trop_color, lw=0.75, 
             label='tropical mean (+- 20 deg) tropopause position')
    ax4.set_ylim([94, 116])
    
    # plot best fit lines and confidence intervals for the annual analysis
    # we already manually performed the best fit to obtain r- and p-values
    # we also want to plot confidence intervals... rather than doing the bootstrapping
    # ourselves, we use seaborn's "regpot" function, and assume that the fit 
    # is consistent with that found above
    if(PLOT_TRENDS):
        sns.regplot(data=e90_90ppb_eq_lev.to_dataframe(), color=e90_color,
                    x='year', y='E90j_90ppb_lev', ci=REG_CI, n_boot=REG_NBOOT, scatter=False, ax=ax4, 
                    label='best fit w/ 95%ci, (r,p) = ({:.2f}, {:.2e})'\
                           .format(e90_90ppb_eq_fit.rvalue, e90_90ppb_eq_fit.pvalue))
        trop_t_eq_mean_df = pd.DataFrame({'year': time, 'trop': trop_t_eq_mean})
        sns.regplot(data=trop_t_eq_mean_df, color=trop_color,
                    x='year', y='trop', ci=REG_CI, n_boot=REG_NBOOT, scatter=False, ax=ax4,
                    label='best fit w/ 95%ci, (r,p) = ({:.2f}, {:.2e})'\
                          .format(trop_eq_fit.rvalue, trop_eq_fit.pvalue))
        
    ax4.set_title('Tropical tropopause evolution')
    ax4.set_xlabel('time [years]')
    ax4.set_ylabel('lev [hPa]')

if(PLOT_ERUPTIONS):
    # draw vertical lines on eruption years, with the thickness corresponding to the 
    # injection magnitude
    for year in eruption_mags_eq.keys():
        if(eruption_mags_eq[year] > 0):
            ax4.axvline(year, color=eruption_color, zorder=0,
                        lw=(eruption_mags_eq_normed[year])*2.5)
    # dummy plot for legend
    ax4.plot([1980,1980], [0,0], color=eruption_color,
             label='tropical eruptions (+-line width = TgS)')

# if we aren't analyzing the annual data, then zoom these fgures into the 
# volcanically active period contining Pinatubio and El Chichon
if(not DO_ANNUAL):
    ax4.set_xlim([1960, 2000])
    
# done; draw legend
ax4.legend(fancybox=False, loc='lower left', ncol=2)
ax4.invert_yaxis()


# -------------------------------------------
# --------------- save figure ---------------
plt.tight_layout()
plt.savefig('figs/SCATTER_{}_trend{}_tropical{}_erupt{}_tb{}_ml{}_wlev{}_wlat{}.png'.format(['monthly', 'annual'][int(DO_ANNUAL)], int(PLOT_TRENDS), int(TROPICAL_MEANS),int(PLOT_ERUPTIONS), TROP_BOUND, MIN_LEV, int(WEIGHT_LEV), int(WEIGHT_LAT)))


if(DO_ANNUAL): exit()

# ---------------------------------------------------------------------------
# --------------- second figure: eruption mags vs. fit slopes ---------------

print('---- plotting post-eruption trends...')
# --------------- create figure ---------------
# figure will have two plots; eruption mag. vs e90_strat and st80_trop slopes,
# and tropical eruption mag vs. e90_90ppb line and tropopause position

fig = plt.figure(figsize=(13,5))
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)
ms=10
fit_mask    = [m > REG_MIN_TGS for m in eruption_fit_mags]
fit_eq_mask = [m > REG_MIN_TGS for m in eruption_eq_fit_mags]
if(not PLOT_TRENDS): scatter_alpha = 1
if(PLOT_TRENDS):     scatter_alpha = 0.5

# whether or not to include ST80 in the mean ppb tendecy plot...
PLOT_ST80 = False

# ---- eruption mag. vs e90_strat and st80_trop slopes
# scale slopes to ppb/month from ppb/year
e90_strat_anom_slopes = np.array([e.slope for e in e90_strat_anom_eruption_fits]) / 12
e90_strat_anom_stderr = np.array([e.stderr for e in e90_strat_anom_eruption_fits]) / 12
st80_trop_anom_slopes = np.array([e.slope for e in st80_trop_anom_eruption_fits]) / 12
st80_trop_anom_stderr = np.array([e.stderr for e in st80_trop_anom_eruption_fits]) / 12
ax1.errorbar(eruption_fit_mags, e90_strat_anom_slopes, yerr=e90_strat_anom_stderr, 
             fmt='^', ms=ms, color=e90_color, alpha=scatter_alpha, 
             label=e90_eruption_fit_label)
if(PLOT_ST80):
    ax1.errorbar(eruption_fit_mags, st80_trop_anom_slopes, yerr=st80_trop_anom_stderr, 
                 fmt='^', ms=ms, color=st80_color, alpha=scatter_alpha, 
                 label=st80_eruption_fit_label)

if(PLOT_TRENDS):
    # add linear regression to this correlation plot, with confidence interval estimated
    # via bootstrapping, and data points weighted inversely by the slope stderr. Only
    # consider eruptions greater than REG_MIN_TGS
    ds = pd.DataFrame({'slope':e90_strat_anom_slopes[fit_mask], 'mag':eruption_fit_mags[fit_mask]}, 
                      columns=['slope', 'mag'])
    sns.regplot(data=ds, color=e90_color, x='mag', y='slope', ci=95, scatter=False, ax=ax1)
    if(PLOT_ST80):
        ds = pd.DataFrame({'slope':st80_trop_anom_slopes[fit_mask], 'mag':eruption_fit_mags[fit_mask]}, 
                          columns=['slope', 'mag'])
        sns.regplot(data=ds, color=st80_color, x='mag', y='slope', ci=95, scatter=False, ax=ax1)
    # dummy plots for legend
    ax1.plot([],[],'-',color=e90_color,label='E90 best fit w/ 95% ci')
    if(PLOT_ST80): ax1.plot([],[],'-',color=st80_color,label='ST80 best fit w/ 95% ci')
    
ax1.set_xlabel('eruption magnitude [TgS]')
ax1.set_ylabel('{}-month mean tracer tendency [ppb/month]'.format(ERM))
ax1.legend(loc='lower right')

# mean concentrations are higher wehen restricting view to tropics
if(not TROPICAL_MEANS):
    ax1.set_ylim([-0.12, 0.12])
else:
    ax1.set_ylim([-0.13, 0.15])
ax1.grid()

# ---- eruption mag. vs e90_90ppb and tropopause height slopes
# scale slopes to hPa/month from hPa/year
e90_90ppb_anom_slopes = np.array([e.slope for e in e90_90ppb_eq_anom_eruption_fits]) / 12
e90_90ppb_anom_stderr = np.array([e.stderr for e in e90_90ppb_eq_anom_eruption_fits]) / 12
trop_anom_slopes      = np.array([e.slope for e in trop_eq_anom_eruption_fits]) / 12
trop_anom_stderr      = np.array([e.stderr for e in trop_eq_anom_eruption_fits]) / 12
ax2.errorbar(eruption_eq_fit_mags, e90_90ppb_anom_slopes, yerr=e90_90ppb_anom_stderr, 
             fmt='^', color=e90_color, ms=ms, alpha=scatter_alpha,
             label='E90 tropical-mean 90ppb position anomaly (+- stderr)')
ax2.errorbar(eruption_eq_fit_mags, trop_anom_slopes, yerr=trop_anom_stderr, 
             fmt='^', color=trop_color, ms=ms, alpha=scatter_alpha,
             label='tropical-mean tropopause position anomaly (+- stderr)')

if(PLOT_TRENDS):
    # add linear regression to this correlation plot, with confidence interval estimated
    # via bootstrapping, and data points weighted inversely by the slope stderr. Only
    # consider eruptions greater than REG_MIN_TGS
    ds = pd.DataFrame({'slope':e90_90ppb_anom_slopes[fit_eq_mask], 
                       'mag':eruption_eq_fit_mags[fit_eq_mask]}, 
                       columns=['slope', 'mag'])
    sns.regplot(data=ds, color=e90_color, x='mag', y='slope', ci=95, scatter=False, ax=ax2)
    ds = pd.DataFrame({'slope':trop_anom_slopes[fit_eq_mask], 
                       'mag':eruption_eq_fit_mags[fit_eq_mask]}, 
                       columns=['slope', 'mag'])
    sns.regplot(data=ds, color=trop_color, x='mag', y='slope', ci=95, scatter=False, ax=ax2)
    # dummy plots for legend
    ax1.plot([],[],'-',color=e90_color,label='E90 90ppb best fit w/ 95% ci')
    ax1.plot([],[],'-',color=trop_color,label='tropopause best fit w/ 95% ci')
    
ax2.set_xlabel('eruption magnitude [TgS]')
ax2.set_ylabel('{}-month mean altitude tendency [hPa/month]'.format(ERM))
ax2.legend(loc='upper right')
ax2.set_ylim([-3, 3])
ax2.grid()
ax2.invert_yaxis()

# -------------------------------------------
# --------------- save figure ---------------
plt.tight_layout()
plt.savefig('figs/TRENDS_{}_tropical{}_erm{:02d}_tb{}_ml{}_wlev{}_wlat{}.png'.format(
            ['monthly', 'annual'][int(DO_ANNUAL)], int(TROPICAL_MEANS), 
            ERM, TROP_BOUND, MIN_LEV, int(WEIGHT_LEV), int(WEIGHT_LAT)))




