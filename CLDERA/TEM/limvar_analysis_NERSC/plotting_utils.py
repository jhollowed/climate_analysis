import pdb
import sys
import copy
import scipy
import warnings
import numpy as np
import xarray as xr
import time as timer
import matplotlib as mpl
from datetime import datetime
from matplotlib import colors
from datetime import timedelta
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.dates as mdates
import matplotlib.ticker as ticker
import matplotlib.patches as patches
from collections.abc import Sequence
from scipy.interpolate import griddata
from scipy.interpolate import interp1d
from contourpy import contour_generator
from matplotlib.dates import DateFormatter
from matplotlib.ticker import FuncFormatter
from matplotlib.ticker import ScalarFormatter
from matplotlib.colors import LinearSegmentedColormap

# local imports
import compute_ensemble_stats as ces

sys.path.insert(1, '/global/homes/j/jhollo/repos/climate_analysis/CLDERA/wavePaperFigs/util')
import nclcmaps as ncm

# --- constants
g = 9.80665    # global average of gravity at mean sea level in m/s^2
H = 7000       # scale height
a = 6.37123e6  # radius of Earth in m
spd = 24*60*60 # seconds per day


# =======================================================================
# =======================================================================


def get_variable(var, mass, freq='10daily', q=None, 
                 pmin=None, pmax=None, latmin=None, latmax=None, tmin=None, tmax=None, 
                 average_lat=True, average_pres=True, average_time=True, overwrite=False, 
                 debug=False, skip_nosrctag=False, return_intersection=True, 
                 return_members=False, pass_var=False):
    '''
    returns a dictionary of the ensemble data, TEM data, TEM budget data, as well as those dataset's 
    counterfactuals, impacts, p-values, and coherence for a specfied variable with a specified 
    slicing. This variable can be stored in either the ensemble mean, TEM, or TEM budget data. They 
    will all be searched for a variable with the matching name. 
    This function calls compute_ensemble_stats.get_data_and_stats() in order to get the processed data. 
    If the data has not yet been processed for this choice of slicing, it will be computed. If it has, 
    it will simply be read and returned.
    The slicing uses plotting_utils.do_slicing() (see docs therein), with the same input arguments, 
    but with average=True. Meaning, the data will be averaged over the specified slice before the 
    stats are computed.

    Parameters
    ----------
    var : str
        the variable name
    mass : int
        the SO2 mass of the ensemble, in Tg. Options are 10 or 15
    freq : str, optional
        which frequency dataset to use. Options are:
        '10daily' (default) : the 10-daily averaged data
        'monthly'           : the monthly-averaged data
    q : str, optional
        tracer name, either 'aoa', 'e90', or None (the default), in which case
        the specified variable is not a traer quantity
    pmin : float, optional
        minimum pressure (highest vertical level)
    pmax : float, optional
        maximum pressure (lowest vertical level)
    latmin : float, optional
        minimum latitude 
    latmax : float, optional
        maximum latitude
    tmin : cftime._cftime.Datetime, optional
        minimum time
    tmax : cftime._cftime.Datetime, optional
        maximum time
    average_lat : bool, optional
        whether or not to average over the latitude slice specified by latmin,latmax
        If True, slicing is done BEFORE integrating, and otherwise, is done AFTER
        Defaults to True
    average_pres : bool, optional
        whether or not to average over the vertical slice specified by pmin,pmax
        If True, slicing is done BEFORE integrating, and otherwise, is done AFTER
        Defaults to False
    average_time : bool, optional
        whether or not to average over the temporal slice specified by tmin,tmax
        Defaults to True
    overwrite : bool, optional
        whether or not to overwrite the data produced by a previous run of this function.
        Defaults to False, in which case the processed data is simply read and returned 
        if it already exists.
    debug : bool, optional
        if True, print out available variables for each dataset. Defaults to False.
    skip_nosrctag : bool, optional
        whether or not to exclude the non-source-tagged ensemble members (ens91-95) from 
        the returned data. Defaults to False.
    return_intersection : bool, optional
        whether or not to return only the portion of data where the data and the 
        counterfactual data overlap in time. Defaults to True.
    return_members : bool, optional
        whether or not to return the member-level data.
        The counterfactual, forced runs, impacts and coherence are all provided with an
        'ens' dimension. Defaults to False.
    pass_var : bool, optional
        whether or not to pass the variable name requested to get_data_and_stats().
        Defaults to False, in which case the data is processed for the entire dataset
        under the slicing criteria specified.
        If True, then only the variable rquested in processed
        Setting to True will be much faster if only one variable is needed.
        If many variables are going to be requested, then setting this to False will be
        faster
        
    Returns
    -------
    dict
        A dictionary containing four xarray Datasets with the keys
        'ensmean', 'cfmean', 'impact', 'pval'
    '''
    
    # ---- get tracer index
    qi = {None:0, 'aoa':1, 'e90':2}
    qi = qi[q]
    
    # ---- get processing args
    args={'overwrite':overwrite, 'pmin':pmin, 'pmax':pmax, 'average_pres':average_pres, 
          'latmin':latmin, 'latmax':latmax, 'average_lat':average_lat,
          'tmin':tmin, 'tmax':tmax, 'average_time':average_time,
          'skip_nosrctag':skip_nosrctag, 'return_intersection':return_intersection, 
          'return_members':return_members}
    if(pass_var): args['var'] = var
    
    # ---- ensemble data
    try:
        if(return_members):
            data, cf, impact, pval, coherence, data_members, cf_members, impact_members = \
                            ces.get_data_and_stats(dataset='ens', mass=mass, freq=freq, qi=0, **args)
        else:
            data, cf, impact, pval, coherence = ces.get_data_and_stats(dataset='ens', mass=mass, freq=freq, qi=0, **args)
        if(debug):
            print('available data vars: {}'.format(q, list(data.data_vars)))
            
        data, cf, impact, pval, coherence = [data[var], cf[var], impact[var], pval[var], coherence[var]]
        data_dict = {'ensmean':data, 'cfmean':cf, 'impact':impact, 'pval':pval, 'coherence':coherence}
        if(return_members):
            data_members, cf_members, impact_members = [data_members[var], cf_members[var], impact_members[var]]
            data_dict.update({'members':data_members, 'cf_members':cf_members, 'impact_members':impact_members})
        return data_dict
    except KeyError: pass
    
    # ---- tem data
    try:
        if(return_members):
            tem_data, tem_cf, tem_impact, tem_pval, tem_coherence, tem_members, tem_cf_members, tem_impact_members = \
                   ces.get_data_and_stats(dataset='tem', mass=mass, freq=freq, qi=qi, **args)
        else:
            tem_data, tem_cf, tem_impact, tem_pval, tem_coherence = \
                              ces.get_data_and_stats(dataset='tem', mass=mass, freq=freq, qi=qi, **args)
        if(debug):
            print('available TEM vars for tracer {}: {}'.format(q, list(tem_data.data_vars)))

        data, cf, impact, pval, coherence = [tem_data[var], tem_cf[var], tem_impact[var],
                                             tem_pval[var], tem_coherence[var]]
        data_dict = {'ensmean':data, 'cfmean':cf, 'impact':impact, 'pval':pval, 'coherence':coherence}
        if(return_members):
            data_members, cf_members, impact_members = [tem_members[var], tem_cf_members[var], 
                                                        tem_impact_members[var]]
            data_dict.update({'members':data_members, 'cf_members':cf_members, 'impact_members':impact_members})
        return data_dict
    except KeyError: pass
    
    # ---- budget data
    try:
        if(return_members):
            budget_data, budget_cf, budget_impact, budget_pval, budget_coherence,\
            budget_members, budget_cf_members, budget_impact_members = \
                                    ces.get_data_and_stats(dataset='budget', mass=mass, freq=freq, qi=qi, **args)
        else:
            budget_data, budget_cf, budget_impact, budget_pval, budget_coherence = \
                                    ces.get_data_and_stats(dataset='budget', mass=mass, freq=freq, qi=qi, **args)
        if(debug):
            print('available TEM budget vars for tracer{}: {}'.format(q, list(budget_data.data_vars)))

        data, cf, impact, pval, coherence = [budget_data[var], budget_cf[var], budget_impact[var],
                                             budget_pval[var], budget_coherence[var]]
        data_dict = {'ensmean':data, 'cfmean':cf, 'impact':impact, 'pval':pval, 'coherence':coherence}
        if(return_members):
            data_members, cf_members, impact_members = [budget_members[var], budget_cf_members[var], 
                                                        budget_impact_members[var]]
            data_dict.update({'members':data_members, 'cf_members':cf_members, 'impact_members':impact_members})
        return data_dict
    except KeyError: pass

    # ---- if variable hasn't been found by here, it doesn't exist
    assert isinstance(data, xr.core.dataarray.DataArray), 'variable {} not found!'.format(var)

# -----------------------------------------------------------------------

def do_slicing(data, pmin=None, pmax=None, latmin=None, latmax=None, 
               tmin=None, tmax=None, average=False):
    '''
    Slices the data optionally on pressure, latitude, and time.

    Parameters
    ----------
    data : xarray DataArray, or list of DataArrays, or dict of DataArrays, or Dataset
        The data to slice. If list or dict, each item must be an xarray DataArray, and 
        the slicing is applied to each
    pmin : float, optional
        minimum pressure (highest vertical level)
    pmax : float, optional
        maximum pressure (lowest vertical level)
    latmin : float, optional
        minimum latitude 
    latmax : float, optional
        maximum latitude
    tmin : cftime._cftime.Datetime, optional
        minimum time
    tmax : cftime._cftime.Datetime, optional
        maximum time
    average : bool, optional
        whether or not to average the data over the specified slice. Defaults to False
    '''

    # --- check inputs
    assert (np.array([pmin, pmax]) == None).sum() != 1,     "must pass both or neither or pmin, pmax"
    assert (np.array([latmin, latmax]) == None).sum() != 1, "must pass both or neither or latmin, latmax"
    assert (np.array([tmin, tmax]) == None).sum() != 1,     "must pass both or neither or tmin, tmax"
    
    assert(type(data) == type(xr.DataArray()) or type(data) == type(xr.Dataset()) or \
           type(data) == list or type(data) == dict),\
           "data must be a DataArray, Dataset, or a list or dict of those"
    if(type(data) == dict):
        data_keys = data.keys()
        if(('impact' in data_keys or 'pval' in data_keys or 'coherence' in data_keys) and average==True):
            raise RuntimeError('It is incorrect to take an average of impacts, p-values, or coherence! '\
                               'Either set average=False, or pass the args of this call to do_slicing() '\
                               'to get_variable() isntead, in which case the averaging will be done '\
                               'BEFORE p-values are computed.')
    else:
        data_keys = None
    
    if(type(data) == list): 
        assert np.sum([type(d) in [type(xr.DataArray()), type(xr.Dataset())] for d in data]) == len(data),\
               "if data is a list, each element must be a DataArray or Dataset"
    elif(type(data) == dict):
        assert np.sum([type(data[k]) in [type(xr.DataArray()), type(xr.Dataset())] for k in data_keys]) == len(data),\
               "if data is a dict, each element must be a DataArray or Dataset"
        data = list(data.values())
    else:
        data = [data]

    # --- do slicing and averaging
    for i in range(len(data)):
       
        # pressure
        if(pmin is not None):
            data[i] = data[i].sel(plev = slice(pmin, pmax))
            if(average):
                # weight by the pressure thickness of each level
                weights = xr.DataArray(np.gradient(data[i].plev), dims=('plev'))
                weights.name = "weights"
                data[i] = data[i].weighted(weights)
                data[i] = data[i].mean(("plev"))
        # latitude
        if(latmin is not None): 
            data[i] = data[i].sel(lat = slice(latmin, latmax))
            if(average):
                # weight by cos(lat)
                weights = np.cos(np.deg2rad(data[i].lat))
                weights.name = "weights"
                data[i] = data[i].weighted(weights)
                data[i] = data[i].mean(("lat"))
        # time
        if(tmin is not None):   
            data[i] = data[i].sel(time = slice(tmin, tmax))
            if(average):
                data[i] = data[i].mean('time')
    
    # --- return
    # if a list was passed, return a list. Otherwise return a DataArray or Dataset
    if(len(data) == 1): 
        return data[0]
    else: 
        if(data_keys is None):
            return data
        else:
            return dict(zip(data_keys, data))
        
# -----------------------------------------------------------------------
    
def get_integrated_tendency(tend, tinit, mass, freq='10daily', q=None,
                            pmin=None, pmax=None, average_pres=False, 
                            latmin=None, latmax=None,average_lat=True,
                            tmin=None, tmax=None, average_time=False, 
                            overwrite=False, skip_nosrctag=False):
    '''
    Returns the ensemble mean, impact, coherence, and pvalues for integrated TEM tendencies
    for either the daily or 10-daily data, with specified data slicing.
    The slicing uses plotting_utils.do_slicing() (see docs therein), with the same input 
    arguments, but with average=True. Meaning, the data will be averaged over the specified 
    slice before the stats are computed. If wanting to take slices without averaging, use 
    plotting_utils.dos_slicing() after loading the data.

    Parameters
    ----------
    tinit : cftime.datetime object
        an initial time from which to integrate the wind and tracer time-tendencies 
        present in the data. If this time is not an exact match to a time coordinate 
        in the data, the initial condition used will be the next nearest time in the 
        dataset.
    mass : int
        the SO2 mass of the ensemble, in Tg. Options are 10 or 15
    freq : str, optional
        which frequency dataset to use. Options are:
        'daily' (default)   : the daily data
        '10daily'           : the 10-daily averaged data
    q : str, optional
        tracer name, either 'aoa', 'e90', or None (the default), in which case
        the specified variable is not a traer quantity
    pmin : float, optional
        minimum pressure (highest vertical level)
    pmax : float, optional
        maximum pressure (lowest vertical level)
    latmin : float, optional
        minimum latitude 
    latmax : float, optional
        maximum latitude
    tmin : cftime._cftime.Datetime, optional
        minimum time
    tmax : cftime._cftime.Datetime, optional
        maximum time
    average_lat : bool, optional
        whether or not to average over the latitude slice specified by latmin,latmax
        If True, slicing is done BEFORE integrating, and otherwise, is done AFTER.
        Defaults to True
    average_space : bool, optional
        whether or not to average over the vertical slice specified by pmin,pmax
        If True, slicing is done BEFORE integrating, and otherwise, is done AFTER.
        Defaults to True
    average_time : bool, optional
        whether or not to average over the temporal slice specified by tmin,tmax
        If True, slicing is done BEFORE integrating, and otherwise, is done AFTER.
        Defaults to True
    overwrite : bool, optional
        whether or not to overwrite the data produced by a previous run of this function.
        Defaults to False, in which case the processed data is simply read and returned 
        if it already exists.
    '''
    # ---- get tracer index
    qi = {None:0, 'aoa':1, 'e90':2}
    qi = qi[q]
    
    # ---- get processing args
    args={'pmin':pmin, 'pmax':pmax, 'latmin':latmin, 'latmax':latmax, 
          'tmin':tmin, 'tmax':tmax, 'average_lat':average_lat, 
          'average_pres':average_pres, 'average_time':average_time,
          'overwrite':overwrite, 'skip_nosrctag':skip_nosrctag}
    
    # ---- get integrated tendencies
    data, cf, impact, pval, coherence = ces.get_integrated_tendencies_and_stats(tinit, mass, freq, qi, **args)
    
    # ---- get requested tendency variable and return
    data, cf, impact, pval, coherence = [data[tend], cf[tend], impact[tend], pval[tend], coherence[tend]]
    return {'ensmean':data, 'cfmean':cf, 'impact':impact, 'pval':pval, 'coherence':coherence}

# -----------------------------------------------------------------------

def sel(data, plev=None, lat=None, time=None, method='nearest'):
    '''
    Selects the data optionally on pressure, latitude, and time.

    Parameters
    ----------
    data : xarray DataArray, or list of DataArrays, or dict of DataArrays
        The data to select. If list or dict, each item must be an xarray DataArray, and 
        the slicing is applied to each
    plev : float, optional
        pressure in hPa
    lat : float, optional
        latitude in degrees
    time : cftime._cftime.Datetime, optional
        time as a date
    method : str, optional
        selection method for xarray.DataArray.sel(). Defaults to 'nearest'
    '''

    # --- check inputs
    assert(type(data) == type(xr.DataArray()) or type(data) == type(xr.Dataset()) or \
           type(data) == list or type(data) == dict),\
           "data must be a DataArray, Dataset, or a list or dict of those"
    if(type(data) == dict):
        data_keys = data.keys()
    else:
        data_keys = None
    
    if(type(data) == list): 
        assert np.sum([type(d) in [type(xr.DataArray()), type(xr.Dataset())] for d in data]) == len(data),\
               "if data is a list, each element must be a DataArray or Dataset"
    elif(type(data) == dict):
        assert np.sum([type(data[k]) in [type(xr.DataArray()), type(xr.Dataset())] for k in data_keys]) == len(data),\
               "if data is a dict, each element must be a DataArray or Dataset"
        data = list(data.values())
    else:
        data = [data]

    # --- do selection
    for i in range(len(data)):
        if(plev is not None):
            data[i] = data[i].sel(plev=plev, method=method)
        if(lat is not None):
            data[i] = data[i].sel(lat=lat, method=method)
        if(time is not None):
            data[i] = data[i].sel(time=time, method=method)
    
    # --- return
    # if a list was passed, return a list. Otherwise return a DataArray or Dataset
    if(len(data) == 1): 
        return data[0]
    else: 
        if(data_keys is None):
            return data
        else:
            return dict(zip(data_keys, data))

# -----------------------------------------------------------------------

def filter_significance(data, sigvar, sigtype='pval', thresh=None, fill=None):
    '''
    Filters an input dataset by signifiance. Any datapoints that are
    not significant are replaced with nans. The signifiance definition
    depends on the arguments sigvar and thresh, which set a treshold
    level in either p-value or coherence.
    
    Parameters
    ----------
    data : xarray Dataset or DataArray
        The data to filter
    sigvar : xarray DataArray
        variable to use for defining significance
    sigtype : str, optional
        specifies what type of data is in sigvar
        Options are 'pvalue' or 'coherence'. Default is 'pval'
    thresh : float, optional
        threshold value for defining significance in the variable sigvar.
        Default for sigtype='pvalue' is 0.05
        Default for sigtype='coherence' is 0.75
    fill : float
        fill value for deleted values, rather than nan
        
    Returns
    -------
    xarray DataArray or Dataset
        The input data dictionary, with all insignificant values replaced
        with nans
    '''
    
    assert sigtype in ['pval', 'coherence'], 'sigtype must be either \'pval\' or \'coherence\''
    
    if(sigtype == 'pval'):
        if(thresh is None): thresh = 0.05
        sig = sigvar < thresh
    if(sigtype == 'coherence'):
        if(thresh is None): thresh = 0.75
        sig = sigvar > thresh
    
    if(fill is None):
        return data.where(sig)
    else:
        return data.where(sig, other=fill)
    
# -----------------------------------------------------------------------

def filter_vector_significance(v1, v2, sigvar1, sigvar2, sigtype='pval', thresh=None, fill=None):
    '''
    Filters a two-component vector dataset by signifiance. Any datapoints that are
    not significant in at least one vector component are replaced with nans. The 
    signifiance definition depends on the arguments sigvar and thresh, which set a 
    treshold level in either p-value or coherence.
    
    Parameters
    ----------
    v1 : xarray Dataset or DataArray
        The first vector component
    v2 : xarray Dataset or DataArray
        The second vector component
    sigvar : xarray DataArray
        variable to use for defining significance of v1
    sigvar : xarray DataArray
        variable to use for defining significance of v2
    sigtype : str, optional
        specifies what type of data is in sigvar1, sigvar2
        Options are 'pvalue' or 'coherence'. Default is 'pvalue'
    thresh : float, optional
        threshold value for defining significance in the variable sigvar.
        Default for sigtype='pvalue' is 0.05
        Default for sigtype='coherence' is 0.75
    fill : float
        fill value for deleted values, rather than nan
        
    Returns
    -------
    xarray DataArray or Dataset
        The input data dictionary, with all insignificant values replaced
        with nans
    '''
    
    assert sigtype in ['pval', 'coherence'], 'sigtype must be either \'pval\' or \'coherence\''
    
    if(sigtype == 'pval'):
        if(thresh is None): thresh = 0.05
        sig1 = sigvar1 < thresh
        sig2 = sigvar2 < thresh
    if(sigtype == 'coherence'):
        if(thresh is None): thresh = 0.75
        sig1 = sigvar1 > thresh
        sig2 = sigvar2 > thresh
    
    sig = np.logical_or(sig1, sig2)
    
    if(fill is None):    
        return [v1.where(sig), v2.where(sig)]
    else:
        return [v1.where(sig, other=fill), v2.where(sig, other=fill)]

# -----------------------------------------------------------------------

def adjust_10daily_integrated_tendency(x, dxdti):
    '''
    Adjusted the 10-daily integrated tendencies such that the initial condition
    differences that were induced by the 10-daily time averaging process are 
    removed from the cumulative sums
    
    Parameters
    ----------
    x : xarray DataArray
        data to sample initial condition from
    dxdti : xarray DataArray
        integrated time tendency of the data x
    '''
    return dxdti - (dxdti - x).isel(time=0).values

# -----------------------------------------------------------------------

def shift_integrated_tendency(x, dxdti, icdate, enddate=None, ndays=None, do_adjustment=True):
    '''
    Shifts an integrated tendencies in time, given an initial condition date
    
    Parameters
    ----------
    x : xarray DataArray
        data to sample initial condition from
    dxdti : xarray DataArray
        integrated time tendency of the data x
    icdate : cftime._cftime.Datetime
        date to use as initial condition
    enddate : cftime._cftime.Datetime, optional
        date at which to truncate the output data
        If passed, ndays must not be passed
    ndays : int, optional
        number of days following the initial date at which to truncate
        the output data.
        If passed, ndays must not be passed
    do_adustment : bool
        whether or not to adjust the data to align at the initial condition.
        See adjust_10daily_integrated_tendency
    '''
    
    # check inputs
    assert(int(ndays is not None) + int(enddate is not None) < 2), "may pass ndays or enddate, not both"
    assert(dxdti.time.min() <= icdate < dxdti.time.max()), "icdate must be within the data time range"
    if(enddate is not None):
        assert(dxdti.time.min() < enddate <= dxdti.time.max()), "enddate must be within the data time range"
        
    # adjust time-averaged 10-daily initial conditions if requested...
    if(do_adjustment):
        dxdti = adjust_10daily_integrated_tendency(x, dxdti)
    
    # "backfill" is used here so that the selection has the same effect as slice()
    targs = {'time':icdate, 'method':'backfill'}
    x0 = x.sel(**targs)
    x_integrated = x0 + dxdti - dxdti.sel(**targs)

    # return only the data after the initial condition date, and truncate if requested    
    if(ndays is not None):
        enddate = to_datetime(icdate) + timedelta(days=ndays)
    x_integrated = x_integrated.sel(time=slice(icdate, enddate))
    return x_integrated

# -----------------------------------------------------------------------

def downsample_latp(data, dslat=None, dsp=None, logp=True, interp_lat=False, interp_plev=False):
    '''
    Downsamples a variable in a latitude-pressure plane

    Parameters
    ----------
    data : xarray DataArray
        The data to downsample.
        Must have coordinates 'plev' and 'lat'
    dslat : float, optional
        latitude downsampling factor. If None (the default), no downsampling 
        is done. Otherwise, the data will be downsampled in latitude by dslat
        times, starting at the dslat element. That is, if dslat=4, then the
        data is subsampled via slice(4, None, 4)
    dsp : float, optional
        pressure downsampling factor. If None (the default), no downsampling 
        is done. Otherwise, do the following:
        For a linear pressure axis,the data will be downsampled in pressure 
        by dsp times, starting at the dslat element. That is, if dsp=4, then 
        the data is subsampled via slice(4, None, 4)
        For a logarithmic pressure axis, the data will be downsampled to a
        total of nlev/dsp points, with indices logarithmically spaced from 0
        to nlev (where nlev is the number of vertical levels in the data)
    logp : bool, optional
        whether or not to downsample the data logarithmically in the pressure
        axis, if interp_plev=False. Defaults to True, in which case more data
        is sampled from smaller pressures than larger ones.
    interp_lat : bool, optional
        If dslat is are passed, then this argument controls the 
        downsampling method. If interp=False, then instead of downsampling
        by index (as described in the argument descriptions above for dslat), 
        xr.coarsen() will be called, which interpolates the data to a coarser 
        grid.
    interp_plev : bool, optional
        If dsp is are passed, then this argument controls the 
        downsampling method. If interp=False, then instead of downsampling
        by index (as described in the argument descriptions above for dsp), 
        xr.coarsen() will be called, which interpolates the data to a coarser 
        grid.
        
    Returns
    -------
    xarray DataArray
        The downsampled data
    '''
    
    if(dslat is not None):
        if(interp_lat):
            data = data.fillna(1e-15).coarsen(lat=dslat, boundary='trim').mean()
            data = data.where(data!=1e-15)
        else:
            data = data.isel(lat=slice(dslat, None, dslat))
    if(dsp is not None):
        if(interp_plev):
            data = data.coarsen(plev=dsp, boundary='trim').mean()
            data = data.where(data!=1e-15)
        else:
            if(not logp):
                data = data.isel(plev=slice(dsp, None, dsp))
            elif(logp):
                idx = np.logspace(0, np.log10(len(data.plev)-dsp), num=int(len(data.plev)/dsp), base=10)
                idx = np.unique(np.round(idx).astype(int))
                data = data.isel(plev=idx)
    return data
    
# -----------------------------------------------------------------------    

def scale_EP_flux_vectors(fig, ax, epfy, epfz, P0=101325, dslat=None, dsp=None, 
                          dsplog=None, interp_lat=False, interp_plev=False, 
                          log_vectors=False):
    '''
    scales the EP flux vectors for graphical display, according to Jucker (2021), 
    their Table 2. The scaling of the vectors depends on the aspect ratio of the
    figure, and so this must not change after these scalings are computed, else
    the displayed vectors will be incorrect. The following important assumptions
    are made:
     - it is assumed that the limits of the data coordiantes are identical to
       the limits of the coordinate axes in the plot. That is, any slicing of 
       the data must be done FIRST
     - it is assumed that pressure decreases upward on the y-axis
     - it is assumed that latitude increases to the right on the x-axis 

    Parameters
    ----------
    fig : matplotlib figure object
        the figure on which the plotting axis is located.
    ax : matplotlib axis object
        The axis on which the vectors will be plotted.
    epfy : xarray DataArray
        the northward component of the EP flux in log-pressure cordinate
        (epfy from PyTEMDiags, or Eq. A13 from Gerber+Manzini 2016).
        Must have coordinates 'plev' in hPa and 'lat' in deg
    epfz : xarray DataArray
        the upward component of the EP flux in log-pressure cordinate
        (epfz from PyTEMDiags, or Eq. A14 from Gerber+Manzini 2016).
        Must have coordinates 'plev' in hPa and 'lat' in deg
    dslat : float, optional
        latitude downsampling factor. If None (the default), no downsampling 
        is done. Otherwise, the data will be downsampled in latitude by dslat
        times, starting at the dslat element. That is, if dslat=4, then the
        data is subsampled via slice(4, None, 4)
    dsp : float, optional
        pressure downsampling factor. If None (the default), no downsampling 
        is done. Otherwise, do the following:
        For a linear pressure axis,the data will be downsampled in pressure 
        by dsp times, starting at the dslat element. That is, if dsp=4, then 
        the data is subsampled via slice(4, None, 4)
        For a logarithmic pressure axis, the data will be downsampled to a
        total of nlev/dsp points, with indices logarithmically spaced from 0
        to nlev (where nlev is the number of vertical levels in the data)
    dsplog : bool, optional
        whether or not to pass logp=True to downsample_latp(). Default is None, 
        in which case this will be set to True if the yscale of the axis is
        'log', and False otherwise.
    interp_lat : bool, optional
        If dslat is are passed, then this argument controls the 
        downsampling method. If interp=False, then instead of downsampling
        by index (as described in the argument descriptions above for dslat), 
        xr.coarsen() will be called, which interpolates the data to a coarser 
        grid. Defaults to True
    interp_plev : bool, optional
        If dsp is are passed, then this argument controls the 
        downsampling method. If interp=False, then instead of downsampling
        by index (as described in the argument descriptions above for dsp), 
        xr.coarsen() will be called, which interpolates the data to a coarser 
        grid. Defaults to False
    P0 : float, optional
        reference pressure, in Pa. Defaults to 1000 hPa
    '''

    # ----- get display dimensions
    # in matplotlib, we need get get the figure size, and then
    # get the portion of that figure that is covered by the axes
    fig_width, fig_height    = fig.get_size_inches()
    _, _, ax_wsize, ax_hsize = ax.get_position().bounds
    X, Y = fig_width * ax_wsize, fig_height * ax_hsize
    
    # ----- get pressure scale of axis
    pscale = ax.get_yscale()
    logp = pscale == 'log'
    if(dsplog is not None):
        logp=dsplog
    
    # ----- downsample vector components
    epfy = downsample_latp(epfy, dslat, dsp, logp, interp_lat, interp_plev)
    epfz = downsample_latp(epfz, dslat, dsp, logp, interp_lat, interp_plev)
        
    # ----- get dimensions
    plev       = epfy.plev * 100 # get pressure coordinate in Pa 
    lat        = epfy.lat
    acoslat    = a*np.cos(np.deg2rad(lat))
    p0,   p1   = plev.max(), plev.min()
    lat0, lat1 = lat.min(),  lat.max()

    # ----- transform from log-pressure to pressure coordinate
    Fphi = epfy * P0/plev
    Fp   = epfz * -P0/H

    # ----- get EP flux components, i.e. Jucker Eq.1-2
    fphi = Fphi / acoslat
    fp   = Fp   / acoslat

    # ----- apply Edmon (1980) scaling, i.e. Jucker Eq. 4
    hFphi = 2*np.pi/g * acoslat**2 * (fphi)
    hFp   = 2*np.pi/g * acoslat**2 * (a * fp)

    # ----- set the type of pressure scaling
    pscale = ax.get_yscale()
    if(pscale == 'linear'): pdiff = p1-p0
    elif(pscale == 'log'):  pdiff = plev * np.log(p1/p0)
    
    # ----- do scaling
    Fx = hFphi * (X/Y) / ((lat1-lat0) * np.pi/180)
    Fy = hFp * 1/pdiff
    
    # ----- log-scale vector lengths if requested
    if(log_vectors):
        Fx, Fy = log_quiver(Fx, Fy)
    
    return Fx, Fy

# -----------------------------------------------------------------------

def streamfunction_gradient_normal(fig, ax, f=None, dfdlat=None, dfdp=None, 
                                   dslat=None, dsp=None, dsplog=None, interp_lat=True, 
                                   interp_plev=False, log_vectors=False):
    '''
    Computes and returns gradient-normal vectors for an input scalar field, or the 
    gradient components of a scalar field
    (intended to compute a vector field parallel to streamfunction contours)
    
    Parameters
    ----------
    fig : matplotlib figure object
        the figure on which the plotting axis is located.
    ax : matplotlib axis object
        The axis on which the vectors will be plotted.
    f : xarray DataArray, optional
        The input scalar field (streamfunction). Must have coordinates 
        'plev' and 'lat'. Must be passed if fglat and fgp are not
        passed.
    dfdlat : xarray DataArray, optional
        The latitude gradient of the scalar field (streamfunction). 
        Must have coordinates 'plev' and 'lat'. Must be passed if f
        is not passed.
    dfdp : xarray DataArray, optional
        The latitude gradient of the scalar field (streamfunction). 
        Must have coordinates 'plev' and 'lat'. Must be passed if 
        f is not passed.
    dslat : float, optional
        latitude downsampling factor. If None (the default), no downsampling 
        is done. Otherwise, the data will be downsampled in latitude by dslat
        times, starting at the dslat element. That is, if dslat=4, then the
        data is subsampled via slice(4, None, 4)
    dsp : float, optional
        pressure downsampling factor. If None (the default), no downsampling 
        is done. Otherwise, do the following:
        For a linear pressure axis,the data will be downsampled in pressure 
        by dsp times, starting at the dslat element. That is, if dsp=4, then 
        the data is subsampled via slice(4, None, 4)
        For a logarithmic pressure axis, the data will be downsampled to a
        total of nlev/dsp points, with indices logarithmically spaced from 0
        to nlev (where nlev is the number of vertical levels in the data)
    dsplog : bool, optional
        whether or not to pass logp=True to downsample_latp(). Default is None, 
        in which case this will be set to True if the yscale of the axis is
        'log', and False otherwise.
    interp_lat : bool, optional
        If dslat is are passed, then this argument controls the 
        downsampling method. If interp=False, then instead of downsampling
        by index (as described in the argument descriptions above for dslat), 
        xr.coarsen() will be called, which interpolates the data to a coarser 
        grid. Defaults to True
    interp_plev : bool, optional
        If dsp is are passed, then this argument controls the 
        downsampling method. If interp=False, then instead of downsampling
        by index (as described in the argument descriptions above for dsp), 
        xr.coarsen() will be called, which interpolates the data to a coarser 
        grid. Defaults to False
    '''
    
    # ---- check args
    if(f is None):
        assert dfdlat is not None and dfdp is not None, "both dfdlat and dfdp must be passed if f is not passed"
    if(f is not None):
        assert dfdlat is None and dfdp is None, "neither dfdlat or dfdp should be passed if f is passed"
    if(dfdlat is not None): gradient_input = True
    else:                  gradient_input = False
    
    # compute gradient in data units if not passed
    if(not gradient_input):
        dfdlat, dfdp = xr.zeros_like(f), xr.zeros_like(f)
        dfdlat.values, dfdp.values = np.gradient(f, lat, plev)
    
    # ----- get display dimensions
    # in matplotlib, we need get get the figure size, and then
    # get the portion of that figure that is covered by the axes
    fig_width, fig_height    = fig.get_size_inches()
    _, _, ax_wsize, ax_hsize = ax.get_position().bounds
    dX, dY = fig_width * ax_wsize, fig_height * ax_hsize
    
    # ----- get pressure scale of axis
    pscale = ax.get_yscale()
    logp = pscale == 'log'
    if(dsplog is None):
        dsplog = logp
    
    # ----- get dimensions
    lat,  plev = dfdlat.lat, dfdlat.plev
    p0,   p1   = np.max(plev.values), np.min(plev.values)
    lat0, lat1 = np.min(lat.values),  np.max(lat.values)
    dlat       = lat1-lat0
    dp         = p1-p0
    
    # ----- convert gradient to display units
    dfdx = dfdlat * dlat
    if(not logp):
        dfdy = dfdp * (dX/dY) * dp
    else:
        y = np.log(p0/plev)/np.log(p0/p1)
        dfdy = dfdp * (dX/dY) * np.log(p1/p0) * p0*(p1/p0)**y
    
    # ----- get normal vectors to gradient
    # sign convention is chosen such that the normal vectors are clockwise around local maxima
    dfdxn, dfdyn =  -dfdy, dfdx
    
    # ----- downsample if requested
    dfdxn = downsample_latp(dfdxn, dslat, dsp, dsplog, interp_lat, interp_plev)
    dfdyn = downsample_latp(dfdyn, dslat, dsp, dsplog, interp_lat, interp_plev)
    
    # ----- log-scale vector lengths if requested
    if(log_vectors):
        dfdxn, dfdyn = log_quiver(dfdxn, dfdyn)
    
    return dfdxn, dfdyn

# -----------------------------------------------------------------------

def regrid_vectors_latp(u, v, usig=None, vsig=None, logp=True):
    '''
    Regrids a vector field to a regular grid in latitude-pressure, 
    either distributing the pressure levels uniformyl in p or log(p), 
    matching the input number of points
    
    Parameters
    ----------
    u,v : xarray DataArray, xarray DataArray
        the vector components in latitude and pressure
    usig,vsig : xarray DataArray, xarray DataArray, optional
        either pvalues or coherence values for u and v, 
        which can optionally also be interpolated
    logp : bool
        whether or not to build the interpolation grid
        uniformly in log-pressure
        
    Returns
    -------
    ui, vi : xarray DataArray, xarray DataArray
        the interpolated vector components
    '''
    
    # ---- scale pressure if plot is log-scaled
    p  = u.plev
    if(logp): p = np.log10(p)
    
    # --- build new pressure grid and interpolator
    pi     = np.linspace(np.min(p.values), np.max(p.values), len(p))
    uinterp = interp1d(p, u.values, axis=u.dims.index('plev'))
    vinterp = interp1d(p, v.values, axis=v.dims.index('plev'))
    
    # ---- decalre interpolated arrays and update coords
    if(logp): picoord = 10**pi
    else: picoord = pi
    ui, vi = xr.zeros_like(u), xr.zeros_like(v)
    ui = ui.assign_coords(plev=picoord)
    vi = vi.assign_coords(plev=picoord)
    
    # ---- interpolate
    ui.values, vi.values = uinterp(pi), vinterp(pi)
    
    # ---- do the same for significance variables if passed
    if(usig is not None and vsig is not None):
        uinterp = interp1d(p, usig.values, axis=u.dims.index('plev'))
        vinterp = interp1d(p, vsig.values, axis=v.dims.index('plev'))
        usigi, vsigi = xr.zeros_like(usig), xr.zeros_like(vsig)
        usigi = usigi.assign_coords(plev=picoord)
        vsigi = vsigi.assign_coords(plev=picoord)
        usigi.values, vsigi.values = uinterp(pi), vinterp(pi)
        return ui, vi, usigi, vsigi
    else:
        return ui, vi
    
# -----------------------------------------------------------------------

def integrate_tendency(tend, ic):
    '''
    Integrates a tendency variable across time, given an initial condition.
    The cumulative sum of the tendency data at time i is deposited at time i+i.
    Thus, the data at the first time position in the return is identical to 
    the initial condition.
    The input data will be trimmed such that any data occuring before the time
    of the initial condition is discarded.
    
    Parameters
    ----------
    tend : xarray DataArray
        the tendency data to integrate
    ic : xarray DataArray
        the initial condition. Must have a 'time' dimension of length 1
        
    Returns
    -------
    xarray DataArray
        the integrated tendency data
    '''
    
    # --- check inputs
    try: tic = np.atleast_1d(ic.time.values)
    except: raise RuntimeError('initial condition must have a \'time\' dimension of length 1')
    if(len(tic) > 1): 
        raise RuntimeError('initial condition should not have a time dimension!')
    
    # --- get initial time and initial condition
    tic = tic[0]
    ic  = ic.drop_vars('time')
    
    # --- trim the data after the integration initialization time
    tend = tend.sel(time=slice(tic, None))
    time = tend.time
    
    # --- integrate
    tend = spd * tend.cumsum('time').assign_coords(time=time)
    tend = ic + tend.shift(time=1, fill_value=0)
    return tend

# -----------------------------------------------------------------------

def plot_significance_filtered_contours(ax, x, y, var, levels, sigvar, sigcrit=0.05, **kwargs):
    '''
    Extracts all contours from vcontour which intersect at or are enclosed by at least one 
    contour from sigcontour, and plots it on the provided axis
    
    Parameters
    ----------
    ax : matplotlib axes object
        the axis on which to plot the result
    varcontours : matplotlib.contour.QuadContourSet
        the variable contours
    sigcontours : matplotlib.contour.QuadContourSet
        the significance contours
    fenc : float, optional
        fraction of variable contour that must be enclosed in significant regions in 
        order to keep. Default is 0.5, meaning that all variable contours for which
        half of their total length is within a significant region will be kept
    color : str, optional
        color for the plotted contours. Default is None, in which case the color 
        is inherited from varcontour
    lw : float, optional
        linewidth for the plotted contours. Default is None, in which case the linewidth 
        is inherited from varcontour
    ls : float, optional
        linestyle for the plotted contours. Default is None, in which case the linestyle
        is inherited from varcontour
    zorder : int, optional
        the zorder for the generated contours. Defaults to 0
    '''
    
    
    cg  = contour_generator(x, y, var)
    scg = contour_generator(x, y, sigvar)
    sigcontours = scg.lines(sigcrit)
    sigpaths    = [Path(sc) for sc in sigcontours]
    
    # manually close significance contours which are open (hit the plot edge)
    for i in range(len(sigpaths)):
        vertices = sigpaths[i].vertices
        if not np.array_equal(vertices[0], vertices[-1]):
            sigpaths[i] = Path(np.vstack([vertices, vertices[0]]))
    
    # find contour enclosures and intersections
    keep_contours = {}
    for level in levels:
        keep_contours[level] = []
        contours = cg.lines(level)
        for contour in contours:
            path = Path(contour)
            for sigpath in sigpaths:
                intersects = sigpath.intersects_path(path, filled=False)
                encloses = sigpath.contains_path(path)
                if(intersects or encloses): 
                    keep_contours[level].append(contour)
                    break
                
    # ----- plot
    for level in keep_contours.keys():
        if(level < 0): ls = '--'
        else: ls = '-'
        for contour in keep_contours[level]:
            xp, yp = contour.T[0], contour.T[1]
            ax.plot(xp, yp, ls=ls, **kwargs)
            
# -----------------------------------------------------------------------

def climatology(data, ctype='monthly'):
    '''
    Compute the climatological time average of the input data
    
    Parameters
    ----------
    data : xarray DataArray
        the data to average
    ctype : str
        The type of climatology to compute. Either:
        'monthly' (default) : computes a monthly climatology (return 
                            has integer time coordinate of length 12)
        'seasonal' : computes a seasonal climatology (return has a string
                                             time coordinate of length 4)
                                             
    Returns
    -------
    xarray DataArray
        the climatologically averaged data
    '''
    
    if(ctype == 'monthly'):
        climo = data.groupby('time.month').mean()
    if(ctype == 'seasonal'):
        month_len = data.time.dt.days_in_month
        month_len = month_len.groupby('time.season')
        weights   = month_len / month_len.sum()
        np.testing.assert_allclose(weights.groupby('time.season').sum().values, np.ones(4))
        data = data * weights
        climo = data.groupby('time.season').sum(dim='time')
    return climo

# -----------------------------------------------------------------------

def seasonal_mean(data):
    '''
    Computes the seasonal mean of the input data. The return will have two
    coordinates attached to the 'time' dimension; 'time' (mathching the 
    input time coordinate), and 'season'. Slicing or selecting in 'time'
    will also trim the 'season' coordinate. To slice or select directly
    on 'season' rather than 'time', run
    return_data.set_xindex('season')
    on the return array return_data
    
    
    Parameters
    ----------
    data : xarray DataArray
        the data to average

    Returns
    -------
    xarray DataArray
        the seasonally averaged data
    '''
    
    # get time coordinate
    time    = data.time
    
    # build year_season dimension so that we can group
    years   = data['time.year'].values
    months  = data['time.month'].values
    seasons = data['time.season'].values
    for i in range(len(years)): 
        if(months[i] == 12): years[i] += 1
    year_season = ['{}_{}'.format(y,s) for y,s in zip(years, seasons)]
    
    # get month lengths for weighting
    month_len = data.time.dt.days_in_month

    # assign year_season as the time coord for the data
    time      = time.assign_coords(year_season=('time', year_season))
    data      = data.assign_coords(year_season=('time', year_season))
    month_len = month_len.assign_coords(year_season=('time', year_season))
    
    # group data by season each year
    time          = time.groupby('year_season').first()
    data_sum      = (data * month_len).groupby('year_season').sum()
    month_len_sum = month_len.groupby('year_season').sum()
    
    # drop incomplete seasons
    counts           = month_len.groupby('year_season').count(dim='time')
    complete_seasons = counts == 3
    time          = time.where(complete_seasons).dropna(dim='year_season')
    data_sum      = data_sum.where(complete_seasons).dropna(dim='year_season')
    month_len_sum = month_len_sum.where(complete_seasons).dropna(dim='year_season')
    
    # compute seasonal averages
    data = data_sum / month_len_sum
    
    # replace the year_season coord with 'time' (cftime datetimes, corresponding to the first
    # month of each season)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        data = data.rename(year_season='time')
    data = data.assign_coords(time=('time', time.values))
    
    # annoyingly, during the operations which grouped the data by year_season, the data groups
    # are sorted alphabetically by the year_season string labels. Undo this by sorting in time
    data = data.sortby('time')
    
    # finally, add a "season" coordinate
    data = data.assign_coords(season=('time', data.time.dt.season.values))
    
    return data

# -----------------------------------------------------------------------

def remove_troposphere(x, trop, buffer=0):
    '''
    Removes (sets to nan) all data below the tropopause for the input variable
    
    Parameters
    ----------
    x : xarray DataArray
        the input data with a pressure coordinte 'plev' in hPa
    trop : xarray DataArray
        the tropopause pressure position in Pa
    buffer : float, optional
        extra puffer around tropopause to remove data, in hPa
        If positive, then extra data above the tropopause is removed
        If negative, then extra data below the tropopause is retained
        Defaults to zero
    '''
    return x.where(trop/100-buffer >= x.coords['plev'])

# -----------------------------------------------------------------------

def confidence_interval(x, std, N, alpha):
    '''
    Computes confidence interval bounds on a mean quantity x
    
    Parameters
    ----------
    x : xarray DataArray
        the mean quantity
    std : xarray DataArray
        the standard devation of the mean quantity
    N : int
        the number of ensemble members (dof is N-1)
    alpha : float
        two-sided threshold in p-value at which to find t-crit
        
    Return
    ------
    2-tuple of xarray DataArray
        the lower and upper confidence interval bounds
    '''
    
    tcrit = scipy.stats.t.ppf(1-alpha/2, N-1)
    return x - (tcrit * std/np.sqrt(N)), x + (tcrit * std/np.sqrt(N))

# -----------------------------------------------------------------------

def to_datetime(times):
    '''
    Converts an DataArray of cftime.datetime objects to datetime objects, which
    matplotlib prefers
    
    Parameters
    ----------
    times : cftime.datetime object, or array of cftime.datetime object
        the time(s) to convert
    
    Returns
    -------
    datetimes : array of datetime objects
        the converted times
    '''
    if(type(times) == type(xr.DataArray())): times=times.values
    if(not is_array(times)): times = [times]
    datetimes = [datetime(date.year, date.month, date.day) for date in times]
    if(len(datetimes) == 1):
        return datetimes[0]
    else:
        return datetimes

# -----------------------------------------------------------------------

def season_timeticks(ax, times, option, which='x', minor=True,
                     include_year=False, year_delimiter=' ', year_on_jan_only=False):
    '''
    Changes the ticks on an input time axis so that they are drawn at specified months
    
    Parameters
    ---------
    taxis : matplotlib axex object
        the axes containing the time axis to modify
    times : array-like of cftime.datetime or datetime objects
        the time coordiantes of the data being plotted
    option : str, or list
        Which ticks to draw. Options are:
        - 'month'    : draw at DJF MAM JJA SON (every month)
        - 'season'   : draw at J A J O         (center month of every season)
        - 'solstice' : draw at J J             (center month of summer, winter)
        - 'equinox'  : draw at A O             (center month of spring, fall)
        - 'winter'   : draw Jan
        - 'summer'   : draw July
        - 'year'     : draw years as integers on January tick, with no month name
    which : str, optional
        'x' or 'y', giving which axis on the plot is the time axis. 
        Deafaults to 'x'
    minor : bool, optional
        whether or not to draw minor ticks on every month, in the case that
        the option arg is not "month"
    include_year : bool, optional
        whether or not to include the year in the tick string. 
        If option='year', this does nothing.
        Defaults to False.
    year_delimiter : str,m optional
        separator string for the month and year name, if include_year=True.
        Defaults to a space.
    year_on_jan_only : bool
        if include_year is True, then setting this option to True will 
        print the year only on January ticks, rather than every major tick
    '''

    # define January month-year formatter
    class JanDateFormatter(mdates.DateFormatter):
        def __call__(self, x, pos=None):
            date = mdates.num2date(x)
            return f"{date.strftime('%b')}{year_delimiter + date.strftime('\'%y') if date.month == 1 else ''}"

    
    # define ticks
    all_months = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    months = {'month'   : all_months, 
              'season'  : [1, 4, 7, 10],
              'solstice': [1, 7],
              'equinox' : [4, 10],
              'winter'  : [1],
              'summer'  : [7],
              'year'    : [1]}
    try:
        months = months[option]
    except TypeError:
        months = option
    
    # configure axis
    assert which in ['x', 'y'], "arg \'which\' must be \'x\' or \'y\'"
    if(which=='x'):
        taxis = ax.xaxis
        tlim = ax.get_xlim()
        limset = ax.set_xlim
    elif(which =='y'): 
        taxis = ax.yaxis
        ax.get_ylim()
        limset = ax.set_ylim
    
    # set ticks
    ticks = []
    minor_ticks = []
    for year in np.unique([t.year for t in times]):
        for month in months:
            ticks.append(datetime(year, month, 1))
        for month in all_months:
            minor_ticks.append(datetime(year, month, 1))
    taxis.set_ticks(ticks)
    if(minor_ticks!=ticks and minor==True):
        taxis.set_ticks(minor_ticks, minor=True)
    limset(tlim)
    
    # add year strings if requested
    if(option=='year'):
        date_format = DateFormatter("%Y")
    else:
        if(include_year):
            if(year_on_jan_only): date_format = JanDateFormatter(None)
            else: date_format = DateFormatter("%b{}'%y".format(year_delimiter))
        else:
            date_format = DateFormatter("%b")
    taxis.set_major_formatter(date_format)    
    
# -----------------------------------------------------------------------
    
def format_paxis(ax, ticks=None):
    '''
    Formats an the pressure axis of an input axis for inverted, log-scaled, scalar formatted ticks
    
    Parameters
    ----------
    ax : matplotlib.axis object
        the axes to modify. It is assumed that the pressure axis is the y-axis
    '''
    
    if(ticks is None):
        ticks=[1000, 300, 100, 30, 10, 3, 1, 0.3, 0.1, 0.03, 0.01]
        tick_str = []
    
    lim = sorted(ax.get_ylim())[::-1]
    ax.invert_yaxis()
    ax.set_yscale('log')
    ax.yaxis.set_major_formatter(ScalarFormatter())
    ax.yaxis.set_ticks(ticks)
    ax.set_ylim(lim)
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: \
                              f'{x:.3f}'.rstrip('0').rstrip('.')))
    
# -----------------------------------------------------------------------

def format_lataxis(ax, which='x', ticks=None):
    '''
    Format ticks on a latitude axis, where positive numbers X, 
    negative numbers -X, and zeros, respectively, become
    X°N, X°S, and 0°
    
    Parameters
    ----------
    ax : matplotlib axes object
        the axes containing the latitude axis
    which : str
        which axis to format. Either 'x' or 'y'
    ticks : list of float
        the tick positions in degrees. Default is None, in which
        case the tick positions present on the input axis are used
    ''' 
    if(which == 'x'):
        if(ticks is not None): ax.xaxis.set_ticks(ticks)
        ax.xaxis.set_major_formatter(ticker.FuncFormatter(format_longitude))
    elif(which == 'y'):
        if(ticks is not None): ax.yaxis.set_ticks(ticks)
        ax.yaxis.set_major_formatter(ticker.FuncFormatter(format_longitude))
def format_longitude(x, pos):
    if x > 0: 
        if(int(x) == x): return '{}'.format(x).split('.')[0] + '°N'
        else:            return '{}'.format(x).rstrip('0') + '°N'
    elif x < 0: 
        if(int(x) == x): return '{}'.format(-x).split('.')[0] + '°S'
        else:            return '{}'.format(-x).rstrip('0') + '°S'
    else: return '0°'

# -----------------------------------------------------------------------

def format_ticks(ax, x='bottom', y='left'):
    '''
    Formats tick marks so that they appear on all sides of the plot, and
    positions the tick labels on two specified sides
    
    Parameters
    ----------
    ax : matplotlib axes object
        the axes containing the latitude axis
    x : str, optional
        where to render the x-axis labels. Either 'bottom' (default) or 'top'
    y : str, optional
        where to render the y-axis labels. Either 'left' (default) or 'right'
    '''
    
    targs = {'top':True,'bottom':True,'left':True,'right':True,'which':'both',
             'labeltop':False,'labelbottom':False,'labelleft':False,'labelright':False}
    if(y == 'left'):     targs['labelleft'] = True
    elif(y == 'right'):  targs['labelright'] = True
    elif(y is not None): raise RuntimeError('argument \'y\' must be \'left\', \'right\', or None')
    if(x == 'top'):      targs['labeltop'] = True
    elif(x == 'bottom'): targs['labelbottom'] = True
    elif(x is not None): raise RuntimeError('argument \'x\' must be \'left\', \'right\', or None')
    ax.tick_params(**targs)
    if(y == 'left'):     ax.yaxis.set_label_position('left')
    elif(y == 'right'):  ax.yaxis.set_label_position('right')
    else:                ax.yaxis.set_label('')
    if(x == 'top'):      ax.xaxis.set_label_position('top')
    elif(x == 'bottom'): ax.xaxis.set_label_position('bottom')
    else:                ax.xaxis.set_label('')

# -----------------------------------------------------------------------

def get_cmap_norm(levels, norm):
    '''
    Construct matplotlib color normalization objects from keywords

    Parameters
    ----------
    levels : numpy array or xarray DataArray
        the contour levels
    norm : str
        string specifying the type of norm. Options are:
        'linear    : linear colormap
        'twoslope' : separate linear color scales are used on each side 
                     of zero, to utilize the full range 
        'uneven'   : the levels are non-equidistant, and colors should 
                     be sampled evenly for each level, reagrdless of the 
                     separation between levels
        'log'      : logarithmically space the colors
        'symlog'   : logarithmically space the colors symmetrically on 
                     each side of zero
    '''
    
    if(norm == 'linear' or norm == 'twoslope'):
        low, high = min(levels), max(levels)
        if(low < 0 and high > 0):
            norm = colors.TwoSlopeNorm(vmin=min(levels), vcenter=0, vmax=max(levels))
        else:
            norm = None
    elif(norm == 'log'):
        norm = colors.LogNorm(vmin=min(levels), vmax=max(levels))
    elif(norm == 'symlog'):
        norm = colors.SymLogNorm(linthresh=min(levels[levels>0]), linscale=1)
    elif(norm == 'uneven'):
        norm = colors.BoundaryNorm(levels, 256, extend='both')
    else:
        raise RuntimeError('norm type \'{}\' not supported'.format(norm))
    return norm

# -----------------------------------------------------------------------

def give_cmap_white_center(cmap_name, center_gap=0.2, center_width=0.15, 
                           steepness=2, demo=False):
    '''
    Smoothly places white at the center of a colormap
    
    Parameters
    ----------
    cmap_name : str
        colormap name
    center_gap : float, optional
        fraction to cut out of the colormap center, in color space. 
        For eaxmple, if center_gap=0.1, then colors bookmarking the 
        start and end of the white fade will be sampled from the colormap
        at (0.5-center_gap/2) and (0.5+center_gap/2).
        For eaxmple, if center_gap=0.1, then the colors between 45% and 55%
        are removed fom the colormap.
        Defaults to 0.2
    center_width : float, optional
        width of fade to and from white, as a fraction of the total colormap
        width. Defaults to 0.15
    steepness : int, optional
        the steepness of the fade to and from white in the colormap center.
        This is an integer that decides the structure of the sequence
        [left_color, (white)*steepness, right_color]
        which is used for interpolation via LinearSegmentedColormap(). Higher
        values will result in a steeper fade to/from white. Defaults to 2.
    demo : bool, optional
        whether or not to demo the colormap, by rendering a plot and calling
        plt.show(). Defaults to False.
    '''
    
    # --- get left and right portions of colormap
    cmap=plt.get_cmap(cmap_name)
    left=cmap(np.linspace(0,0.5-center_gap/2,256))
    right=cmap(np.linspace(0.5+center_gap/2,1,256))

    # --- construct the center fade from each inner color, and white
    center_color = list(mpl.colors.to_rgb('w'))
    center = LinearSegmentedColormap.from_list('left_extension', 
             [left[-1]] + [center_color]*steepness + [right[0]], N=(256*2)*center_width)
    
    # --- combine the left, center, and right colormap portions
    cmap_out = np.vstack([left, center(np.linspace(0,1,center.N)), right])
    cmap_out = LinearSegmentedColormap.from_list('{}_with_white_center'.format(cmap_name), cmap_out)

    # --- plot demo
    if(demo):
        fig, ax  = plt.subplots(2, 1, figsize=(6, 2))
        gradient = np.linspace(0, 1, 512).reshape(1, -1)
        for axi in ax: axi.set_axis_off()
        ax[0].imshow(gradient, aspect='auto', cmap=cmap)
        ax[1].imshow(gradient, aspect='auto', cmap=cmap_out)
        ax[0].set_title(cmap_name, fontsize=11)
        ax[1].set_title('{} with white center added'.format(cmap_name), fontsize=11)
        plt.tight_layout()
        plt.show()
        
    return cmap_out

# -----------------------------------------------------------------------

def cbarfmt(v, pos):
    # Format the float with up to 2 decimal places
    vfmt = "{:.2f}".format(v)
    # Remove trailing zeros and the decimal point if not needed
    vfmt = vfmt.rstrip('0').rstrip('.')
    # if after this formatting, the number is zero, must have been smaller than 0.01
    if(float(vfmt) == 0 and v != 0):
        vfmt = '{:.2e}'.format(v).replace('e-0','e-').replace('e+0','e+')
        # Remove trailing zeros and the decimal point if not needed
        val, exp = vfmt.split('e')
        val = val.rstrip('0').rstrip('.')
        vfmt = val + 'e' + exp
    #done
    return vfmt
def pretty_format_colorbar():
    '''
    Returns a colorbar FuncFormatter, which can be used as the "format" input
    argument to matplotlib.figure.Figure.colorbar, which formats all values grater 
    than 0.01 as decimal values with two decimal digits, with trailing zeros removed.
    All values less than 0.01 are formatted in scientific notation, with leading 
    zeros on the exponent removed, and trailing zeros on the value removed.

    Parameters
    ----------
    fmt : str
        string format for values. If not supplied, 
    '''
    return(FuncFormatter(cbarfmt))
            
# -----------------------------------------------------------------------

def log_quiver(u, v):
    '''
    Log-scale vector components, maintaining sign
    
    Parameters
    ----------
    u : array-like
        the horizontal vector component
    v : array-lioke
        the vertical vector component
    '''
    length = np.hypot(u, v)
    min_length = np.min(length.where(length!=0, drop=True))
    unorm = u/length
    vnorm = v/length
    ulog = np.log10(length/min_length) * unorm
    vlog = np.log10(length/min_length) * vnorm
    return ulog, vlog

# -----------------------------------------------------------------------

def fake_streamplot(contours, ax, spacing=50, length=30, min_length=10, color='k', lw=1, head_width=3, head_length=3):
    
    collections = contours.collections
    for collection in collections:
        for path in collection.get_paths():
            # Get vertices of the contour line
            vertices = path.vertices
            vcodes    = path.codes
            if len(vertices) < 2:  # Skip too-short paths
                continue

            # Compute tangential vectors
            for i in range(1, len(vertices)):
                if(i%spacing != 0):
                    continue

                # get a section of the contour before the arrow
                path_len = length
                path_min = min_length
                path_len = min(len(vertices)-i, path_len)
                if(path_len < path_min): continue
                verts = vertices[i-1:i+path_len]
                codes  = vcodes[i-1:i+path_len]
                verts = verts[:next((i for i, code in enumerate(codes) if code not in [1,2]), len(verts))]
                codes = codes[:next((i for i, code in enumerate(codes) if code not in [1,2]), len(codes))]
                codes[0] = Path.MOVETO

                # Place an arrow
                new_path = Path(verts, codes)
                patch = patches.PathPatch(new_path, facecolor='none', lw=lw)
                ax.add_patch(patch)

                p2 = vertices[i - 1]
                p1 = vertices[i]
                dx, dy = p2 - p1

                # Normalize the direction
                mag = np.sqrt(dx**2 + dy**2)
                if length == 0:
                    continue
                dx /= mag
                dy /= mag

                #ax.arrow(
                #    p1[0], p1[1], dx * scale, dy * scale,  # Scale arrow length
                #    head_width=head_width, head_length=head_length, fc=color, ec=color
                #)
                arrow = patches.FancyArrowPatch((0.2, 0.2), (0.8, 0.8),
                               transform=ax.transAxes,  # Use axes coordinates
                               arrowstyle="->",
                               color="red")
                #ax.quiver([p1[0]], p1[1], dx, dy, angles='xy', scale=1, width=0.08, pivot='tip', headwidth=head_width, headlength=head_length)

# -----------------------------------------------------------------------

def is_array(obj):
    '''
    Check if input object is an array-like object
    '''
    return hasattr(obj, '__len__') and not isinstance(obj, (str, bytes))

# -----------------------------------------------------------------------

def symlog(data):
    '''
    Take a symmetric logarithm of the input data, i.e. scale the data on each
    side of zero
    
    Parameters
    ----------
    data : array-like
        the data to log-scale
    '''
    data = np.sign(data) * np.log10(np.abs(data))
    return data

# -----------------------------------------------------------------------

def symlogspace(a, b, s=None):
    '''
    Return a symmtric logspace, i.e. a signed logspace on each side of zero
    
    Parameters
    ----------
    data : array-like
        the data to log-scale
    '''
    p = np.logspace(a, b, s)
    n = -np.logspace(a, b, s)[::-1]
    return np.hstack([n, [0], p])