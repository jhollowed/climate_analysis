'''
This module provides various functions for computing various climate metrics, indices, and
performing other common tasks on climate or weather data
'''

import os
import sys
import pdb
import glob
import numpy as np
import xarray as xr
import Nio, Ngl
from cftime import DatetimeNoLeap

# ==========================================================================================

def check_data_inputs(data):
    acceptable_types = [xr.core.dataset.Dataset, xr.core.dataarray.DataArray, str]
    if(type(data) not in acceptable_types):
        raise RuntimeError('data must be provided as an xarray Dataset, DataArray, '\
                           'or a string to the file')
    if(type(data) == str):
        data = xr.open_dataset(data)
    return data


# -------------------------------------------------------------


def compute_climatology(data, ncdf_out=None):
    '''
    Compute the monthly climatology on the data contained in file f

    Parameters
    ----------
    data : xarray Dataset, xarray DataArray, or string
        An xarray object containing the data, or path to the file as a string. 
        Must include times with at least monthly resolution.
    ncdf_out : string, optional
        The file to write the resulting climatology data out to. Default is None,
        in which case nothing is written out, and the result is instead returned 
        as an xarray object.

    Returns
    -------
        climatology : xarray data object
            The monthly climatology from the input data
    '''

    # check inputs
    data = check_data_inputs(data)

    # compute climatology 
    climatology = data.groupby('time.month').mean('time')
    if(ncdf_out is not None):
        climatology.to_netcdf(ncdf_out, format='NETCDF4')
    else:
        return climatology


# -------------------------------------------------------------


def QBO_index(data, p=30, u_var_name='ua', p_var_name='plev'):
    '''
    Computes the QBO index as a time series from the input data. Note, lat/lon
    dimensions are assumed to be named exactly 'lat', 'lon'
    
    Parameters
    ----------
    data : xarray Dataset, xarray DataArray, or string
        An xarray object containing the data, or path to the file as a string. 
    p : float, optional
        Pressure level at which to take the zonal-mean zonal wind, in hPa. The standard
        choices for this index at 30 hPa (the default), and 50 hPa. Defaults to 30hPa. If
        this pressure level not present in the data, it will be linearly interpolated
    u_var_name : string
        Varibale name to expect for the zonal wind. Defaults to the CMIP6 default, which is 'ua'
    p_var_name : string
        Varibale name to expect for the pressure level. Defaults to the CMIP6 default, 
        which is 'plev'

    Returns
    -------
    u_qbo : DataArray
        A 1D time series of the QBO index
    '''

    data = check_data_inputs(data)

    # define lat range for meridional averaging (QBO isolation) as +- 5 deg
    latslice = slice(-5, 5)
    
    # open data, select on meridional range and average, take zonal mean
    u = data[u_var_name].sel({'lat':latslice})
    u = u.mean('lat')
    u = u.mean('lon')

    # take zonal wind at indexing pressure level (in Pa)
    p = p*100
    if(p in u['plev']):
        u_qbo = u.sel({p_var_name:3000})
    else:
        u_qbo = u.interp({p_var_name:3000})
    return u_qbo
    
    
# -------------------------------------------------------------


def Nino34_index(data, sst_var_name='tos', time_samples_mon=1):
    '''
    Computes the Nin 3.4 index as a time series from the input data. Note, lat/lon 
    dimensions are assumed to be named exactly 'lat', 'lon'
    
    Parameters
    ----------
    data : xarray Dataset, xarray DataArray, or string
        An xarray object containing the data, or path to the file as a string. 
    sst_var_name : string
        Varibale name to expect for the sea surface temperature. Defaults to 
        the CMIP6 default, which is 'tos'
    time_samples_mon : int
        Number of time sample represented in the data per month (to facilitate
        5 month rolling average). Defaults to 1, in which case the data is
        assumed to be monthly. If e.g. the data is daily, this arg should be 30

    Returns
    -------
    sst_nino : DataArray
        A 1D time series of the El Nino 3.4 index
    '''

    data = check_data_inputs(data)

    # define index averaging region as [5N-5S, 170W-120W]
    latslice = slice(-5, 5)
    lonslice = slice(120, 170)

    # take unweighted average SST in the region (unweighted since region is meridionally small)
    data_region = data.sel({'lat':latslice, 'lon':lonslice})
    sst = data_region[sst_var_name]
    sst = sst.mean('lat')
    sst = sst.mean('lon')
    
    # take 5 month running average
    sst_nino = sst.rolling(time=time_samples_mon).mean()
    return sst_nino
    

# -------------------------------------------------------------
    
    


















