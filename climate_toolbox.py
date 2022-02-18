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
from cftime import DatetimeNoLeap

sys.path.append('/glade/u/home/jhollowed/repos/ncl_exports')
from wrappers import dpres_hybrid_ccm

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
    '''

    data = check_data_inputs(data)

    # define lat range for meridional averaging (QBO isolation) as +- 5 deg
    latslice = slice(-5, 5)
    
    # open data, select on meridional range and average, take zonal mean
    u = data['u'].sel({'lat':latslice})
    u = u.mean('lat')
    u = u.mean('lon')

    # take zonal wind at indexing pressure level (in Pa)
    p = p*100
    if(p in u['plev']):
        u = u.sel({'plev':3000})
    else:
        u = u.interp({'plev':3000})
    return u
    
    
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
    '''

    data = check_data_inputs(data)

    # define index averaging region as [5N-5S, 170W-120W]
    latslice = slice(-5, 5)
    lonslice = slice(120, 170)

    # take unweighted average SST in the region (unweighted since region is meridionally small)
    data_region = (xr.open_dataset(data)).sel({'lat':latslice, 'lon':lonslice})
    sst = data_region[sst_var_name]
    sst = sst.mean('lat')
    sst = sst.mean('lon')
    
    # take 5 month running average
    sst = sst.rolling(time=time_samples_mon).mean()
    return sst
    

# -------------------------------------------------------------
    
    
def l2_norm(dat, varname='U', norm_type=15, compdat=None, compMean=False):
    '''
    Computes l2 norms from Jablonowski+Williamson 2006

    Parameters
    ----------
    dat : xarray Dataset
        the input data
    varname : string
        Variable to take the norm of
    norm_type: int
        Either 14 or 15, giving the equation number from JW06
        - Default is 15, in which case the norm is taken between the 
        zonal mean of the field, and the initial zonal mean of the 
        field (this measures the departure from stability)
        - If 14, the norm is taken between the zonal mean of the field, 
        and the 3D field (this measures the departure from symmetry)
    compdat : xarray Dataset
        Another set of input data to compare the input against.
        Default is None, in which case the norm is taken as described
        in Jablonowski+06, respecting norm_type. 
        If provided, norm is taken between the entire 3D fields of each
        file as a time series (i.e. norm_type will be ignored)
        It will be assumed that the comparison data set is on the same 
        horizontal and vertical grid as the input (horizontal only if looking
        at 2D quantities), and output at the same timesteps
    compMean : bool
        If True, and compdat is given, copare the zonally averaged fields
        from each dataset, rather than the 3D field
    '''

    var = dat[varname]
    u = dat['U']
    ps = dat['PS']
    hyai = dat['hyai']
    hybi = dat['hybi']
    P0 = 100000.0         #Pa, to match ps units
    dp = dpres_hybrid_ccm(ps, P0, hyai, hybi).values
    dp = xr.DataArray(dp, dims=u.dims, coords=u.coords, name='dp')
    
    # approx horizontal weighting
    rad = np.pi/180
    lat = dat['lat']
    wy  = np.cos(lat*rad) #proxy [ sin(j+1/2)-sin(j-1/2) ]

    # get time samples and zonal-mean var
    ntime = len(dat['time'])
    varm = var.mean('lon')

    # read comparison dataset if given
    if(compdat is not None):
        norm_type = 0
        compvar = compdat[varname]
        comp_varm = var.mean('lon')
        
    # compute the norm...
    
    # =============== EQ 14 ===============
    if(norm_type == 14):
        
        norm = np.zeros(ntime)
        for i in range(ntime):
            # the broadcasting here works by matching dimensions by their names, 
            # which importantly comes from the input netcdf file
            diff2 = (var[i] - varm[i])**2
            num = np.sum(diff2 * wy * dp[i])
            den = np.sum(wy * dp[i])
            norm[i] = np.sqrt(num/den)    
        
    # =============== EQ 14 ===============
    elif(norm_type == 15):
        
        norm = np.zeros(ntime)
        for i in range(ntime):
            # the broadcasting here works by matching dimensions by their names, 
            # which importantly comes from the input netcdf file
            diff2 = (varm[i] - varm[0])**2
            num = np.sum(diff2 * wy * dp[i])
            den = np.sum(wy * dp[i])
            norm[i] = np.sqrt(num/den)
    
    # =============== COMPARE ===============
    elif(norm_type == 0):
        
        norm = np.zeros(ntime)
        for i in range(ntime):
            # the broadcasting here works by matching dimensions by their names, 
            # which importantly comes from the input netcdf file
            if(compMean):
                diff2 = (varm[i] - compvarm[i])**2
            else:
                diff2 = (var[i] - compvar[i])**2
            num = np.sum(diff2 * wy * dp[i])
            den = np.sum(wy * dp[i])
            norm[i] = np.sqrt(num/den)
    
    return norm



