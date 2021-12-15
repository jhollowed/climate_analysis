import pdb
import glob
import Nio, Ngl
import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from cftime import DatetimeNoLeap

CMIP6 = '/nfs/turbo/cjablono2/hollowed/CMIP6/E3SM_AMIP/pinatubo_period'
FILES = np.array(glob.glob('{}/*.nc'.format(CMIP6)))
PINATUBO_PERIOD = [1980, 2004]


# ===============================================================================================
# ===============================================================================================


def isolate_barnes_period_data():
    '''
    Gather data for all CMIP6 variables in the period of 1980-2004, write
    out to single NetCDF file
    '''

    # create empty array of netCDF variables as xarray DataSets
    realizations = np.array([int(f.split('amip_r')[-1].split('i')[0]) for f in FILES])
    variables = np.array([f.split('/')[-1].split('_')[0] for f in FILES])
    var_names = np.unique(variables)

    # cast pinatubo period to cftime objects
    pt = [DatetimeNoLeap(year, 1, 1) for year in PINATUBO_PERIOD]
   
    # loop over realiozations, variables
    for i in range(3):
        rlz_mask = realizations == i+1
        allvars_alltime = np.empty(len(var_names), dtype=xr.core.dataset.Dataset)
        
        for j in range(len(var_names)):
            print('gathering data for variable {}'.format(var_names[j]))
            
            # for each realization, concatenate contents of time intervals to dataset
            var_mask = variables == var_names[j]
            file_mask = np.logical_and(rlz_mask, var_mask)
            files = FILES[file_mask]
            
            dsets = np.empty(len(files), dtype=xr.core.dataset.Dataset)  
            for k in range(len(files)):
                dsets[k] = xr.open_dataset(files[k])
            var_alltime = xr.concat(dsets, dim='time')
           
            # isolate Pinatubo period, and append this dataset to collection of all datasets 
            # for this realization
            # ...
            # apparently xarray time slices in cftime format need to be exact;
            # getting the earliest/latest dates appearing in this dataset in the pinatubo period:
            start = sorted(var_alltime['time'].values[var_alltime['time.year'] == 1980])[0]
            end = sorted(var_alltime['time'].values[var_alltime['time.year'] == 2004])[-1]
            allvars_alltime[j] = var_alltime.sel(time=slice(start, end))
        
        # write out one NetCDF file per this realization
        dest = '{}/allVars_pinatuboPeriod_Amon_E3SM-1-0_amip_r{}.nc'.format(CMIP6, i+1) 
        print('writing out to {}'.format(dest))
        out = xr.merge(allvars_alltime)
        out.to_netcdf(dest, format="NETCDF4_CLASSIC")

        
# ===============================================================================================



            

        
# ===============================================================================================
# ===============================================================================================


if __name__ == '__main__':

    if(len(glob.glob('{}/allVars*'.format(CMIP6))) == 0):
        isolate_barnes_period_data()
            
                
            
