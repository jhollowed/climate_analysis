import sys
import pdb
import glob
import Nio, Ngl
import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from cftime import DatetimeNoLeap


# ===============================================================================================
# ===============================================================================================


def isolate_barnes_period_data():
    '''
    Gather data for all CMIP6 variables in the period of 1980-2004, write out to 
    file per-variable

    Error messages containing 'Killed' means the progam has run out of memory; try
    again on an interactive slurm job
    '''

    # gather metadata from file names
    realizations = np.array([int(f.split('_r')[-1].split('i')[0]) for f in CMIP_FILES])
    variables = np.array([f.split('/')[-1].split('_')[0] for f in CMIP_FILES])
    var_names = np.unique(variables)
    
    # construct mask on files which contain any data in the pinatubo period
    # (data corresponding to the precise period will be extracted from these files later)
    times = [f.split('/')[-1].split('_')[-1].split('.nc')[0].split('-') for f in CMIP_FILES]
    for time in times:
        time[0] = int(time[0][:-2])
        time[1] = int(time[1][:-2])
    time_mask = [1 in np.searchsorted(time, PINATUBO_PERIOD) for time in times]

    # loop over realiozations, variables
    for i in range(3):
        rlz_mask = realizations == i+1
        print('working on realization {}'.format(i+1))
        
        for j in range(len(var_names)):
            print('gathering data for variable {}'.format(var_names[j]))
            
            # identify files belogning to this variable, this realization, relevant time period
            var_mask = variables == var_names[j]
            file_mask = rlz_mask & var_mask & time_mask
            
            # this sort is important;
            # time slice bounds need to be unique values present in the dataset unless 
            # the sliced dimension is sorted monotonically. Sorting the data after reading,
            # which will require lots of memory, or we can just mak sure that we read them 
            # in in temporally ascending order. Here we assume that the filenames are 
            # identical except for the date range at the end, in which case this is all we 
            # need to do
            files = sorted(CMIP_FILES[file_mask])
 
            # concatenate contents of relevant time intervals to dataset...
            dsets = np.empty(len(files), dtype=xr.core.dataset.Dataset)
            for k in range(len(files)):
                dsets[k] = xr.open_dataset(files[k])
            var_alltime = xr.concat(dsets, dim='time')
           
            # isolate Pinatubo period
            var_out = var_alltime.sel(time=slice(str(PINATUBO_PERIOD[0]), str(PINATUBO_PERIOD[1])))
            
            # write out NetCDF file
            start = var_out['time'].values[0]
            end = var_out['time'].values[-1]
            startstr = '{}{}'.format(start.year, f'{start.month:02}')
            endstr = '{}{}'.format(end.year, f'{end.month:02}')
            dest = '{}/{}_Amon_E3SM-1-0_{}_r{}_{}-{}.nc'.format(PINATUBO, var_names[j], MODEL_TYPE,
                                                                i+1, startstr, endstr) 
            print('writing out to {}'.format(dest))
            var_out.to_netcdf(dest, format="NETCDF4_CLASSIC")

        
# ===============================================================================================


def compute_climatology():
    '''
    Computes the monthly climatology of the CMIP6 data for the pinatubo period, and writes out to
    a NetCDF file
    '''
    
    # gather metadata from file names
    realizations = np.array([int(f.split('/')[-1].split('_r')[-1].split('_')[0]) for f in PINATUBO_FILES])
    variables = np.array([f.split('/')[-1].split('_')[0] for f in PINATUBO_FILES])
    var_names = np.unique(variables)

    for i in range(3):
        i += 1
        rlz_mask = realizations == i
        print('working on realization {}'.format(i))
        
        for j in range(len(var_names)):
            print('computing monthly climatology for variable {}'.format(var_names[j]))
            
            # get this variable for this realization, open file
            var_mask = variables == var_names[j]
            file_mask = rlz_mask & var_mask
            files = PINATUBO_FILES[file_mask]
            if(len(files) != 1):
                raise RuntimeError('More than one file found for {} in realization r{}:\n{}'.format(
                                    var_names[j], i, files))
            dset = xr.open_dataset(files[0])
           
            # compute climatology, write out
            dest = '{}{}_climatology_{}'.format(files[0].split(var_names[j])[0],
                                               var_names[j], files[0].split('Amon_')[-1]) 
            climatology = dset.groupby('time.month').mean('time')
            climatology.to_netcdf(dest, format='NETCDF4_CLASSIC')


# ===============================================================================================
            

def compute_anomalies():
    '''
    Compute the monthly anomalies of variable x with respect to the mean climatology.
    Write out result per each variable, per each ensemble member, and the mean of the
    three members per each variable
    '''
    
    # gather metadata from file names
    realizations = np.array([int(f.split('/')[-1].split('_r')[-1].split('_')[0]) for f in 
                             PINATUBO_FILES])
    variables = np.array([f.split('/')[-1].split('_')[0] for f in PINATUBO_FILES])
    var_names = np.unique(variables)

    for j in range(len(var_names)):
        print('computing monthly anomalies for variable {}'.format(var_names[j]))
        for i in range(3):
            i += 1
            rlz_mask = realizations == i
            print('working on realization {}'.format(i))
            
            # get this variable for this realization, and climatology
            var_mask = variables == var_names[j]
            file_mask = rlz_mask & var_mask
            files = PINATUBO_FILES[file_mask]
            if(len(files) != 2):
                raise RuntimeError('Expected 2 files, climatology and monthly average time series,'\
                                   'for {} in realization r{}, but found {}:\n{}'.format(
                                    var_names[j], i, len(files), files))
            
            # read in
            mask = np.array(['climatology' in f for f in files])
            fclimatology = files[mask][0]
            fmonthly = files[~mask][0]
            climatology = xr.open_dataset(fclimatology)
            monthly = xr.open_dataset(fmonthly)
            
            pdb.set_trace()

            # compute anomalies, update mean
            anomaly = monthly - climatology
            
            # write out anomaly
            monthly.close()
            climatology.close()
            dest = '{}{}_anomaly_{}'.format(fmonthly.split(var_names[j])[0],
                                            var_names[j], fmonthly.split('Amon_')[-1])
            anomaly.to_netcdf(dest, format='NETCDF4_CLASSIC')

            # update ensemble mean
            if(i == 1):
                ensemble_mean = anomaly
            else:
                ensemble_mean += anomaly
                if(i == 3):
                    ensemble_mean /= 3
                    dest = '{}{}_meanAnomaly_{}'.format(fmonthly.split(var_names[j])[0], var_names[j], 
                                                        fmonthly.split('_')[-1])
                    ensemble_mean.to_netcdf(dest, format='NETCDF4_CLASSIC')
                    
             
# ===============================================================================================
# ===============================================================================================


if __name__ == '__main__':
    
    MODEL_TYPE = sys.argv[1] # either 'AMIP' or 'historical'
    CMIP6 = '/nfs/turbo/cjablono2/hollowed/CMIP6/E3SM_{}/pinatubo_period'.format(MODEL_TYPE)
    PINATUBO = '/nfs/turbo/cjablono2/hollowed/CMIP6/E3SM_{}/pinatubo_period_exact'.format(MODEL_TYPE)
    CMIP_FILES = np.array(glob.glob('{}/*.nc'.format(CMIP6)))
    PINATUBO_FILES = np.array(glob.glob('{}/*.nc'.format(PINATUBO)))
    PINATUBO_PERIOD = [1980, 2004]

    if(len(glob.glob('{}/*'.format(PINATUBO))) == 0):
        isolate_barnes_period_data()
    if(len(glob.glob('{}/*climatology*'.format(PINATUBO))) == 0):
        compute_climatology()
    if(len(glob.glob('{}/*anomaly*'.format(PINATUBO))) == 0):
        compute_anomalies()
