'''
Joe Hollowed 2024

This script performs zonal averaging of a single limvar ensemble member. 
Command line arguments are:

Usage
-----
python ./limvar_zonalmeans.py [Nens] [histnum] [massMag] [tmini] [tmaxi] [overwrite] [dry]

Parameters
----------
Nens : int
    the ensemble member specified as an integer
histnum : int
    history file number to use
    0 = monthly data
    1 = daily data
massMag : int
    the mass magnitude of the ensemble. Valid values are:
    0, 1, 3, 5, 7, 10, 13, 15
    if massMag = 0, this sepcifies the counterfactual ensemble
tmini : int
    index of the first history file to include in the calculation
tmaxi : int
    index of the last history file to include in the calculation
overwrite : int
    whether or not to overwrite the zonal mean data if it already exists
    at the specified outpit file. 0 = False, 1 = True, 
    If False, then the script skips the interpolation for any file that 
    already exists
dry : int
    whether or not to do a "dry run" of the function, printing information
    about the execution but skipping the actual averaging. 0 = False, 1 = True
add_tropopause : bool, optional
    If this script has already been run and produced zonal-mean data outputs, 
    this argument can be set to add tropopause information from the native data
    to the datasets, rather than rerunning the procedure for all of the variables

Note that tmini, tmaxi specify only the *index* of the identified time files. 
Each file may not cover the same amount of model time, depending on the histnum. 
However, for the limvar ensembles, both the monthly (histnum=0) and daily (histnum=1) 
history files include a single month of time, and thus a given combination of 
(tmini,tmaxi) will specify the same period of simulated time for either choice of 
histnum. 
'''
import os.path
import sys
import pdb
import glob
import dask
import numpy as np
import xarray as xr
import PyTEMDiags as pt
from datetime import datetime
from combine_zinterp_native import combine_interp_native_data

outdir = '/gpfs/cldera/data/e3sm/e3sm-cldera/CLDERA-historical/limvar_src_tag/tem_processed/interp_to_pressure_zonalMeans'
interploc   = '/gpfs/cldera/data/e3sm/e3sm-cldera/CLDERA-historical/limvar_src_tag/tem_processed/interp_to_pressure_3D'
mapdir = '/ascldap/users/jphollo/repos/PyTEMDiags/maps'

Nens         = int(sys.argv[1])
histnum      = int(sys.argv[2])
massMag      = int(sys.argv[3])
tmini        = int(sys.argv[4])
tmaxi        = int(sys.argv[5])
overwrite    = bool(int(sys.argv[6]))
dry          = bool(int(sys.argv[7]))
try:
    add_tropopause = bool(int(sys.argv[8]))
except IndexError:
    add_tropopause = None
print('args: {}'.format(sys.argv))

if(histnum == 0):
    temloc = '/ascldap/users/jphollo/data/limvar/limvar_monthly'
if(histnum == 1):
    temloc = '/ascldap/users/jphollo/data/limvar/limvar_daily'

cfb = massMag == 0  # counterfactual flag
L   = 45            # spherical harmonic max order for TEM zonal averaging
zm_debug = False    # whether or not to turn on debugging output from the zonal averager

# ----------------------------------------------------------------

# get ensemble member data
print('locating data...')
if(not cfb):
    enshist = sorted(glob.glob('{}/*{}Tg*ens{}.eam*h{}*'.format(
                               interploc, massMag, Nens, histnum)))[tmini:tmaxi+1]
else:
    enshist = sorted(glob.glob('{}/*ens{}.cf.*h{}*'.format(
                               interploc, Nens, histnum)))[tmini:tmaxi+1]

# do interpolation
print('-------- working on ENS {}, times {}-{}, {} Tg --------'.format(Nens, tmini, tmaxi, massMag))
print('first file is: {}'.format(enshist[0].split('/')[-1]))
print('last file is: {}'.format(enshist[-1].split('/')[-1]))

for i in range(len(enshist)):
    
    fname = enshist[i].split('/')[-1]
    print('\n ----- working on {}'.format(fname))
    if(dry): exit(0)
    

    # find TEM file (doesn't matter which one; needed so that we can know the correct 
    # output latitudes for the zonal mean)
    temfile = sorted(glob.glob('{}/*'.format(temloc)))[0]
    temfname = temfile.split('/')[-1]
    print('TEM reference file is {}'.format(temfname))
    if(not os.path.isfile(temfile)):
        raise RuntimeError('TEM file does not exist! Aborting...')
    
    outfile = '{}/{}_zonalmeans.nc'.format(outdir, fname.split('.nc')[0])
    if(os.path.isfile(outfile) and not overwrite):
        if(add_tropopause is None):
            print('output exists; skipping...')
            continue
        else:
            print('output exists, but add_tropopause is True; '\
                  'reading existing zonal mean data...')
            zmdata = xr.open_dataset(outfile)
        
    # open TEM data, get output lat
    temdata = xr.open_dataset(temfile)
    lat_out   = temdata['lat']

    # if add_tropopause is true and overwrite is false, assume the script has already 
    # been run, and we can add the tropopause data to the existing file at 'outfile'
    if(add_tropopause and not overwrite):
       
        # get parent dataset
        pd  = xr.open_dataset(xr.open_dataset(enshist[i]).attrs['parent_dataset'])
        lat = pd['lat']
        # create zonal averager
        ZM = pt.sph_zonal_averager(lat, lat_out, L, grid_name='180x360', 
                                   save_dest=mapdir, debug=zm_debug)
        ZM.sph_compute_matrices()
        
        zmdata['TROP_P'] = ZM.sph_zonal_mean(pd['TROP_P'].T)
        zmdata['TROP_Z'] = ZM.sph_zonal_mean(pd['TROP_Z'].T)
        zmdata.attrs['history'] = '{}, tropopause data added on {}'.format(
                                       zmdata.attrs['history'], 
                                       datetime.today().strftime('%Y-%m-%d'))
        # write out
        print('writing to {}...'.format(outfile.split('/')[-1]))
        zmdata.to_netcdf(outfile)
        continue

    # else, do the zonal averaging of the native dataset
    if(not add_tropopause):
        # call helper function to combine vertically interpolated and native data
        data    = combine_interp_native_data(enshist[i])  

        # get data information
        p0        = float(data['P0'].values)
        lat       = data['lat']
        plev      = data['plev']
        time      = data['time']
        time_bnds = data['time_bnds']
        attrs     = data.attrs

        # create zonal averager
        ZM = pt.sph_zonal_averager(lat, lat_out, L, grid_name='180x360', 
                                   save_dest=mapdir, debug=zm_debug)
        ZM.sph_compute_matrices()
       
        # zonally average
        print('---- doing zonal averaging...')
        print('u...')
        zm_ua  = ZM.sph_zonal_mean(data['U'].T)
        print('v...')
        zm_va  = ZM.sph_zonal_mean(data['V'].T)
        print('t...')
        zm_ta  = ZM.sph_zonal_mean(data['T'].T)
        print('wap...')
        zm_wap = ZM.sph_zonal_mean(data['OMEGA'].T)
        print('aoa...')
        zm_aoa = ZM.sph_zonal_mean(data['AOA'].T)
        print('e90...')
        zm_e90 = ZM.sph_zonal_mean(data['E90j'].T)
        print('trop_p...')
        zm_tropp = ZM.sph_zonal_mean(data['TROP_P'].T)
        print('trop_z...')
        zm_tropz = ZM.sph_zonal_mean(data['TROP_P'].T)
        
        # if operating on monthly data, also zonally average gravity wave forcing
        if(histnum == 0):
            print('gw1...')
            zm_gw1 = ZM.sph_zonal_mean(data['BUTGWSPEC'].T)
            print('gw2...')
            zm_gw2 = ZM.sph_zonal_mean(data['UTGWORO'].T)
            print('gw3...')
            zm_gw3 = ZM.sph_zonal_mean(data['UTGWSPEC'].T)
            
            # merge all varaibles to Dataset
            data = xr.merge([zm_ua, zm_va, zm_ta, zm_wap, zm_aoa, zm_e90, 
                             zm_tropp, zm_tropz, zm_gw1, zm_gw2, zm_gw3])

        # if operating on daily data, compute tendencies
        if(histnum == 1):
            print('---- computing tendencies...')
            zm_utend = zm_ua.diff('time', label='lower') # in m/s/day
            zm_utend = zm_utend/24/60/60  # in m/s/s
            zm_utend.name = 'UTEND'
            zm_utend.attrs['long_name'] = 'U tendency'
            zm_utend.attrs['units']     = 'm/s/s'
            
            zm_aoatend = zm_aoa.diff('time', label='lower') # in day/day
            zm_aoatend = zm_aoatend/24/60/60  # in day/s
            zm_aoatend.name = 'AOATEND'
            zm_aoatend.attrs['long_name'] = 'AOA tendency'
            zm_aoatend.attrs['units']     = 'day/s'
           
            zm_e90tend = zm_e90.diff('time', label='lower') # in kg/kg/days
            zm_e90tend = zm_e90tend/24/60/60  # in kg/kg/s
            zm_e90tend.name = 'E90TEND'
            zm_e90tend.attrs['long_name'] = 'E90 tendency'
            zm_e90tend.attrs['units']     = 'kg/kg/s'
        
            # merge all varaibles to Dataset
            data = xr.merge([zm_ua, zm_va, zm_ta, zm_wap, zm_aoa, zm_e90, 
                             zm_utend, zm_aoatend, zm_e90tend])
        
        # sanity check for nans. 180*72 nans are expected for each of the tendency 
        # variables UTEND, AOATEND, and E90TEND, since the tendency at the final
        # timestep is not computed
        numnan = sum([int(np.isnan(data[dv]).sum()) for dv in data.data_vars])
        if(numnan > (180*72*3)):
            pdb.set_trace()
            raise RuntimeError('nans found in interpolated variables! Debug')
        
        # add other variables
        data['plev']      = plev
        data['lat']       = lat_out
        data['time']      = time
        data['time_bnds'] = time_bnds
        data['P0']        = p0
        data.attrs        = attrs
        data.attrs['history'] = '{}, zonally averaged to specified latitude grid on {}'.format(
                                 attrs['history'], datetime.today().strftime('%Y-%m-%d'))
        # write out
        print('writing to {}...'.format(outfile.split('/')[-1]))
        data.to_netcdf(outfile) 
