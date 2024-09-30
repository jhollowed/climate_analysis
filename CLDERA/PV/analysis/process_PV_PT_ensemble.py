# Joe Hollowed
# Sandia National Labs 2023
#
# README
#

import os
import pdb
import glob
import numpy as np
import xarray as xr
import climate_toolbox as ctb
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import ScalarFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
import wrappers
import climate_toolbox as ctb

# -------------------------------------------------------------------------------------

# ---------- locate data ----------
fig_dir = '/global/homes/j/jhollo/repos/climate_analysis/CLDERA/PV/analysis/historical_figs'
ens_dir = '/pscratch/sd/j/jhollo/E3SM/E3SMv2_cases/historical_cases/'
tmp_dir = '{}/post-processed'.format(ens_dir)
histnum = 1 # 1 for average daily, 2 for instant. daily

# get counterfactual and injection runs
ens_members       = sorted(glob.glob('{}/*ens[1-9]'.format(ens_dir)))
ens_native_files  = sorted(glob.glob('{}/*ens[1-9]/run/*pmcpu.limvar*eam.h{}*00.nc'.format(
                                      ens_dir, histnum)))
cf_native_files   = sorted(glob.glob('{}/*ens[1-9].cf/run/*pmcpu.limvar*eam.h{}*00.nc'.format(
                                      ens_dir, histnum)))

do_isentropes = False
lev_spec = ['_plev', '_thetalev']
sfx = lev_spec[do_isentropes]

t = len(ens_members)
concat_vars = ['PT', 'PT_TRCR', 'PV', 'PV_TRCR']
all_vars    = ['PT', 'PT_TRCR', 'PT_INCONSISTENCY', 'PV', 'PV_TRCR', 'PV_INCONSISTENCY']
tmp_meandat       = '{}/dat_native_ensmean{}.nc'.format(tmp_dir, sfx)
tmp_meandatcf     = '{}/cfdat_native_ensmean{}.nc'.format(tmp_dir, sfx)
tmp_mean_cfdiff   = '{}/dat_native_ensmean_cfdiff{}.nc'.format(tmp_dir, sfx)

# loop through ensemble, take ensemble means, interpolate to isentropes
for i in range(N):
    
    print('\n----- N={}'.format(i+1))
    tmp_dat        = '{}/dat_native_ens{}{}.nc'.format(tmp_dir, i+1, sfx)
    tmp_cfdat      = '{}/cfdat_native_ens{}{}.nc'.format(tmp_dir, i+1, sfx)
    tmp_dat_cfdiff = '{}/dat_native_ens{}_cfdiff{}.nc'.format(tmp_dir, i+1, sfx)


    # -------- concatenate
    # all data files for this ens member need to be opened and concatenated
    # (since only one month is stored on each file)
    # If this has already been done, open the result from file. Otherwise, 
    # perform the concatenation    
    force_concat = False
    try:
        if(force_concat): raise FileNotFoundError 
        dat = xr.open_dataset(tmp_dat)
        cfdat = xr.open_dataset(tmp_cfdat)
        print('read concatenated data...')
    except FileNotFoundError:
        try:
            if(not do_isentropes): raise FileNotFoundError
            dat = xr.open_dataset(tmp_dat.replace(lev_spec[1], lev_spec[0]))
            cfdat = xr.open_dataset(tmp_cfdat.replace(lev_spec[1], lev_spec[0]))
            print('read concatenated data...')
        except FileNotFoundError: 
            print('concatenating data...')
            dat_mask = ['ens{}/'.format(i+1) in f for f in ens_native_files]
            cdat_mask = ['ens{}.cf/'.format(i+1) in f for f in cf_native_files]
            this_ens_native_files = np.array(ens_native_files)[dat_mask]
            this_cf_native_files = np.array(cf_native_files)[cdat_mask]
        
            dat = xr.open_dataset(this_ens_native_files[0])
            cfdat = xr.open_dataset(this_cf_native_files[0])
            for j in range(sum(dat_mask) - 1):
                print('--- concat file {}'.format(j))
                dat = xr.concat([dat, xr.open_dataset(this_ens_native_files[j+1])],
                                dim='time', data_vars = concat_vars)
                cfdat = xr.concat([cfdat, xr.open_dataset(this_cf_native_files[j+1])], 
                                  dim='time', data_vars = concat_vars)
        dat.to_netcdf(tmp_dat)
        cfdat.to_netcdf(tmp_cfdat)
        print('wrote concatenated data...')

   
    # -------- interpolate data to isentropes
    # Conservation of PV, PT is only expected on isentropes, so let us interpolate the
    # data from pressure levels, to isentropic levels, if do_isentropes is True
    if not do_isentropes:
        pass
    else:
        force_interp = False
        try:
            if(force_interp): raise FileNotFoundError 
            _ = dat['PV']['theta']
            print('data already interpolated to isentropes...')
        except KeyError:
            print('interpolating data to isentropes...')
            pdim = dat['PV'].dims.index('lev')
            levels = np.linspace(250, 850, 61)
            for var_name in all_vars:
                print('--- interpolating {}'.format(var_name))

                # tmp...
                ctb.interp_to_isentropes(dat['PT'], dat['PT'], levels)

                dat[var_name]   = wrappers.isentrope_interp(tmp_dat, var_name, levels, pdim)
                if(var_name == 'PT_TRCR'): pdb.set_trace()
                cfdat[var_name] = wrappers.isentrope_interp(tmp_cfdat, var_name, levels, pdim)
            dat.to_netcdf(tmp_dat)
            cfdat.to_netcdf(tmp_cfdat)


    # -------- compute PV,PT inconsistency
    # we need to compute the difference between the diagnostic and tracer quantities 
    # for both PV and PT; do this and then append a new variable to the netcdf datasets 
    # created in the concatenation loop above
    # If this has already been done, skip
    force_inconsistency = False
    try:
        if(force_inconsistency): raise KeyError 
        _ = dat['PV_INCONSISTENCY']
        print('pv,pt inconsistency already computed...')
    except KeyError:
        print('computing pv,pt inconsistency...'.format(i+1))
        pv_inconsistency = dat['PV_TRCR'] - dat['PV']
        pt_inconsistency = dat['PT_TRCR'] - dat['PT']
        dat['PV_INCONSISTENCY'] = pv_inconsistency
        dat['PT_INCONSISTENCY'] = pt_inconsistency
        cf_pv_inconsistency = cfdat['PV_TRCR'] - cfdat['PV']
        cf_pt_inconsistency = cfdat['PT_TRCR'] - cfdat['PT']
        cfdat['PV_INCONSISTENCY'] = cf_pv_inconsistency
        cfdat['PT_INCONSISTENCY'] = cf_pt_inconsistency
        os.remove(tmp_dat)
        os.remove(tmp_cfdat)
        dat.to_netcdf(tmp_dat)
        cfdat.to_netcdf(tmp_cfdat)


    # -------- take forced, counterfactual differences
    # not only do we want to see the inconsistency between the diagnostic and tracer PV, 
    # but we also want to see the difference in that inconsistency between each forced run, 
    # and it's corresponding counterfactual
    force_diff = False
    try:
        if(force_diff): raise FileNotFoundError 
        dat_cfdiff = xr.open_dataset(tmp_dat_cfdiff)
        print('read pv,pt diff...')
    except FileNotFoundError:
        print('diffing ens, cf ens...'.format(i+1, i+1))
        dat_cfdiff = dat[all_vars] - cfdat[all_vars]
        dat_cfdiff.to_netcdf(tmp_dat_cfdiff)
        
            
    # -------- average
    # take ensemble average, write to file
    # If this has already been done, skip the averaging and open the result from file
    force_mean = False
    try:
        if(force_mean): raise FileNotFoundError 
        ensmean        = xr.open_dataset(tmp_meandat)
        ensmean_cf     = xr.open_dataset(tmp_meandatcf)
        ensmean_cfdiff = xr.open_dataset(tmp_mean_cfdiff)
        print('read ens mean...')
    except FileNotFoundError:
        print('summing for averaging...'.format(i))
        if(i == 0):
            ensmean        = xr.zeros_like(dat[all_vars])
            ensmean_cf     = xr.zeros_like(cfdat[all_vars])
            ensmean_cfdiff = xr.zeros_like(dat_cfdiff[all_vars])
        ensmean        = ensmean + dat[all_vars]
        ensmean_cf     = ensmean_cf + cfdat[all_vars]
        ensmean_cfdiff = ensmean_cfdiff + dat_cfdiff[all_vars]
        if(i == N-1):
            ensmean        = ensmean / N
            ensmean_cf     = ensmean_cf / N
            ensmean_cfdiff = ensmean_cfdiff / N
            ensmean.to_netcdf(tmp_meandat)
            ensmean_cf.to_netcdf(tmp_meandatcf)
            ensmean_cfdiff.to_netcdf(tmp_mean_cfdiff)
            print('wrote averaged data...')


