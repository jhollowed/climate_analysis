import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import glob
import pdb
import scipy
import sys


monthly_data = '/pscratch/sd/j/jhollo/E3SM/historical_data/limvar/limvar_monthly'
tembudget_data = '/pscratch/sd/j/jhollo/E3SM/historical_data/limvar/limvar_monthly_tembudget'
outdir = '/pscratch/sd/j/jhollo/E3SM/historical_data/limvar/analysis'


# -------------------------------------------------------


overwrite=0

print('reading data')

try:
    impact         = xr.load_dataset('{}/impact.nc'.format(outdir))
    dat_ensmean    = xr.load_dataset('{}/data_ensmean.nc'.format(outdir))
    cf_ensmean     = xr.load_dataset('{}/cf_ensmean.nc'.format(outdir))
    impact_ensmean = xr.load_dataset('{}/impact_ensmean.nc'.format(outdir))
    tstat          = xr.load_dataset('{}/tstat.nc'.format(outdir))
    pval           = xr.load_dataset('{}/pval.nc'.format(outdir))
    if(overwrite): raise FileNotFoundError
    print('data read from files...')

except FileNotFoundError:
    data   = sorted(glob.glob('{}/*10Tg*zonalmeans*'.format(monthly_data)))
    data   = [xr.open_dataset(d).drop_vars(['P0', 'time_bnds']) for d in data]
    cf     = sorted(glob.glob('{}/*cf*zonalmeans*'.format(monthly_data)))
    cf     = [xr.open_dataset(d).drop_vars(['P0', 'time_bnds']) for d in cf]

    print('merging data')
    data = xr.concat(data, dim='ens')
    cf   = xr.concat(cf, dim='ens')

    print('getting impact')
    impact = data - cf
    impact.to_netcdf('{}/impact.nc'.format(outdir))

    print('getting ensemble means')
    data_ensmean   = data.mean('ens') 
    cf_ensmean     = cf.mean('ens') 
    impact_ensmean = impact.mean('ens')
    data_ensmean.to_netcdf('{}/data_ensmean.nc'.format(outdir))
    cf_ensmean.to_netcdf('{}/cf_ensmean.nc'.format(outdir))
    impact_ensmean.to_netcdf('{}/impact_ensmean.nc'.format(outdir))

    print('getting ttest')
    tstat = xr.zeros_like(data_ensmean)
    pval  = xr.zeros_like(data_ensmean)
    
    for var in list(data.data_vars):
        # axis 0 is ens
        tstat_var, pval_var = scipy.stats.ttest_rel(data[var], cf[var], axis=0)
        tstat[var].values = tstat_var
        pval[var].values  = pval_var
    
    tstat.to_netcdf('{}/tstat.nc'.format(outdir))
    pval.to_netcdf('{}/pval.nc'.format(outdir))


# -------------------------------------------------------


overwrite=0

for qi in range(3):
    qstr  = ['','_TRACER-AOA','_TRACER-E90j'][int(qi)]
    print('------ working on tracer {}'.format(qi))
    print('reading tem data')
    try:
        tem_impact         = xr.load_dataset('{}/tem_impact{}.nc'.format( outdir, qstr))
        tem_dat_ensmean    = xr.load_dataset('{}/tem_data_ensmean{}.nc'.format( outdir, qstr))
        tem_cf_ensmean     = xr.load_dataset('{}/tem_cf_ensmean{}.nc'.format( outdir, qstr))
        tem_impact_ensmean = xr.load_dataset('{}/tem_impact_ensmean{}.nc'.format( outdir, qstr))
        tem_tstat          = xr.load_dataset('{}/tem_tstat{}.nc'.format( outdir, qstr))
        tem_pval           = xr.load_dataset('{}/tem_pval{}.nc'.format( outdir, qstr))
        if(overwrite): raise FileNotFoundError
        print('data read from files...')

    except FileNotFoundError:
        tem_data   = sorted(glob.glob('{}/*10Tg*TEM*L45{}_monthlymean.nc'.format(monthly_data, qstr)))
        tem_data   = [xr.open_dataset(d).drop_vars('time_bnds') for d in tem_data]
        tem_cf     = sorted(glob.glob('{}/*cf*TEM*L45{}_monthlymean.nc'.format(monthly_data, qstr)))
        tem_cf     = [xr.open_dataset(d).drop_vars('time_bnds') for d in tem_cf]

        print('merging tem data')
        tem_data = xr.concat(tem_data, dim='ens')
        tem_cf   = xr.concat(tem_cf, dim='ens')

        print('getting tem impact')
        tem_impact = tem_data - tem_cf
        tem_impact.to_netcdf('{}/tem_impact{}.nc'.format(outdir, qstr))

        print('getting tem ensemble means')
        tem_data_ensmean   = tem_data.mean('ens') 
        tem_cf_ensmean     = tem_cf.mean('ens') 
        tem_impact_ensmean = tem_impact.mean('ens')
        tem_data_ensmean.to_netcdf('{}/tem_data_ensmean{}.nc'.format(outdir, qstr))
        tem_cf_ensmean.to_netcdf('{}/tem_cf_ensmean{}.nc'.format(outdir, qstr))
        tem_impact_ensmean.to_netcdf('{}/tem_impact_ensmean{}.nc'.format(outdir, qstr))

        print('getting ttest')
        tem_tstat = xr.zeros_like(tem_data_ensmean)
        tem_pval  = xr.zeros_like(tem_data_ensmean)
        
        for var in list(tem_data.data_vars):
            # axis 0 is ens
            tstat_var, pval_var = scipy.stats.ttest_rel(tem_data[var], tem_cf[var], axis=0)
            tem_tstat[var].values = tstat_var
            tem_pval[var].values  = pval_var
        
        tem_tstat.to_netcdf('{}/tem_tstat{}.nc'.format(outdir, qstr))
        tem_pval.to_netcdf('{}/tem_pval{}.nc'.format(outdir, qstr))


# -------------------------------------------------------


overwrite=1

for qi in range(3):
    qstr  = ['','_TRACER-AOA','_TRACER-E90j'][int(qi)]
    print('------ working on tracer {}'.format(qi))
    print('reading budget data')
    try:
        budget_impact         = xr.load_dataset('{}/budget_impact{}.nc'.format( outdir, qstr))
        budget_dat_ensmean    = xr.load_dataset('{}/budget_data_ensmean{}.nc'.format( outdir, qstr))
        budget_cf_ensmean     = xr.load_dataset('{}/budget_cf_ensmean{}.nc'.format( outdir, qstr))
        budget_impact_ensmean = xr.load_dataset('{}/budget_impact_ensmean{}.nc'.format( outdir, qstr))
        budget_tstat          = xr.load_dataset('{}/budget_tstat{}.nc'.format( outdir, qstr))
        budget_pval           = xr.load_dataset('{}/budget_pval{}.nc'.format( outdir, qstr))
        if(overwrite): raise FileNotFoundError
        print('data read from files...')

    except FileNotFoundError:
        budget_data   = sorted(glob.glob('{}/*10Tg*TEMBudget{}.nc'.format(tembudget_data, qstr)))
        budget_data   = [xr.open_dataset(d) for d in budget_data]
        budget_cf     = sorted(glob.glob('{}/*cf*TEMBudget{}.nc'.format(tembudget_data, qstr)))
        budget_cf     = [xr.open_dataset(d) for d in budget_cf]

        print('merging budget data')
        budget_data = xr.concat(budget_data, dim='ens')
        budget_cf   = xr.concat(budget_cf, dim='ens')

        print('getting budget impact')
        budget_impact = budget_data - budget_cf
        budget_impact.to_netcdf('{}/budget_impact{}.nc'.format(outdir, qstr))

        print('getting budget ensemble means')
        budget_data_ensmean   = budget_data.mean('ens') 
        budget_cf_ensmean     = budget_cf.mean('ens') 
        budget_impact_ensmean = budget_impact.mean('ens')
        budget_data_ensmean.to_netcdf('{}/budget_data_ensmean{}.nc'.format(outdir, qstr))
        budget_cf_ensmean.to_netcdf('{}/budget_cf_ensmean{}.nc'.format(outdir, qstr))
        budget_impact_ensmean.to_netcdf('{}/budget_impact_ensmean{}.nc'.format(outdir, qstr))

        print('getting ttest')
        budget_tstat = xr.zeros_like(budget_data_ensmean)
        budget_pval  = xr.zeros_like(budget_data_ensmean)
        
        for var in list(budget_data.data_vars):
            # axis 0 is ens
            tstat_var, pval_var = scipy.stats.ttest_rel(budget_data[var], budget_cf[var], axis=0)
            budget_tstat[var].values = tstat_var
            budget_pval[var].values  = pval_var
        
        budget_tstat.to_netcdf('{}/budget_tstat{}.nc'.format(outdir, qstr))
        budget_pval.to_netcdf('{}/budget_pval{}.nc'.format(outdir, qstr))

# ------------------------------------------------------------------

