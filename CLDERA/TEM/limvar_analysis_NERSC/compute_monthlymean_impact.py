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

print('\n-------- working on ensemble data...')
print('reading data')
try:
    ensmean_read, cf_ensmean_read, impact_ensmean_read = 0,0,0
    impact_read, tstat_read, coherence_read = 0,0,0
    if(overwrite): raise FileNotFoundError    
    impact              = xr.load_dataset('{}/impact.nc'.format(outdir))
    impact_read         = 1
    ensmean             = xr.load_dataset('{}/data_ensmean.nc'.format(outdir))
    ensmean_read        = 1
    cf_ensmean          = xr.load_dataset('{}/cf_ensmean.nc'.format(outdir))
    cf_ensmean_read     = 1
    impact_ensmean      = xr.load_dataset('{}/impact_ensmean.nc'.format(outdir))
    impact_ensmean_read = 1
    tstat               = xr.load_dataset('{}/tstat.nc'.format(outdir))
    pval                = xr.load_dataset('{}/pval.nc'.format(outdir))
    tstat_read          = 1
    coherence           = xr.load_dataset('{}/impact_coherence.nc'.format(outdir))
    coherence_read      = 1
    print('data read from files...')

except FileNotFoundError:

    data   = sorted(glob.glob('{}/*10Tg*zonalmeans*'.format(monthly_data)))
    data   = [xr.open_dataset(d).drop_vars(['P0', 'time_bnds']) for d in data]
    cf     = sorted(glob.glob('{}/*cf*zonalmeans*'.format(monthly_data)))
    cf     = [xr.open_dataset(d).drop_vars(['P0', 'time_bnds']) for d in cf]
    N      = len(data)

    print('merging data')
    print('ensemble members found: {}'.format(N))
    data = xr.concat(data, dim='ens')
    cf   = xr.concat(cf, dim='ens')
    print('data shape after member concat: {}'.format(data['U'].shape))
    print('cf shape after member concat: {}'.format(cf['U'].shape))

    # ---------- impact
    if(impact_read == 0):
        print('getting impact')
        impact = data - cf
        impact.to_netcdf('{}/impact.nc'.format(outdir))

    # ---------- means
    if(ensmean_read == 0):
        print('getting data ensemble mean')
        ensmean = data.mean('ens') 
        ensmean.to_netcdf('{}/data_ensmean.nc'.format(outdir))
    if(cf_ensmean_read == 0):
        print('getting cf ensemble mean')
        cf_ensmean = cf.mean('ens') 
        cf_ensmean.to_netcdf('{}/cf_ensmean.nc'.format(outdir))
    if(impact_ensmean_read == 0):
        print('getting impact ensemble mean')
        impact_ensmean = impact.mean('ens')
        impact_ensmean.to_netcdf('{}/impact_ensmean.nc'.format(outdir))

    # ---------- significance
    if(tstat_read == 0):
        print('getting ttest')
        tstat = xr.zeros_like(ensmean)
        pval  = xr.zeros_like(ensmean) 
        for var in list(data.data_vars):
            # axis 0 is ens
            tstat_var, pval_var = scipy.stats.ttest_rel(data[var], cf[var], axis=0)
            tstat[var].values   = tstat_var
            pval[var].values    = pval_var
        tstat.to_netcdf('{}/tstat.nc'.format(outdir))
        pval.to_netcdf('{}/pval.nc'.format(outdir))

    # ---------- coherence
    if(coherence_read == 0):
        print('getting coherence')
        coherence = np.sign(impact) == np.sign(impact_ensmean)
        coherence = coherence.sum(dim='ens') / N
        coherence.to_netcdf('{}/impact_coherence.nc'.format(outdir))
    


# -------------------------------------------------------


overwrite=0

print('\n-------- working on TEM data...')
for qi in range(3):
    qstr  = ['','_TRACER-AOA','_TRACER-E90j'][int(qi)]
    print('------ working on tracer {}'.format(qi))
    print('reading tem data')
    try:
        tem_ensmean_read, tem_cf_ensmean_read, tem_impact_ensmean_read = 0,0,0
        tem_impact_read, tem_tstat_read, tem_coherence_read = 0,0,0
        if(overwrite): raise FileNotFoundError
        tem_impact              = xr.load_dataset('{}/tem_impact{}.nc'.format(outdir, qstr))
        tem_impact_read         = 1
        tem_ensmean             = xr.load_dataset('{}/tem_ensmean{}.nc'.format(outdir, qstr))
        tem_ensmean_read        = 1
        tem_cf_ensmean          = xr.load_dataset('{}/tem_cf_ensmean{}.nc'.format(outdir, qstr))
        tem_cf_ensmean_read     = 1
        tem_impact_ensmean      = xr.load_dataset('{}/tem_impact_ensmean{}.nc'.format(outdir, qstr))
        tem_impact_ensmean_read = 1
        tem_tstat               = xr.load_dataset('{}/tem_tstat{}.nc'.format(outdir, qstr))
        tem_pval                = xr.load_dataset('{}/tem_pval{}.nc'.format(outdir, qstr))
        tem_tstat_read          = 1
        tem_coherence           = xr.load_dataset('{}/tem_impact_coherence{}.nc'.format(outdir, qstr))
        tem_coherence_read      = 1
        print('data read from files...')

    except FileNotFoundError:
        tem    = sorted(glob.glob('{}/*10Tg*TEM*L45{}_monthlymean.nc'.format(monthly_data, qstr)))
        tem    = [xr.open_dataset(d).drop_vars('time_bnds') for d in tem]
        tem_cf = sorted(glob.glob('{}/*cf*TEM*L45{}_monthlymean.nc'.format(monthly_data, qstr)))
        tem_cf = [xr.open_dataset(d).drop_vars('time_bnds') for d in tem_cf]
        N      = len(tem)

        print('merging tem data')
        tem    = xr.concat(tem, dim='ens')
        tem_cf = xr.concat(tem_cf, dim='ens')
        if(qi == 0):
            print('tem data shape after member concat: {}'.format(tem['vtem'].shape))
            print('tem cf shape after member concat: {}'.format(tem_cf['vtem'].shape))
        else:
            print('tem data shape after member concat: {}'.format(tem['etfy'].shape))
            print('tem cf shape after member concat: {}'.format(tem_cf['etfy'].shape))

        # ---------- impact
        if(tem_impact_read == 0):
            print('getting tem impact')
            tem_impact = tem - tem_cf
            tem_impact.to_netcdf('{}/tem_impact{}.nc'.format(outdir, qstr))

        # ---------- means
        if(tem_ensmean_read == 0):
            print('getting tem ensemble mean')
            tem_ensmean = tem.mean('ens') 
            tem_ensmean.to_netcdf('{}/tem_ensmean{}.nc'.format(outdir, qstr))
        if(tem_cf_ensmean_read == 0):
            print('getting tem cf ensemble mean')
            tem_cf_ensmean = tem_cf.mean('ens') 
            tem_cf_ensmean.to_netcdf('{}/tem_cf_ensmean{}.nc'.format(outdir, qstr))
        if(tem_impact_ensmean_read == 0):
            print('getting tem impact ensemble mean')
            tem_impact_ensmean = tem_impact.mean('ens')
            tem_impact_ensmean.to_netcdf('{}/tem_impact_ensmean{}.nc'.format(outdir, qstr))

        # ---------- significance
        if(tem_tstat_read):
            print('getting ttest')
            tem_tstat = xr.zeros_like(tem_ensmean)
            tem_pval  = xr.zeros_like(tem_ensmean) 
            for var in list(tem.data_vars):
                # axis 0 is ens
                tstat_var, pval_var   = scipy.stats.ttest_rel(tem[var], tem_cf[var], axis=0)
                tem_tstat[var].values = tstat_var
                tem_pval[var].values  = pval_var
            tem_tstat.to_netcdf('{}/tem_tstat{}.nc'.format(outdir, qstr))
            tem_pval.to_netcdf('{}/tem_pval{}.nc'.format(outdir, qstr))
    
        # ---------- coherence
        if(tem_coherence_read == 0):
            print('getting coherence')
            coherence = np.sign(tem_impact) == np.sign(tem_impact_ensmean)
            coherence = coherence.sum(dim='ens') / N
            coherence.to_netcdf('{}/tem_impact_coherence.nc'.format(outdir))


# -------------------------------------------------------


overwrite=0

print('\n-------- working on TEM budget data...')
for qi in range(3):
    qstr  = ['','_TRACER-AOA','_TRACER-E90j'][int(qi)]
    print('------ working on tracer {}'.format(qi))
    print('reading budget data')
    try:
        budget_ensmean_read, budget_cf_ensmean_read, budget_impact_ensmean_read = 0,0,0
        budget_impact_read, budget_tstat_read, budget_coherence_read = 0,0,0
        if(overwrite): raise FileNotFoundError
        budget_impact              = xr.load_dataset('{}/budget_impact{}.nc'.format(outdir, qstr))
        budget_impact_read         = 1
        budget_ensmean             = xr.load_dataset('{}/budget_ensmean{}.nc'.format(outdir, qstr))
        budget_ensmean_read        = 1
        budget_cf_ensmean          = xr.load_dataset('{}/budget_cf_ensmean{}.nc'.format(outdir, qstr))
        budget_cf_ensmean_read     = 1
        budget_impact_ensmean      = xr.load_dataset('{}/budget_impact_ensmean{}.nc'.format(outdir, qstr))
        budget_impact_ensmean_read = 1
        budget_tstat               = xr.load_dataset('{}/budget_tstat{}.nc'.format(outdir, qstr))
        budget_pval                = xr.load_dataset('{}/budget_pval{}.nc'.format(outdir, qstr))
        budget_tstat_read          = 1
        budget_coherence      = xr.load_dataset('{}/budget_impact_coherence{}.nc'.format(outdir, qstr))
        budget_coherence_read = 1
        print('data read from files...')

    except FileNotFoundError:
        budget    = sorted(glob.glob('{}/*10Tg*TEMBudget{}.nc'.format(tembudget_data, qstr)))
        budget    = [xr.open_dataset(d) for d in budget]
        budget_cf = sorted(glob.glob('{}/*cf*TEMBudget{}.nc'.format(tembudget_data, qstr)))
        budget_cf = [xr.open_dataset(d) for d in budget_cf]
        N      = len(budget)

        print('merging budget data')
        budget    = xr.concat(budget, dim='ens')
        budget_cf = xr.concat(budget_cf, dim='ens')
        if(qi == 0):
            print('budget data shape after member concat: {}'.format(budget['UTTOTAL'].shape))
            print('budget cf shape after member concat: {}'.format(budget_cf['UTTOTAL'].shape))
        else:
            print('budget data shape after member concat: {}'.format(budget['QTTOTAL'].shape))
            print('budget cf shape after member concat: {}'.format(budget_cf['QTTOTAL'].shape))

        # ---------- impact
        if(budget_impact_read == 0):
            print('getting budget impact')
            budget_impact = budget - budget_cf
            budget_impact.to_netcdf('{}/budget_impact{}.nc'.format(outdir, qstr))

        # ---------- impact
        if(budget_ensmean_read == 0):
            print('getting budget ensemble mean')
            budget_ensmean = budget.mean('ens') 
            budget_ensmean.to_netcdf('{}/budget_ensmean{}.nc'.format(outdir, qstr))
        if(budget_cf_ensmean_read == 0):
            print('getting budget cfa ensemble mean')
            budget_cf_ensmean = budget_cf.mean('ens') 
            budget_cf_ensmean.to_netcdf('{}/budget_cf_ensmean{}.nc'.format(outdir, qstr))
        if(budget_impact_ensmean_read == 0):
            print('getting budget impacta ensemble mean')
            budget_impact_ensmean = budget_impact.mean('ens')
            budget_impact_ensmean.to_netcdf('{}/budget_impact_ensmean{}.nc'.format(outdir, qstr))
 
        # ---------- significance
        if(budget_tstat_read == 0):
            print('getting ttest')
            budget_tstat = xr.zeros_like(budget_ensmean)
            budget_pval  = xr.zeros_like(budget_ensmean) 
            for var in list(budget.data_vars):
                # axis 0 is ens
                tstat_var, pval_var      = scipy.stats.ttest_rel(budget[var], budget_cf[var], axis=0)
                budget_tstat[var].values = tstat_var
                budget_pval[var].values  = pval_var
            budget_tstat.to_netcdf('{}/budget_tstat{}.nc'.format(outdir, qstr))
            budget_pval.to_netcdf('{}/budget_pval{}.nc'.format(outdir, qstr))
        
        # ---------- coherence
        if(budget_coherence_read == 0):
            print('getting coherence')
            coherence = np.sign(budget_impact) == np.sign(budget_impact_ensmean)
            coherence = coherence.sum(dim='ens') / N
            coherence.to_netcdf('{}/budget_impact_coherence.nc'.format(outdir))

# ------------------------------------------------------------------

