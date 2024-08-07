import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import glob
import pdb
import scipy
import sys


#latband_data = '/pscratch/sd/j/jhollo/E3SM/historical_data/limvar/limvar_latbands'
#tembudget_data = '/pscratch/sd/j/jhollo/E3SM/historical_data/limvar/limvar_latbands_tembudget'
#bands = ['SHpole', 'SHmid', 'tropics', 'eq', 'NHmid', 'NHpole']

latband_data = '/pscratch/sd/j/jhollo/E3SM/historical_data/limvar/limvar_latbands_10daily'
tembudget_data = '/pscratch/sd/j/jhollo/E3SM/historical_data/limvar/limvar_latbands_tembudget_10daily'
bands = ['SHpole_10daily', 'SHmid_10daily', 'tropics_10daily', 'eq_10daily', 
         'NHmid_10daily', 'NHpole_10daily']

outdir = '/pscratch/sd/j/jhollo/E3SM/historical_data/limvar/analysis'

# -------------------------------------------------------

for band in bands:

    overwrite = 0
    skip = 0
    fmt = [outdir, band]
    
    print('++++++++++ working on band {} ++++++++++'.format(band))
    try:
        ensmean_read, cf_ensmean_read, impact_ensmean_read = 0,0,0
        impact_read, tstat_read, coherence_read = 0,0,0
        if(overwrite or skip): raise FileNotFoundError    
        print('reading data')
        impact              = xr.load_dataset('{}/impact_{}.nc'.format(*fmt))
        impact_read         = 1
        ensmean             = xr.load_dataset('{}/ensmean_{}.nc'.format(*fmt))
        ensmean_read        = 1
        cf_ensmean          = xr.load_dataset('{}/cf_ensmean_{}.nc'.format(*fmt))
        cf_ensmean_read     = 1
        impact_ensmean      = xr.load_dataset('{}/impact_ensmean_{}.nc'.format(*fmt))
        impact_ensmean_read = 1
        tstat               = xr.load_dataset('{}/tstat_{}.nc'.format(*fmt))
        pval                = xr.load_dataset('{}/pval_{}.nc'.format(*fmt))
        tstat_read          = 1
        coherence           = xr.load_dataset('{}/impact_coherence_{}.nc'.format(*fmt))
        coherence_read      = 1
        print('data read from files...')

    except FileNotFoundError:
        if(not skip):
            print('reading data')
            data   = sorted(glob.glob('{}/*10Tg*zonalmeans*{}.nc'.format(latband_data, band)))
            data   = [xr.open_dataset(d).drop_vars(['P0']) for d in data]
            cf     = sorted(glob.glob('{}/*cf*zonalmeans*{}.nc'.format(latband_data, band)))
            cf     = [xr.open_dataset(d).drop_vars(['P0']) for d in cf]
            N      = len(data)

            print('merging data')
            print('ensemble members found: {}'.format(len(data)))
            data = xr.concat(data, dim='ens')
            cf   = xr.concat(cf, dim='ens')
            print('data shape after member concat: {}'.format(data['U'].shape))
            print('cf shape after member concat: {}'.format(cf['U'].shape))

            # ---------- impact
            if(impact_read == 0):
                print('getting impact')
                impact = data - cf
                impact.to_netcdf('{}/impact_{}.nc'.format(*fmt))

            # ---------- means
            if(ensmean_read == 0):
                print('getting data ensemble mean')
                ensmean   = data.mean('ens') 
                ensmean.to_netcdf('{}/ensmean_{}.nc'.format(*fmt))
            if(cf_ensmean_read == 0):
                print('getting cf ensemble mean')
                cf_ensmean     = cf.mean('ens') 
                cf_ensmean.to_netcdf('{}/cf_ensmean_{}.nc'.format(*fmt))
            if(impact_ensmean_read == 0):
                print('getting impact ensemble mean')
                impact_ensmean = impact.mean('ens')
                impact_ensmean.to_netcdf('{}/impact_ensmean_{}.nc'.format(*fmt))

            # ---------- significance
            if(tstat_read == 0):
                print('getting ttest')
                tstat = xr.zeros_like(ensmean)
                pval  = xr.zeros_like(ensmean) 
                for var in list(data.data_vars):
                    # axis 0 is ens
                    print('working on var {} ({}/{})'.format(var, 
                         list(data.data_vars).index(var)+1, len(list(data.data_vars))))
                    tstat_var, pval_var = scipy.stats.ttest_rel(data[var], cf[var], axis=0)
                    tstat[var].values = tstat_var
                    pval[var].values  = pval_var
                print('writing tstat...')
                tstat.to_netcdf('{}/tstat_{}.nc'.format(*fmt))
                print('writing pval...')
                pval.to_netcdf('{}/pval_{}.nc'.format(*fmt))
            
            # ---------- coherence
            if(coherence_read == 0):
                print('getting coherence')
                coherence = np.sign(impact) == np.sign(impact_ensmean)
                coherence = coherence.sum(dim='ens') / N
                coherence.to_netcdf('{}/impact_coherence_{}.nc'.format(*fmt))
        else: 
            print('skipping')


    # -------------------------------------------------------


    overwrite = 0
    skip = 0

    for qi in range(3):
        qstr  = ['','_TRACER-AOA','_TRACER-E90j'][int(qi)]
        fmt = [outdir, qstr, band]
        print('------ working on tracer {}'.format(qi))
        try:
            tem_ensmean_read, tem_cf_ensmean_read, tem_impact_ensmean_read = 0,0,0
            tem_impact_read, tem_tstat_read, tem_coherence_read = 0,0,0
            if(overwrite or skip): raise FileNotFoundError    
            print('reading tem data')
            tem_impact              = xr.load_dataset('{}/tem_impact{}_{}.nc'.format(*fmt))
            tem_impact_read         = 1
            tem_ensmean             = xr.load_dataset('{}/tem_ensmean{}_{}.nc'.format(*fmt))
            tem_ensmean_read        = 1
            tem_cf_ensmean          = xr.load_dataset('{}/tem_cf_ensmean{}_{}.nc'.format(*fmt))
            tem_cf_ensmean_read     = 1
            tem_impact_ensmean      = xr.load_dataset('{}/tem_impact_ensmean{}_{}.nc'.format(*fmt))
            tem_impact_ensmean_read = 1
            tem_tstat               = xr.load_dataset('{}/tem_tstat{}_{}.nc'.format(*fmt))
            tem_pval                = xr.load_dataset('{}/tem_pval{}_{}.nc'.format(*fmt))
            tem_coherence_read      = 1
            coherence               = xr.load_dataset('{}/tem_impact_coherence{}_{}.nc'.format(*fmt))
            print('data read from files...')

        except FileNotFoundError:
            if(not skip):
                print('reading tem data')
                tem_data   = sorted(glob.glob('{}/*10Tg*TEM*L45{}_latband_{}.nc'.format(
                                                              latband_data, qstr, band)))
                tem_data   = [xr.open_dataset(d) for d in tem_data]
                tem_cf     = sorted(glob.glob('{}/*cf*TEM*L45{}_latband_{}.nc'.format(
                                                            latband_data, qstr, band)))
                tem_cf     = [xr.open_dataset(d) for d in tem_cf]
                N          = len(tem_data)

                print('merging tem data')
                tem_data = xr.concat(tem_data, dim='ens')
                tem_cf   = xr.concat(tem_cf, dim='ens')
                if(qi == 0):
                    print('tem data shape after member concat: {}'.format(tem_data['vtem'].shape))
                    print('tem cf shape after member concat: {}'.format(tem_cf['vtem'].shape))
                else:
                    print('tem data shape after member concat: {}'.format(tem_data['etfy'].shape))
                    print('tem cf shape after member concat: {}'.format(tem_cf['etfy'].shape))

                # ---------- impact
                if(tem_impact_read == 0):
                    print('getting tem impact')
                    tem_impact = tem_data - tem_cf
                    tem_impact.to_netcdf('{}/tem_impact{}_{}.nc'.format(*fmt))

                # ---------- means
                if(tem_ensmean_read == 0):
                    print('getting tem ensemble mean')
                    tem_ensmean   = tem_data.mean('ens') 
                    tem_ensmean.to_netcdf('{}/tem_ensmean{}_{}.nc'.format(*fmt))
                if(tem_cf_ensmean_read == 0):
                    print('getting tem cf ensemble mean')
                    tem_cf_ensmean     = tem_cf.mean('ens') 
                    tem_cf_ensmean.to_netcdf('{}/tem_cf_ensmean{}_{}.nc'.format(*fmt))
                if(tem_impact_ensmean_read == 0):
                    print('getting tem impact ensemble mean')
                    tem_impact_ensmean = tem_impact.mean('ens')
                    tem_impact_ensmean.to_netcdf('{}/tem_impact_ensmean{}_{}.nc'.format(*fmt))

                # ---------- significance
                if(tem_tstat_read == 0):
                    print('getting ttest')
                    tem_tstat = xr.zeros_like(tem_ensmean)
                    tem_pval  = xr.zeros_like(tem_ensmean) 
                    for var in list(tem_data.data_vars):
                        # axis 0 is ens
                        print('working on var {} ({}/{})'.format(var, 
                             list(tem_data.data_vars).index(var)+1, len(list(tem_data.data_vars))))
                        tstat_var, pval_var = scipy.stats.ttest_rel(tem_data[var], tem_cf[var], axis=0)
                        tem_tstat[var].values = tstat_var
                        tem_pval[var].values  = pval_var
                    print('writing tstat...')
                    tem_tstat.to_netcdf('{}/tem_tstat{}_{}.nc'.format(*fmt))
                    print('writing pval...')
                    tem_pval.to_netcdf('{}/tem_pval{}_{}.nc'.format(*fmt))
            
                # ---------- coherence
                if(tem_coherence_read == 0):
                    print('getting coherence')
                    tem_coherence = np.sign(tem_impact) == np.sign(tem_impact_ensmean)
                    tem_coherence = tem_coherence.sum(dim='ens') / N
                    coherence.to_netcdf('{}/tem_impact_coherence{}_{}.nc'.format(*fmt))
            else: 
                print('skipping')


    # -------------------------------------------------------


    overwrite = 0
    skip = 0

    for qi in range(3):
        qstr  = ['','_TRACER-AOA','_TRACER-E90j'][int(qi)]
        fmt = [outdir, qstr, band]
        print('------ working on tracer {}'.format(qi))
        try:
            budget_ensmean_read, budget_cf_ensmean_read, budget_impact_ensmean_read = 0,0,0
            budget_impact_read, budget_tstat_read, budget_coherence_read = 0,0,0
            if(overwrite or skip): raise FileNotFoundError    
            print('reading budget data') 
            budget_impact              = xr.load_dataset('{}/budget_impact{}_{}.nc'.format(*fmt))
            budget_impact_read         = 1
            budget_ensmean             = xr.load_dataset('{}/budget_ensmean{}_{}.nc'.format(*fmt))
            budget_ensmean_read        = 1
            budget_cf_ensmean          = xr.load_dataset('{}/budget_cf_ensmean{}_{}.nc'.format(*fmt))
            budget_cf_ensmean_read     = 1
            budget_impact_ensmean      = xr.load_dataset('{}/budget_impact_ensmean{}_{}.nc'.format(*fmt))
            budget_impact_ensmean_read = 1
            budget_tstat               = xr.load_dataset('{}/budget_tstat{}_{}.nc'.format(*fmt))
            budget_pval                = xr.load_dataset('{}/budget_pval{}_{}.nc'.format(*fmt))
            budget_tstat_read          = 1
            budget_coherence      = xr.load_dataset('{}/budget_impact_coherence{}_{}.nc'.format(*fmt))
            budget_coherence_read = 1
            print('data read from files...')

        except FileNotFoundError:
            if(not skip):
                budget_data   = sorted(glob.glob('{}/*10Tg*{}*TEMBudget{}.nc'.format(
                                                      tembudget_data, band, qstr)))
                budget_data   = [xr.open_dataset(d) for d in budget_data]
                budget_cf     = sorted(glob.glob('{}/*cf*{}*TEMBudget{}.nc'.format(
                                                    tembudget_data, band, qstr)))
                budget_cf     = [xr.open_dataset(d) for d in budget_cf]
                N      = len(budget_data)

                print('merging budget data')
                budget_data = xr.concat(budget_data, dim='ens')
                budget_cf   = xr.concat(budget_cf, dim='ens')
                if(qi == 0):
                    print('budget data shape after member concat: {}'.format(
                                                budget_data['UTTOTAL_NOGW'].shape))
                    print('budget cf shape after member concat: {}'.format(
                                                budget_cf['UTTOTAL_NOGW'].shape))
                else:
                    print('budget data shape after member concat: {}'.format(
                                                budget_data['QTTOTAL'].shape))
                    print('budget cf shape after member concat: {}'.format(
                                                budget_cf['QTTOTAL'].shape))

                # ---------- impact
                if(budget_impact_read == 0):
                    print('getting budget impact')
                    budget_impact = budget_data - budget_cf
                    budget_impact.to_netcdf('{}/budget_impact{}_{}.nc'.format(*fmt))

                # ---------- means
                if(budget_ensmean_read == 0):
                    print('getting budget ensemble mean')
                    budget_ensmean   = budget_data.mean('ens') 
                    budget_ensmean.to_netcdf('{}/budget_ensmean{}_{}.nc'.format(*fmt))
                if(budget_cf_ensmean_read == 0):
                    print('getting budget cf ensemble mean')
                    budget_cf_ensmean     = budget_cf.mean('ens') 
                    budget_cf_ensmean.to_netcdf('{}/budget_cf_ensmean{}_{}.nc'.format(*fmt))
                if(budget_impact_ensmean_read == 0):
                    print('getting budget impact ensemble mean')
                    budget_impact_ensmean = budget_impact.mean('ens')
                    budget_impact_ensmean.to_netcdf('{}/budget_impact_ensmean{}_{}.nc'.format(*fmt))

                # ---------- significance
                if(budget_tstat_read == 0):
                    print('getting ttest')
                    budget_tstat = xr.zeros_like(budget_ensmean)
                    budget_pval  = xr.zeros_like(budget_ensmean) 
                    for var in list(budget_data.data_vars):
                        # axis 0 is ens
                        print('working on var {} ({}/{})'.format(var, 
                             list(budget_data.data_vars).index(var)+1, len(list(budget_data.data_vars))))
                        tstat_var, pval_var = scipy.stats.ttest_rel(budget_data[var], 
                                                                    budget_cf[var], axis=0)
                        budget_tstat[var].values = tstat_var
                        budget_pval[var].values  = pval_var
                    print('writing tstat...')
                    budget_tstat.to_netcdf('{}/budget_tstat{}_{}.nc'.format(*fmt))
                    print('writing pval...')
                    budget_pval.to_netcdf('{}/budget_pval{}_{}.nc'.format(*fmt))
            
                # ---------- coherence
                if(budget_coherence_read == 0):
                    print('getting coherence')
                    budget_coherence = np.sign(budget_impact) == np.sign(budget_impact_ensmean)
                    budget_coherence = budget_coherence.sum(dim='ens') / N
                    budget_coherence.to_netcdf('{}/budget_impact_coherence{}_{}.nc'.format(*fmt))
            else:
                print('skipping')

    # ------------------------------------------------------------------

