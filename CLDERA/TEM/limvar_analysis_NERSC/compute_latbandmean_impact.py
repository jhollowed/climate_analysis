import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import glob
import pdb
import scipy
import sys


latband_data = '/pscratch/sd/j/jhollo/E3SM/historical_data/limvar/limvar_latbands'
tembudget_data = '/pscratch/sd/j/jhollo/E3SM/historical_data/limvar/limvar_latbands_tembudget'
outdir = '/pscratch/sd/j/jhollo/E3SM/historical_data/limvar/analysis'

bands = ['SHpole', 'SHmid', 'tropics', 'NHmid', 'NHpole']

# -------------------------------------------------------

for band in bands:

    print('++++++++++ working on band {} ++++++++++'.format(band))

    overwrite=0
    skip = 1


    try:
        if(overwrite or skip): raise FileNotFoundError
        print('reading data')
        impact         = xr.load_dataset('{}/impact_{}.nc'.format(outdir, band))
        dat_ensmean    = xr.load_dataset('{}/data_ensmean_{}.nc'.format(outdir, band))
        cf_ensmean     = xr.load_dataset('{}/cf_ensmean_{}.nc'.format(outdir, band))
        impact_ensmean = xr.load_dataset('{}/impact_ensmean_{}.nc'.format(outdir, band))
        tstat          = xr.load_dataset('{}/tstat_{}.nc'.format(outdir, band))
        pval           = xr.load_dataset('{}/pval_{}.nc'.format(outdir, band))
        print('data read from files...')

    except FileNotFoundError:
        if(not skip):
            print('reading data')
            data   = sorted(glob.glob('{}/*10Tg*zonalmeans*{}.nc'.format(latband_data, band)))
            data   = [xr.open_dataset(d).drop_vars(['P0']) for d in data]
            cf     = sorted(glob.glob('{}/*cf*zonalmeans*{}.nc'.format(latband_data, band)))
            cf     = [xr.open_dataset(d).drop_vars(['P0']) for d in cf]

            print('merging data')
            data = xr.concat(data, dim='ens')
            cf   = xr.concat(cf, dim='ens')

            print('getting impact')
            impact = data - cf
            impact.to_netcdf('{}/impact_{}.nc'.format(outdir, band))

            print('getting ensemble means')
            data_ensmean   = data.mean('ens') 
            cf_ensmean     = cf.mean('ens') 
            impact_ensmean = impact.mean('ens')
            data_ensmean.to_netcdf('{}/data_ensmean_{}.nc'.format(outdir, band))
            cf_ensmean.to_netcdf('{}/cf_ensmean_{}.nc'.format(outdir, band))
            impact_ensmean.to_netcdf('{}/impact_ensmean_{}.nc'.format(outdir, band))

            print('getting ttest')
            tstat = xr.zeros_like(data_ensmean)
            pval  = xr.zeros_like(data_ensmean)
            
            for var in list(data.data_vars):
                # axis 0 is ens
                print('working on var {} ({}/{})'.format(var, 
                     list(data.data_vars).index(var)+1, len(list(data.data_vars))))
                tstat_var, pval_var = scipy.stats.ttest_rel(data[var], cf[var], axis=0)
                tstat[var].values = tstat_var
                pval[var].values  = pval_var
            
            print('writing tstat...')
            tstat.to_netcdf('{}/tstat_{}.nc'.format(outdir, band))
            print('writing pval...')
            pval.to_netcdf('{}/pval_{}.nc'.format(outdir, band))
        else: 
            print('skipping')


    # -------------------------------------------------------


    overwrite=0
    skip = 1

    for qi in range(3):
        qstr  = ['','_TRACER-AOA','_TRACER-E90j'][int(qi)]
        print('------ working on tracer {}'.format(qi))
        try:
            if(overwrite or skip): raise FileNotFoundError
            print('reading tem data')
            tem_impact         = xr.load_dataset('{}/tem_impact{}_{}.nc'.format(
                                                             outdir, qstr, band))
            tem_dat_ensmean    = xr.load_dataset('{}/tem_data_ensmean{}_{}.nc'.format(
                                                                   outdir, qstr, band))
            tem_cf_ensmean     = xr.load_dataset('{}/tem_cf_ensmean{}_{}.nc'.format(
                                                                 outdir, qstr, band))
            tem_impact_ensmean = xr.load_dataset('{}/tem_impact_ensmean{}_{}.nc'.format(
                                                                     outdir, qstr, band))
            tem_tstat          = xr.load_dataset('{}/tem_tstat{}_{}.nc'.format(
                                                            outdir, qstr, band))
            tem_pval           = xr.load_dataset('{}/tem_pval{}_{}.nc'.format(
                                                           outdir, qstr, band))
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

                print('merging tem data')
                tem_data = xr.concat(tem_data, dim='ens')
                tem_cf   = xr.concat(tem_cf, dim='ens')

                print('getting tem impact')
                tem_impact = tem_data - tem_cf
                tem_impact.to_netcdf('{}/tem_impact{}_{}.nc'.format(outdir, qstr, band))

                print('getting tem ensemble means')
                tem_data_ensmean   = tem_data.mean('ens') 
                tem_cf_ensmean     = tem_cf.mean('ens') 
                tem_impact_ensmean = tem_impact.mean('ens')
                tem_data_ensmean.to_netcdf('{}/tem_data_ensmean{}_{}.nc'.format(outdir, qstr, band))
                tem_cf_ensmean.to_netcdf('{}/tem_cf_ensmean{}_{}.nc'.format(outdir, qstr, band))
                tem_impact_ensmean.to_netcdf('{}/tem_impact_ensmean{}_{}.nc'.format(outdir, qstr, band))

                print('getting ttest')
                tem_tstat = xr.zeros_like(tem_data_ensmean)
                tem_pval  = xr.zeros_like(tem_data_ensmean)
                
                for var in list(tem_data.data_vars):
                    # axis 0 is ens
                    print('working on var {} ({}/{})'.format(var, 
                         list(tem_data.data_vars).index(var)+1, len(list(tem_data.data_vars))))
                    tstat_var, pval_var = scipy.stats.ttest_rel(tem_data[var], tem_cf[var], axis=0)
                    tem_tstat[var].values = tstat_var
                    tem_pval[var].values  = pval_var
                
                print('writing tstat...')
                tem_tstat.to_netcdf('{}/tem_tstat{}_{}.nc'.format(outdir, qstr, band))
                print('writing pval...')
                tem_pval.to_netcdf('{}/tem_pval{}_{}.nc'.format(outdir, qstr, band))
            else: 
                print('skipping')


    # -------------------------------------------------------


    overwrite = 1
    skip = 0

    for qi in range(3):
        qstr  = ['','_TRACER-AOA','_TRACER-E90j'][int(qi)]
        print('------ working on tracer {}'.format(qi))
        print('reading budget data')
        try:
            if(overwrite or skip): raise FileNotFoundError
            budget_impact         = xr.load_dataset('{}/budget_impact{}_{}.nc'.format(
                                                                outdir, qstr, band))
            budget_dat_ensmean    = xr.load_dataset('{}/budget_data_ensmean{}_{}.nc'.format(
                                                                      outdir, qstr, band))
            budget_cf_ensmean     = xr.load_dataset('{}/budget_cf_ensmean{}_{}.nc'.format(
                                                                    outdir, qstr, band))
            budget_impact_ensmean = xr.load_dataset('{}/budget_impact_ensmean{}_{}.nc'.format(
                                                                        outdir, qstr, band))
            budget_tstat          = xr.load_dataset('{}/budget_tstat{}_{}.nc'.format(
                                                               outdir, qstr, band))
            budget_pval           = xr.load_dataset('{}/budget_pval{}_{}.nc'.format(
                                                              outdir, qstr, band))
            print('data read from files...')

        except FileNotFoundError:
            if(not skip):
                budget_data   = sorted(glob.glob('{}/*10Tg*{}*TEMBudget{}.nc'.format(
                                                      tembudget_data, band, qstr)))
                budget_data   = [xr.open_dataset(d) for d in budget_data]
                budget_cf     = sorted(glob.glob('{}/*cf*{}*TEMBudget{}.nc'.format(
                                                    tembudget_data, band, qstr)))
                budget_cf     = [xr.open_dataset(d) for d in budget_cf]

                print('merging budget data')
                budget_data = xr.concat(budget_data, dim='ens')
                budget_cf   = xr.concat(budget_cf, dim='ens')

                print('getting budget impact')
                budget_impact = budget_data - budget_cf
                budget_impact.to_netcdf('{}/budget_impact{}_{}.nc'.format(outdir, qstr, band))

                print('getting budget ensemble means')
                budget_data_ensmean   = budget_data.mean('ens') 
                budget_cf_ensmean     = budget_cf.mean('ens') 
                budget_impact_ensmean = budget_impact.mean('ens')
                budget_data_ensmean.to_netcdf('{}/budget_data_ensmean{}_{}.nc'.format(
                                                                   outdir, qstr, band))
                budget_cf_ensmean.to_netcdf('{}/budget_cf_ensmean{}_{}.nc'.format(
                                                               outdir, qstr, band))
                budget_impact_ensmean.to_netcdf('{}/budget_impact_ensmean{}_{}.nc'.format(
                                                                       outdir, qstr, band))

                print('getting ttest')
                budget_tstat = xr.zeros_like(budget_data_ensmean)
                budget_pval  = xr.zeros_like(budget_data_ensmean)
                
                for var in list(budget_data.data_vars):
                    # axis 0 is ens
                    print('working on var {} ({}/{})'.format(var, 
                         list(budget_data.data_vars).index(var)+1, len(list(budget_data.data_vars))))
                    tstat_var, pval_var = scipy.stats.ttest_rel(budget_data[var], budget_cf[var], axis=0)
                    budget_tstat[var].values = tstat_var
                    budget_pval[var].values  = pval_var
               
                print('writing tstat...')
                budget_tstat.to_netcdf('{}/budget_tstat{}_{}.nc'.format(outdir, qstr, band))
                print('writing pval...')
                budget_pval.to_netcdf('{}/budget_pval{}_{}.nc'.format(outdir, qstr, band))
            else:
                print('skipping')

    # ------------------------------------------------------------------

