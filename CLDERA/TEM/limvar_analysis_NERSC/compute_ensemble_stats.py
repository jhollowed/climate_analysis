import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import glob
import pdb
import scipy
import sys
#local imports
import plotting_utils as putil

datadir   = '/pscratch/sd/j/jhollo/E3SM/historical_data/limvar/limvar_10daily'
outdir = '/pscratch/sd/j/jhollo/E3SM/historical_data/limvar/processed_data'


# -------------------------------------------------------

def get_10daily_stats(dataset, qi=0, overwrite=False, pmin=None, pmax=None, 
                      latmin=None, latmax=None, tmin=None, tmax=None, year=None, 
                      month=None, skip_nosrctag=False):
    '''
    Computes the ensemble mean, impact, coherence, and pvalues for the 10-daily data, with
    specified data slicing.
    The slicing uses plotting_utils.do_slicing() (see docs therein), with the same input arguments, 
    but with average=True. Meaning, the data will be averaged over the specified slice before the 
    stats are computed. 

    Parameters
    ----------
    dataset : str
        which dataset to process. Options are:
        'ens'   : the atmospheric ensemble data
        'tem'    : the TEM ensemble data
        'budget' : the TEM budget data
    qi : int
        Tracer index for the TEM and TEM budget data. Options are:
        0 : no tracer
        1 : AOA
        2:  E90
    overwrite : bool, optional
        whether or not to overwrite the data produced by a previous run of this function.
        Defaults to False, in which case the processed data is simply read and returned 
        if it already exists.
    pmin : float, optional
        minimum pressure (highest vertical level)
    pmax : float, optional
        maximum pressure (lowest vertical level)
    latmin : float, optional
        minimum latitude 
    latmax : float, optional
        maximum latitude
    tmin : cftime._cftime.Datetime, optional
        minimum time
    tmax : cftime._cftime.Datetime, optional
        maximum time
    year : int, optional
        year to slice
    month : int, optional
        month to slice
    '''
    
    # --- get tracer string
    qstr  = ['','_TRACER-AOA','_TRACER-E90j'][int(qi)]
    # turn this option off for tracers, since the non-source-tagged data doesn't exist for tracers
    if(qi != 0): skip_nosrtag=False 
    
    # --- build string to append to the end of the filenames to specify the slicing
    ttos = lambda time: str(time).split(' ')[0] # convert cftime datetime to a string, with time of day removed
    slice_args = {'pmin':pmin, 'pmax':pmax, 'latmin':latmin, 'latmax':latmax, 
                  'tmin':tmin, 'tmax':tmax, 'year':year, 'month':month, 'average':True}
    slice_strs = {'plev':['_plev{}-{}'.format(pmin,pmax), ''][pmin is None], 
                  'lat' :['_lat{}-{}'.format(latmin,latmax), ''][latmin is None], 
                  'time' :['_{}--{}'.format(ttos(tmin),ttos(tmax)), ''][tmin is None]}
    if(year is not None and month is None):
        slice_strs['time'] = '_{}'.format(year)
    if(year is not None and month is not None):
        slice_strs['time'] = '_{}-{}'.format(year,month)
    skip_nosrctag_str = ['_srctagonly', ''][skip_nosrctag is False]
    sfx = '{}{}{}{}'.format(slice_strs['plev'], slice_strs['lat'], 
                            slice_strs['time'], skip_nosrctag_str)
    #print('string added to filenames will be: {}'.format(sfx))

    if(dataset == 'ens'):
        # --- begin processing
        
        impact_str         = '{}/impact{}.nc'.format(outdir, sfx)
        ensmean_str        = '{}/data_ensmean{}.nc'.format(outdir, sfx)
        cf_ensmean_str     = '{}/cf_ensmean{}.nc'.format(outdir, sfx)
        impact_ensmean_str = '{}/impact_ensmean{}.nc'.format(outdir, sfx)
        #tstat_str          = '{}/tstat{}.nc'.format(outdir, sfx)
        pval_str           = '{}/pval{}.nc'.format(outdir, sfx)
        coherence_str      = '{}/impact_coherence{}.nc'.format(outdir, sfx)
        
        try:
            ensmean_read, cf_ensmean_read, impact_ensmean_read = 0,0,0
            impact_read, tstat_read, coherence_read = 0,0,0
            if(overwrite): 
                xr.backends.file_manager.FILE_CACHE.clear()
                raise FileNotFoundError    
            impact              = xr.open_dataset(impact_str)
            impact_read         = 1
            ensmean             = xr.open_dataset(ensmean_str)
            ensmean_read        = 1
            cf_ensmean          = xr.open_dataset(cf_ensmean_str)
            cf_ensmean_read     = 1
            impact_ensmean      = xr.open_dataset(impact_ensmean_str)
            impact_ensmean_read = 1
            #tstat               = xr.open_dataset(tstat_str)
            pval                = xr.open_dataset(pval_str)
            tstat_read          = 1
            coherence           = xr.open_dataset(coherence_str)
            coherence_read      = 1
            #print('data read from files')
            return ensmean, cf_ensmean, impact_ensmean, pval, coherence

        except FileNotFoundError:
            
            print('\n-------- processing ensemble zonal mean data...')

            # select both the source-tagged and non-source-tagged data
            data   = sorted(glob.glob('{}/*ens*[0-9].eam*zonalmeans_10daily.nc'.format(datadir)))
            cf     = sorted(glob.glob('{}/*cf*zonalmeans_10daily.nc'.format(datadir)))
            # remove the non-source-tagged data from the search results if requested
            if(skip_nosrctag):
                mask = [int(d.split('/')[-1].split('ens')[-1].split('.')[0]) < 90 for d in data]
                data = list(np.array(data)[mask])
                mask = [int(d.split('/')[-1].split('ens')[-1].split('.')[0]) < 90 for d in cf]
                cf   = list(np.array(cf)[mask])
            N = len(data)

            # remove time_bnds and P0
            for i in range(len(data)):
                data[i] = xr.open_dataset(data[i])
                if 'time_bnds' in data[i].data_vars: data[i] = data[i].drop_vars(['time_bnds'])
                if 'P0' in data[i].data_vars:        data[i] = data[i].drop_vars(['P0'])
            for i in range(len(cf)):
                cf[i] = xr.open_dataset(cf[i])
                if 'time_bnds' in cf[i].data_vars: cf[i] = cf[i].drop_vars(['time_bnds'])
                if 'P0' in cf[i].data_vars:        cf[i] = cf[i].drop_vars(['P0'])

            # do slicing
            print('doing slicing and averaging...')
            for i in range(len(data)):
                data[i] = putil.do_slicing(data[i], **slice_args)
            for i in range(len(data)):
                cf[i] = putil.do_slicing(cf[i], **slice_args)

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
                impact.to_netcdf(impact_str)

            # ---------- means
            if(ensmean_read == 0):
                print('getting data ensemble mean')
                ensmean = data.mean('ens') 
                ensmean.to_netcdf(ensmean_str)
            if(cf_ensmean_read == 0):
                print('getting cf ensemble mean')
                cf_ensmean = cf.mean('ens') 
                cf_ensmean.to_netcdf(cf_ensmean_str)
            if(impact_ensmean_read == 0):
                print('getting impact ensemble mean')
                impact_ensmean = impact.mean('ens')
                impact_ensmean.to_netcdf(impact_ensmean_str)

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
                #tstat.to_netcdf(tstat_str) # not outputting t-stat for now
                pval.to_netcdf(pval_str)

            # ---------- coherence
            if(coherence_read == 0):
                print('getting coherence')
                coherence = np.sign(impact) == np.sign(impact_ensmean)
                coherence = coherence.sum(dim='ens') / N
                coherence.to_netcdf(coherence_str)
                
            # done, return
            return ensmean, cf_ensmean, impact_ensmean, pval, coherence

    # -------------------------------------------------------

    if(dataset == 'tem'):
        
        impact_str         = '{}/tem_impact{}{}.nc'.format(outdir, qstr, sfx)
        ensmean_str        = '{}/tem_ensmean{}{}.nc'.format(outdir, qstr, sfx)
        cf_ensmean_str     = '{}/tem_cf_ensmean{}{}.nc'.format(outdir, qstr, sfx)
        impact_ensmean_str = '{}/tem_impact_ensmean{}{}.nc'.format(outdir, qstr, sfx)
        #tstat_str          = '{}/tem_tstat{}{}.nc'.format(outdir, qstr, sfx)
        pval_str           = '{}/tem_pval{}{}.nc'.format(outdir, qstr, sfx)
        coherence_str      = '{}/tem_impact_coherence{}{}.nc'.format(outdir, qstr, sfx)
        
        try:
            tem_ensmean_read, tem_cf_ensmean_read, tem_impact_ensmean_read = 0,0,0
            tem_impact_read, tem_tstat_read, tem_coherence_read = 0,0,0
            if(overwrite): 
                xr.backends.file_manager.FILE_CACHE.clear()
                raise FileNotFoundError
            tem_impact              = xr.open_dataset(impact_str)
            tem_impact_read         = 1
            tem_ensmean             = xr.open_dataset(ensmean_str)
            tem_ensmean_read        = 1
            tem_cf_ensmean          = xr.open_dataset(cf_ensmean_str)
            tem_cf_ensmean_read     = 1
            tem_impact_ensmean      = xr.open_dataset(impact_ensmean_str)
            tem_impact_ensmean_read = 1
            #tem_tstat               = xr.open_dataset(tstat_str)
            tem_pval                = xr.open_dataset(pval_str)
            tem_tstat_read          = 1
            tem_coherence           = xr.open_dataset(coherence_str)
            tem_coherence_read      = 1
            #print('data read from files')
            return tem_ensmean, tem_cf_ensmean, tem_impact_ensmean, tem_pval, tem_coherence

        except FileNotFoundError:
            
            print('\n-------- processing TEM data for tracer {}...'.format(qi))
            
            # selects only the source-tagged data for qi>0
            tem    = sorted(glob.glob('{}/*ens*[0-9].eam*TEM*L45{}_10daily.nc'.format(datadir,qstr)))
            tem_cf = sorted(glob.glob('{}/*cf*TEM*L45{}_10daily.nc'.format(datadir,qstr)))
            # remove the non-source-tagged data from the search results if requested
            if(skip_nosrctag):
                mask   = [int(d.split('/')[-1].split('ens')[-1].split('.')[0]) < 90 for d in tem]
                tem    = list(np.array(data)[mask])
                mask   = [int(d.split('/')[-1].split('ens')[-1].split('.')[0]) < 90 for d in tem_cf]
                tem_cf = list(np.array(cf)[mask])
            N = len(tem)

            # remove time_bnds and P0
            for i in range(len(tem)):
                tem[i] = xr.open_dataset(tem[i])
                if 'time_bnds' in tem[i].data_vars: tem[i] = tem[i].drop_vars(['time_bnds'])
                if 'P0' in tem[i].data_vars:       tem[i] = tem[i].drop_vars(['P0'])
            for i in range(len(tem_cf)):
                tem_cf[i] = xr.open_dataset(tem_cf[i])
                if 'time_bnds' in tem_cf[i].data_vars: tem_cf[i] = tem_cf[i].drop_vars(['time_bnds'])
                if 'P0' in tem_cf[i].data_vars:        tem_cf[i] = tem_cf[i].drop_vars(['P0'])

            # do slicing
            print('doing slicing and averaging...')
            for i in range(len(tem)):
                tem[i] = putil.do_slicing(tem[i], **slice_args)
            for i in range(len(tem_cf)):
                tem_cf[i] = putil.do_slicing(tem_cf[i], **slice_args)

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
                tem_impact.to_netcdf(impact_str)

            # ---------- means
            if(tem_ensmean_read == 0):
                print('getting tem ensemble mean')
                tem_ensmean = tem.mean('ens') 
                tem_ensmean.to_netcdf(ensmean_str)
            if(tem_cf_ensmean_read == 0):
                print('getting tem cf ensemble mean')
                tem_cf_ensmean = tem_cf.mean('ens') 
                tem_cf_ensmean.to_netcdf(cf_ensmean_str)
            if(tem_impact_ensmean_read == 0):
                print('getting tem impact ensemble mean')
                tem_impact_ensmean = tem_impact.mean('ens')
                tem_impact_ensmean.to_netcdf(impact_ensmean_str)

            # ---------- significance
            if(tem_tstat_read == 0):
                print('getting ttest')
                tem_tstat = xr.zeros_like(tem_ensmean)
                tem_pval  = xr.zeros_like(tem_ensmean) 
                for var in list(tem.data_vars):
                    # axis 0 is ens
                    tstat_var, pval_var   = scipy.stats.ttest_rel(tem[var], tem_cf[var], axis=0)
                    tem_tstat[var].values = tstat_var
                    tem_pval[var].values  = pval_var
                #tem_tstat.to_netcdf(tstat_str)
                tem_pval.to_netcdf(pval_str)

            # ---------- coherence
            if(tem_coherence_read == 0):
                print('getting coherence')
                tem_coherence = np.sign(tem_impact) == np.sign(tem_impact_ensmean)
                tem_coherence = tem_coherence.sum(dim='ens') / N
                tem_coherence.to_netcdf(coherence_str)

            # done, return
            return tem_ensmean, tem_cf_ensmean, tem_impact_ensmean, tem_pval, tem_coherence

    # -------------------------------------------------------
    
    if(dataset == 'budget'):
        
        impact_str         = '{}/budget_impact{}{}.nc'.format(outdir, qstr, sfx)
        ensmean_str        = '{}/budget_ensmean{}{}.nc'.format(outdir, qstr, sfx)
        cf_ensmean_str     = '{}/budget_cf_ensmean{}{}.nc'.format(outdir, qstr, sfx)
        impact_ensmean_str = '{}/budget_impact_ensmean{}{}.nc'.format(outdir, qstr, sfx)
        #tstat_str          = '{}/budget_tstat{}{}.nc'.format(outdir, qstr, sfx)
        pval_str           = '{}/budget_pval{}{}.nc'.format(outdir, qstr, sfx)
        coherence_str      = '{}/budget_impact_coherence{}{}.nc'.format(outdir, qstr, sfx)
        
        try:
            budget_ensmean_read, budget_cf_ensmean_read, budget_impact_ensmean_read = 0,0,0
            budget_impact_read, budget_tstat_read, budget_coherence_read = 0,0,0
            if(overwrite): 
                xr.backends.file_manager.FILE_CACHE.clear()
                raise FileNotFoundError
            budget_impact              = xr.open_dataset(impact_str)
            budget_impact_read         = 1
            budget_ensmean             = xr.open_dataset(ensmean_str)
            budget_ensmean_read        = 1
            budget_cf_ensmean          = xr.open_dataset(cf_ensmean_str)
            budget_cf_ensmean_read     = 1
            budget_impact_ensmean      = xr.open_dataset(impact_ensmean_str)
            budget_impact_ensmean_read = 1
            #budget_tstat               = xr.open_dataset(tstat_str)
            budget_pval                = xr.open_dataset(pval_str)
            budget_tstat_read          = 1
            budget_coherence           = xr.open_dataset(coherence_str)
            budget_coherence_read      = 1
            #print('data read from files')
            return budget_ensmean, budget_cf_ensmean, budget_impact_ensmean, budget_pval, budget_coherence

        except FileNotFoundError:
            
            print('\n-------- processing TEM budget data for tracer {}...'.format(qi))

            # selects only the source-tagged data for qi>0
            budget    = sorted(glob.glob('{}/*ens*[0-9].eam*TEMBudget{}_10daily.nc'.format(datadir, qstr)))
            budget    = [xr.open_dataset(d) for d in budget]
            budget_cf = sorted(glob.glob('{}/*cf*TEMBudget{}_10daily.nc'.format(datadir, qstr)))
            budget_cf = [xr.open_dataset(d) for d in budget_cf]
            # remove the non-source-tagged data from the search results if requested
            if(skip_nosrctag):
                mask      = [int(d.split('/')[-1].split('ens')[-1].split('.')[0]) < 90 for d in budget]
                budget    = list(np.array(data)[mask])
                mask      = [int(d.split('/')[-1].split('ens')[-1].split('.')[0]) < 90 for d in budget_cf]
                budget_cf = list(np.array(cf)[mask])
            N = len(budget)
            
            # do slicing
            print('doing slicing and averaging...')
            for i in range(len(budget)):
                budget[i] = putil.do_slicing(budget[i], **slice_args)
            for i in range(len(budget_cf)):
                budget_cf[i] = putil.do_slicing(budget_cf[i], **slice_args)

            print('merging budget data')
            budget    = xr.concat(budget, dim='ens')
            budget_cf = xr.concat(budget_cf, dim='ens')
            if(qi == 0):
                print('budget data shape after member concat: {}'.format(budget['UTRESVEL'].shape))
                print('budget cf shape after member concat: {}'.format(budget_cf['UTRESVEL'].shape))
            else:
                print('budget data shape after member concat: {}'.format(budget['QTRESVEL'].shape))
                print('budget cf shape after member concat: {}'.format(budget_cf['QTRESVEL'].shape))

            # ---------- impact
            if(budget_impact_read == 0):
                print('getting budget impact')
                budget_impact = budget - budget_cf
                budget_impact.to_netcdf(impact_str)

            # ---------- means
            if(budget_ensmean_read == 0):
                print('getting budget ensemble mean')
                budget_ensmean = budget.mean('ens') 
                budget_ensmean.to_netcdf(ensmean_str)
            if(budget_cf_ensmean_read == 0):
                print('getting budget cfa ensemble mean')
                budget_cf_ensmean = budget_cf.mean('ens') 
                budget_cf_ensmean.to_netcdf(cf_ensmean_str)
            if(budget_impact_ensmean_read == 0):
                print('getting budget impact ensemble mean')
                budget_impact_ensmean = budget_impact.mean('ens')
                budget_impact_ensmean.to_netcdf(impact_ensmean_str)

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
                #budget_tstat.to_netcdf(tstat_str)
                budget_pval.to_netcdf(pval_str)

            # ---------- coherence
            if(budget_coherence_read == 0):
                print('getting coherence')
                budget_coherence = np.sign(budget_impact) == np.sign(budget_impact_ensmean)
                budget_coherence = budget_coherence.sum(dim='ens') / N
                budget_coherence.to_netcdf(coherence_str)

            # done, return
            return budget_ensmean, budget_cf_ensmean, budget_impact_ensmean, budget_pval, budget_coherence
