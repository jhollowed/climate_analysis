import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import glob
import pdb
import scipy
import sys
#local imports
import plotting_utils as putil

# data locations
dailydir   = '/pscratch/sd/j/jhollo/E3SM/historical_data/limvar/limvar_daily'
daily10dir = '/pscratch/sd/j/jhollo/E3SM/historical_data/limvar/limvar_10daily'
monthlydir = '/pscratch/sd/j/jhollo/E3SM/historical_data/limvar/limvar_monthly'
outdir     = '/pscratch/sd/j/jhollo/E3SM/historical_data/limvar/processed_data'

# --- constants
spd = 24*60*60 # seconds per day

# --- functions

# convert cftime datetime to a string, with time of day removed
ttos = lambda time: str(time).split(' ')[0]

# remove non-source-tagged directories from a list of data directories
def remove_nosrctag(dirs):
    ens  = [int(d.split('/')[-1].split('ens')[-1].split('.')[0]) for d in dirs]
    mask = [e < 90 for e in ens]
    dirs = np.array(dirs)[mask]
    return list(dirs)

# filter data directories by mass ensemble
def filter_mass(dirs, mass):
    Tg   = [d.split('/')[-1].split('ens')[0].split('.')[-2].strip('Tg') for d in dirs]
    mask = np.zeros(len(Tg), dtype=bool)
    for i,Tgi in enumerate(Tg):
        try:
            mask[i] = float(Tgi) == mass
        except ValueError:
            mask[i] = True
    dirs = np.array(dirs)[mask]
    if(mass != 10):
        remove_nosrctag(dirs)
    return list(dirs)

# data directories
def get_data_dirs(freq, mass, qstr, skip_nosrctag=False):
    if(freq == 'daily'):
        data      = sorted(glob.glob('{}/*ens*[0-9].eam*zonalmeans.nc'.format(dailydir)))
        cf        = sorted(glob.glob('{}/*cf*zonalmeans.nc'.format(dailydir)))
        tem       = sorted(glob.glob('{}/*ens*[0-9].eam*TEM*L45{}.nc'.format(dailydir,qstr)))
        tem_cf    = sorted(glob.glob('{}/*cf*TEM*L45{}.nc'.format(dailydir,qstr)))
        budget    = sorted(glob.glob('{}/*ens*[0-9].eam*TEMBudget{}.nc'.format(dailydir,qstr)))
        budget_cf = sorted(glob.glob('{}/*cf*TEMBudget{}.nc'.format(dailydir,qstr)))
    elif(freq == '10daily'):
        data      = sorted(glob.glob('{}/*ens*[0-9].eam*zonalmeans_10daily.nc'.format(daily10dir)))
        cf        = sorted(glob.glob('{}/*cf*zonalmeans_10daily.nc'.format(daily10dir)))
        tem       = sorted(glob.glob('{}/*ens*[0-9].eam*TEM*L45{}_10daily.nc'.format(daily10dir,qstr)))
        tem_cf    = sorted(glob.glob('{}/*cf*TEM*L45{}_10daily.nc'.format(daily10dir,qstr)))
        budget    = sorted(glob.glob('{}/*ens*[0-9].eam*TEMBudget{}_10daily.nc'.format(daily10dir,qstr)))
        budget_cf = sorted(glob.glob('{}/*cf*TEMBudget{}_10daily.nc'.format(daily10dir,qstr)))
    elif(freq == 'monthly'):
        data      = sorted(glob.glob('{}/*ens*[0-9].eam*zonalmeans_monthlymean.nc'.format(monthlydir)))
        cf        = sorted(glob.glob('{}/*cf*zonalmeans_monthlymean.nc'.format(monthlydir)))
        tem       = sorted(glob.glob('{}/*ens*[0-9].eam*TEM*L45{}_monthlymean.nc'.format(monthlydir,qstr)))
        tem_cf    = sorted(glob.glob('{}/*cf*TEM*L45{}_monthlymean.nc'.format(monthlydir,qstr)))
        budget    = sorted(glob.glob('{}/*ens*[0-9].eam*monthlymean*TEMBudget{}.nc'.format(monthlydir,qstr)))
        budget_cf = sorted(glob.glob('{}/*cf*monthlymean*TEMBudget{}.nc'.format(monthlydir,qstr)))
        
    alldata = [data, cf, tem, tem_cf, budget, budget_cf]
    
    # filter by mass and source-tagging
    alldata = tuple([filter_mass(d, mass) for d in alldata])
    
    if(skip_nosrctag or mass != 10):
        return tuple([remove_nosrctag(d) for d in alldata])
    else:
        return tuple(alldata)


# ---------------------------------------------------------------------------------------------------------


def get_data_and_stats(dataset, mass, freq, qi, overwrite=False, 
                       pmin=None, pmax=None, latmin=None, latmax=None, tmin=None, tmax=None, 
                       average_pres=True, average_lat=True, average_time=True, 
                       skip_nosrctag=False, return_intersection=False, return_members=False):
    '''
    Computes the ensemble mean, impact, coherence, and pvalues for either the daily, 10-daily  
    or monthly data, with specified data slicing.
    The slicing uses plotting_utils.do_slicing() (see docs therein), with the same input 
    arguments.

    Parameters
    ----------
    dataset : str
        which dataset to process. Options are:
        'ens'    : the atmospheric ensemble data
        'tem'    : the TEM ensemble data
        'budget' : the TEM budget data
    mass : int
        the SO2 mass of the ensemble in Tg. Options are 10 or 15Tg
    freq : str
        which frequency dataset to use. Options are:
        'daily'   : the daily data
        '10daily' : the 10-daily averaged data
        'monthly' : the monthly-averaged data
    qi : int
        Tracer index for the TEM and TEM budget data. Options are:
        0 : no tracer
        1 : AOA
        2 : E90
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
    average_lat : bool, optional
        whether or not to average over the latitude slice specified by latmin,latmax
        If True, slicing is done BEFORE integrating, and otherwise, is done AFTER
        Defaults to True
    average_pres : bool, optional
        whether or not to average over the vertical slice specified by pmin,pmax
        If True, slicing is done BEFORE integrating, and otherwise, is done AFTER
        Defaults to False
    average_time : bool, optional
        whether or not to average over the temporal slice specified by tmin,tmax
        Defaults to True
    skip_nosrctag : bool, optional
        whether or not to exclude the non-source-tagged ensemble members (ens91-95) from 
        the returned data. Defaults to False.
    return_intersection : bool, optional
        whether or not to return only the portion of data where the data and the 
        counterfactual data overlap in time. Defaults to True.
    return_members : bool, optional
        whether or not to return the member-level data.
        The counterfactual, forced runs, impacts and coherence are all provided with an
        'ens' dimension
    '''

    # --- check args
    assert dataset in ['ens', 'tem', 'budget'], 'dataset must be \'ens\', \'tem\', or \'budget\', not {}'.format(dataset)
    assert freq in ['daily', '10daily', 'monthly'], 'freq must be \'daily\', \'10daily\', or \'monthly\', not {}'.format(freq)
    assert qi in [0, 1, 2], 'qi must be 0, 1, or 2, not {}'.format(qi)
    
    # --- get tracer string
    qstr  = ['','_TRACER-AOA','_TRACER-E90j'][int(qi)]
    qname = [qs.splot('-')[-1] for qs in qstr]
    # turn this option off for tracers, since the non-source-tagged data doesn't exist for tracers
    if(qi != 0): assert not skip_nosrctag, 'if qi>0, then skip_nosrctag must not be False! '\
                                           'Tracer data not defined for the non-sourced-tagged ensemble'

    # --- check slicing for monthly data
    if(freq == 'monthly' and tmin is not None):
        if(tmin.day != 1 or tmax.day != 1):
            raise RuntimeError('if freq = \'mothly\', then the time bounds tmin,tmax '\
                               'must have day=1 (the first of the month)')
    
    # --- build string to append to the end of the filenames to specify the slicing
    tavg_str = ['','Avg'][int(average_time)]
    pavg_str = ['','Avg'][int(average_pres)]
    lavg_str = ['','Avg'][int(average_lat)]
    time_slice_args = {'tmin':tmin, 'tmax':tmax, 'average':average_time}
    lat_slice_args = {'latmin':latmin, 'latmax':latmax, 'average':average_lat}
    pres_slice_args = {'pmin':pmin, 'pmax':pmax, 'average':average_pres}
    slice_strs = {'plev':['_{}plev{}-{}'.format(pavg_str,pmin,pmax), ''][pmin is None], 
                  'lat' :['_{}lat{}-{}'.format(lavg_str,latmin,latmax), ''][latmin is None], 
                  'time' :['_{}{}--{}'.format(tavg_str,ttos(tmin),ttos(tmax)), ''][tmin is None]}
    skip_nosrctag_str = ['_srctagonly', ''][skip_nosrctag is False]
    intersect_str     = ['_nointersect', ''][return_intersection is True]
    sfx = '{}{}{}_{}Tg_{}{}{}'.format(slice_strs['plev'], slice_strs['lat'], slice_strs['time'], 
                                      mass, freq, skip_nosrctag_str, intersect_str)
   
    # --- get data directories
    data, cf, tem, tem_cf, budget, budget_cf = get_data_dirs(freq, mass, qstr, skip_nosrctag)
    N= len(data)
    
    # ======================= begin processing =======================

    if(dataset == 'ens'):
        impact_str         = '{}/impact{}.nc'.format(outdir, sfx)
        ensmean_str        = '{}/data_ensmean{}.nc'.format(outdir, sfx)
        cf_ensmean_str     = '{}/cf_ensmean{}.nc'.format(outdir, sfx)
        impact_ensmean_str = '{}/impact_ensmean{}.nc'.format(outdir, sfx)
        pval_str           = '{}/pval{}.nc'.format(outdir, sfx)
        coherence_str      = '{}/coherence{}.nc'.format(outdir, sfx)
        
        means_read = False
        try:
            ensmean_read, cf_ensmean_read, impact_ensmean_read = 0,0,0
            impact_read, tstat_read, coherence_read = 0,0,0
            if(overwrite): 
                xr.backends.file_manager.FILE_CACHE.clear()
                raise FileNotFoundError    
            impact              = xr.load_dataset(impact_str)
            impact_read         = 1
            ensmean             = xr.load_dataset(ensmean_str)
            ensmean_read        = 1
            cf_ensmean          = xr.load_dataset(cf_ensmean_str)
            cf_ensmean_read     = 1
            impact_ensmean      = xr.load_dataset(impact_ensmean_str)
            impact_ensmean_read = 1
            pval                = xr.load_dataset(pval_str)
            tstat_read          = 1
            coherence           = xr.load_dataset(coherence_str)
            coherence_read      = 1
            means_read = True
            if(return_members): raise FileNotFoundError

        except FileNotFoundError:
            
            print('\n-------- processing ensemble zonal mean data...')

            # remove time_bnds and P0
            for i in range(N):
                data[i] = xr.load_dataset(data[i])
                if 'time_bnds' in data[i].data_vars: data[i] = data[i].drop_vars(['time_bnds'])
                if 'P0' in data[i].data_vars:        data[i] = data[i].drop_vars(['P0'])
                if 'UTEND' in data[i].data_vars:     data[i] = data[i].drop_vars(['UTEND'])
                if 'AOATEND' in data[i].data_vars:   data[i] = data[i].drop_vars(['AOATEND'])
                if 'E90TEND' in data[i].data_vars:   data[i] = data[i].drop_vars(['E90TEND'])
            for i in range(N):
                cf[i] = xr.load_dataset(cf[i])
                if 'time_bnds' in cf[i].data_vars: cf[i] = cf[i].drop_vars(['time_bnds'])
                if 'P0' in cf[i].data_vars:        cf[i] = cf[i].drop_vars(['P0'])
                if 'UTEND' in cf[i].data_vars:     cf[i] = cf[i].drop_vars(['UTEND'])
                if 'AOATEND' in cf[i].data_vars:   cf[i] = cf[i].drop_vars(['AOATEND'])
                if 'E90TEND' in cf[i].data_vars:   cf[i] = cf[i].drop_vars(['E90TEND'])
                
            # ensure that the data and counterfactual data have the same shapeensmean_su
            data = [data[i].transpose(*tuple(data[i].dims)) for i in range(len(data))]
            cf   = [cf[i].transpose(*tuple(data[i].dims)) for i in range(len(cf))]

            # merge ensemble members
            print('merging data')
            print('ensemble members found: {}'.format(N))
            data = xr.concat(data, dim='ens')
            cf   = xr.concat(cf, dim='ens')
            
            # check for time intersection of datasets
            # impact and significance can only be computed over regions of overlap between the 
            # data and the counterfactual. Take their intersection in time before continuing.
            # (that is, the impact, pval, and coherence returned will be defined only for time
            # coordinates where both the data and counterfactual data are defined)
            itime = np.intersect1d(data['time'].values, cf['time'].values)
            idata = data.sel(time=itime)
            icf   = cf.sel(time=itime)
            if(return_intersection):
                data = idata
                cf   = icf

            # do time and space slicing if requested
            if(tmin is not None):
                print('doing temporal slicing and averaging...')
                idata = putil.do_slicing(idata, **time_slice_args)
                icf   = putil.do_slicing(icf, **time_slice_args)
                if(return_intersection):
                    data = idata
                    cf   = icf
                else:
                    data  = putil.do_slicing(data, **time_slice_args)
                    cf    = putil.do_slicing(cf, **time_slice_args)
            if(latmin is not None):
                print('doing latitude slicing and averaging...')
                idata = putil.do_slicing(idata, **lat_slice_args)
                icf   = putil.do_slicing(icf, **lat_slice_args)
                if(return_intersection):
                    data = idata
                    cf   = icf
                else:
                    data  = putil.do_slicing(data, **lat_slice_args)
                    cf    = putil.do_slicing(cf, **lat_slice_args)
            if(pmin is not None):
                print('doing vertical slicing and averaging...')
                idata = putil.do_slicing(idata, **pres_slice_args)
                icf   = putil.do_slicing(icf, **pres_slice_args)
                if(return_intersection):
                    data = idata
                    cf   = icf
                else:
                    data  = putil.do_slicing(data, **pres_slice_args)
                    cf    = putil.do_slicing(cf, **pres_slice_args)
            print('data shape after member concat and slicing: {}'.format(data['U'].shape))
            print('cf shape after member concat and slicing: {}'.format(cf['U'].shape))

            # ---------- impact
            if(impact_read == 0):
                print('getting impact')
                impact = idata - icf
                impact.to_netcdf(impact_str)
            
            if(not means_read):
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
                    tstat = xr.zeros_like(idata.isel(ens=0, drop=True))
                    pval  = xr.zeros_like(idata.isel(ens=0, drop=True)) 
                    axis  = idata['U'].get_axis_num('ens')
                    for i, var in enumerate(list(idata.data_vars)):
                        print('{}/{}...'.format(i+1, len(idata.data_vars)), end='\r')
                        tstat_var, pval_var = scipy.stats.ttest_rel(idata[var], icf[var], axis=axis)
                        tstat[var].values   = tstat_var
                        pval[var].values    = pval_var
                    pval.to_netcdf(pval_str)

                # ---------- coherence
                if(coherence_read == 0):
                    print('getting coherence')
                    coherence = np.sign(impact) == np.sign(impact_ensmean)
                    coherence = coherence.sum(dim='ens') / N
                    coherence.to_netcdf(coherence_str)
            print('done')
                
        # done, return
        if(return_members):
            return ensmean, cf_ensmean, impact_ensmean, pval, coherence,\
                   data, cf, impact
        else:
            return ensmean, cf_ensmean, impact_ensmean, pval, coherence

    # -------------------------------------------------------

    if(dataset == 'tem'):
        impact_str         = '{}/tem_impact{}{}.nc'.format(outdir, qstr, sfx)
        ensmean_str        = '{}/tem_ensmean{}{}.nc'.format(outdir, qstr, sfx)
        cf_ensmean_str     = '{}/tem_cf_ensmean{}{}.nc'.format(outdir, qstr, sfx)
        impact_ensmean_str = '{}/tem_impact_ensmean{}{}.nc'.format(outdir, qstr, sfx)
        pval_str           = '{}/tem_pval{}{}.nc'.format(outdir, qstr, sfx)
        coherence_str      = '{}/tem_coherence{}{}.nc'.format(outdir, qstr, sfx)
        
        means_read = False
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
            tem_pval                = xr.open_dataset(pval_str)
            tem_tstat_read          = 1
            tem_coherence           = xr.open_dataset(coherence_str)
            tem_coherence_read      = 1
            means_read = True
            if(return_members): raise FileNotFoundError

        except FileNotFoundError:
            
            print('\n-------- processing TEM data for tracer {}...'.format(qi))
            xr.backends.file_manager.FILE_CACHE.clear()

            # read and remove time_bnds, P0
            for i in range(len(tem)):
                tem[i] = xr.load_dataset(tem[i])
                if 'time_bnds' in tem[i].data_vars: tem[i] = tem[i].drop_vars(['time_bnds'])
                if 'P0' in tem[i].data_vars:       tem[i] = tem[i].drop_vars(['P0'])
            for i in range(len(tem_cf)):
                tem_cf[i] = xr.load_dataset(tem_cf[i])
                if 'time_bnds' in tem_cf[i].data_vars: tem_cf[i] = tem_cf[i].drop_vars(['time_bnds'])
                if 'P0' in tem_cf[i].data_vars:        tem_cf[i] = tem_cf[i].drop_vars(['P0'])
                
            # ensure that the data and counterfactual data have the same shape
            tem = [tem[i].transpose(*tuple(tem[i].dims)) for i in range(len(tem))]
            tem_cf = [tem_cf[i].transpose(*tuple(tem[i].dims)) for i in range(len(tem_cf))]

            # merge ensemble members
            print('merging tem data')
            tem    = xr.concat(tem, dim='ens')
            tem_cf = xr.concat(tem_cf, dim='ens')
            if(qi == 0):
                print('tem data shape after member concat: {}'.format(tem['vtem'].shape))
                print('tem cf shape after member concat: {}'.format(tem_cf['vtem'].shape))
            else:
                print('tem data shape after member concat: {}'.format(tem['etfy'].shape))
                print('tem cf shape after member concat: {}'.format(tem_cf['etfy'].shape))
                
            # get data and conterfactual temporal intersection
            itime = np.intersect1d(tem['time'].values, tem_cf['time'].values)
            item    = tem.sel(time=itime)
            item_cf = tem_cf.sel(time=itime)
            if(return_intersection):
                tem    = item
                tem_cf = item_cf

            # do slicing if requested
            if(tmin is not None):
                print('doing temporal slicing and averaging...')
                item    = putil.do_slicing(item, **time_slice_args)
                item_cf = putil.do_slicing(item_cf, **time_slice_args)
                if(return_intersection):
                    tem = item
                    tem_cf = item_cf
                else:
                    tem     = putil.do_slicing(tem, **time_slice_args)
                    tem_cf  = putil.do_slicing(tem_cf, **time_slice_args)
            if(latmin is not None):
                print('doing latitude slicing and averaging...')
                item    = putil.do_slicing(item, **lat_slice_args)
                item_cf = putil.do_slicing(item_cf, **lat_slice_args)
                if(return_intersection):
                    tem = item
                    tem_cf = item_cf
                else:
                    tem     = putil.do_slicing(tem, **lat_slice_args)
                    tem_cf  = putil.do_slicing(tem_cf, **lat_slice_args)
            if(pmin is not None):
                print('doing vertical slicing and averaging...')
                item    = putil.do_slicing(item, **pres_slice_args)
                item_cf = putil.do_slicing(item_cf, **pres_slice_args)
                if(return_intersection):
                    tem = item
                    tem_cf = item_cf
                else:
                    tem     = putil.do_slicing(tem, **pres__slice_args)
                    tem_cf  = putil.do_slicing(tem_cf, **pres_slice_args)
            print('tem data shape after member concat and slicing: {}'.format(tem['epdiv'].shape))
            print('tem cf shape after member concat: and slicing {}'.format(tem_cf['epdiv'].shape))

            # ---------- impact
            if(tem_impact_read == 0):
                print('getting tem impact')
                tem_impact = item - item_cf
                tem_impact.to_netcdf(impact_str)

            if(not means_read):
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
                    tem_tstat = xr.zeros_like(item.isel(ens=0, drop=True))
                    tem_pval  = xr.zeros_like(item.isel(ens=0, drop=True)) 
                    axis = item['utendepfd'].get_axis_num('ens')
                    for i, var in enumerate(list(item.data_vars)):
                        print('{}/{}...'.format(i+1, len(item.data_vars)), end='\r')
                        tstat_var, pval_var   = scipy.stats.ttest_rel(item[var], item_cf[var], axis=axis)
                        tem_tstat[var].values = tstat_var
                        tem_pval[var].values  = pval_var
                    tem_pval.to_netcdf(pval_str)

                # ---------- coherence
                if(tem_coherence_read == 0):
                    print('getting coherence')
                    tem_coherence = np.sign(tem_impact) == np.sign(tem_impact_ensmean)
                    tem_coherence = tem_coherence.sum(dim='ens') / N
                    tem_coherence.to_netcdf(coherence_str)
            print('done')

        # done, return
        if(return_members):
            return tem_ensmean, tem_cf_ensmean, tem_impact_ensmean, tem_pval, tem_coherence,\
                   tem, tem_cf, tem_impact
        else:
            return tem_ensmean, tem_cf_ensmean, tem_impact_ensmean, tem_pval, tem_coherence

    # -------------------------------------------------------
    
    if(dataset == 'budget'):
        impact_str         = '{}/budget_impact{}{}.nc'.format(outdir, qstr, sfx)
        ensmean_str        = '{}/budget_ensmean{}{}.nc'.format(outdir, qstr, sfx)
        cf_ensmean_str     = '{}/budget_cf_ensmean{}{}.nc'.format(outdir, qstr, sfx)
        impact_ensmean_str = '{}/budget_impact_ensmean{}{}.nc'.format(outdir, qstr, sfx)
        pval_str           = '{}/budget_pval{}{}.nc'.format(outdir, qstr, sfx)
        coherence_str      = '{}/budget_coherence{}{}.nc'.format(outdir, qstr, sfx)
        
        means_read = False
        try:
            budget_ensmean_read, budget_cf_ensmean_read, budget_impact_ensmean_read = 0,0,0
            budget_impact_read, budget_tstat_read, budget_coherence_read = 0,0,0
            if(overwrite): 
                xr.backends.file_manager.FILE_CACHE.clear()
                raise FileNotFoundError
            budget_impact              = xr.load_dataset(impact_str)
            budget_impact_read         = 1
            budget_ensmean             = xr.load_dataset(ensmean_str)
            budget_ensmean_read        = 1
            budget_cf_ensmean          = xr.load_dataset(cf_ensmean_str)
            budget_cf_ensmean_read     = 1
            budget_impact_ensmean      = xr.load_dataset(impact_ensmean_str)
            budget_impact_ensmean_read = 1
            budget_pval                = xr.load_dataset(pval_str)
            budget_tstat_read          = 1
            budget_coherence           = xr.load_dataset(coherence_str)
            budget_coherence_read      = 1
            means_read = True
            if(return_members): raise FileNotFoundError

        except FileNotFoundError:
            
            print('\n-------- processing TEM budget data for tracer {}...'.format(qi))
            
            # read
            print('doing slicing and averaging...')
            for i in range(len(budget)):
                budget[i] = xr.load_dataset(budget[i])
            for i in range(len(budget_cf)):
                budget_cf[i] = xr.load_dataset(budget_cf[i])
                
            # ensure that the data and counterfactual data have the same shape
            budget    = [budget[i].transpose(*tuple(budget[i].dims)) for i in range(len(budget))]
            budget_cf = [budget_cf[i].transpose(*tuple(budget[i].dims)) for i in range(len(budget_cf))]

            # merge ensemble members
            print('merging budget data')
            budget    = xr.concat(budget, dim='ens')
            budget_cf = xr.concat(budget_cf, dim='ens')
            if(qi == 0):
                print('budget data shape after member concat: {}'.format(budget['utendresvel'].shape))
                print('budget cf shape after member concat: {}'.format(budget_cf['utendresvel'].shape))
            else:
                print('budget data shape after member concat: {}'.format(budget['qtendresvel'].shape))
                print('budget cf shape after member concat: {}'.format(budget_cf['qtendresvel'].shape))
                
            # get data and conterfactual temporal intersection
            itime = np.intersect1d(budget['time'].values, budget_cf['time'].values)
            ibudget    = budget.sel(time=itime)
            ibudget_cf = budget_cf.sel(time=itime)
            if(return_intersection):
                budget    = ibudget
                budget_cf = ibudget_cf
                
            # do slicing if requested
            if(tmin is not None):
                print('doing temporal slicing and averaging...')
                ibudget    = putil.do_slicing(ibudget, **time_slice_args)
                ibudget_cf = putil.do_slicing(ibudget_cf, **time_slice_args)
                if(return_intersection):
                    budget = ibudget
                    budget_cf = ibudget_cf
                else:
                    budget     = putil.do_slicing(budget, **time_slice_args)
                    budget_cf  = putil.do_slicing(budget_cf, **time_slice_args)
            if(latmin is not None):
                print('doing latitude slicing and averaging...')
                ibudget    = putil.do_slicing(ibudget, **lat_slice_args)
                ibudget_cf = putil.do_slicing(ibudget_cf, **lat_slice_args)
                if(return_intersection):
                    budget = ibudget
                    budget_cf = ibudget_cf
                else:
                    budget     = putil.do_slicing(budget, **lat_slice_args)
                    budget_cf  = putil.do_slicing(budget_cf, **lat_slice_args)
            if(pmin is not None):
                print('doing spatial slicing and averaging...')
                ibudget    = putil.do_slicing(ibudget, **pres_slice_args)
                ibudget_cf = putil.do_slicing(ibudget_cf, **pres_slice_args)
                if(return_intersection):
                    budget = ibudget
                    budget_cf = ibudget_cf
                else:
                    budget     = putil.do_slicing(budget, **pres_slice_args)
                    budget_cf  = putil.do_slicing(budget_cf, **pres_slice_args)
            print('budget data shape after member concat and slicing: {}'.format(budget['utendresvel'].shape))
            print('budget cf shape after member concat: and slicing {}'.format(budget_cf['utendresvel'].shape))

            # ---------- impact
            if(budget_impact_read == 0):
                print('getting budget impact')
                budget_impact = ibudget - ibudget_cf
                budget_impact.to_netcdf(impact_str)
            
            if(not means_read):
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
                    budget_tstat = xr.zeros_like(ibudget.isel(ens=0, drop=True))
                    budget_pval  = xr.zeros_like(ibudget.isel(ens=0, drop=True))
                    axis = ibudget['utendresvel'].get_axis_num('ens')
                    for i, var in enumerate(list(ibudget.data_vars)):
                        tstat_var, pval_var      = scipy.stats.ttest_rel(ibudget[var], ibudget_cf[var], axis=axis)
                        budget_tstat[var].values = tstat_var
                        budget_pval[var].values  = pval_var
                    budget_pval.to_netcdf(pval_str)

                # ---------- coherence
                if(budget_coherence_read == 0):
                    print('getting coherence')
                    budget_coherence = np.sign(budget_impact) == np.sign(budget_impact_ensmean)
                    budget_coherence = budget_coherence.sum(dim='ens') / N
                    budget_coherence.to_netcdf(coherence_str)
            print('done')

        # done, return
        if(return_members):
            return budget_ensmean, budget_cf_ensmean, budget_impact_ensmean, budget_pval, budget_coherence,\
                   budget, budget_cf, budget_impact
        else:
            return budget_ensmean, budget_cf_ensmean, budget_impact_ensmean, budget_pval, budget_coherence

        
# ---------------------------------------------------------------------------------------------------------


def get_integrated_tendencies_and_stats(tinit, mass, freq, qi, overwrite=False, skip_nosrctag=False,
                                        pmin=None, pmax=None, latmin=None, latmax=None, 
                                        tmin=None, tmax=None, average_lat=True, average_pres=False,
                                        average_time=False):
    '''
    Computes the ensemble mean, impact, coherence, and pvalues for integrated TEM tendencies
    for either the daily or 10-daily data, with specified data slicing.
    The slicing uses plotting_utils.do_slicing() (see docs therein), with the same input 
    arguments, but with average=True. Meaning, the data will be averaged over the specified 
    slice before the stats are computed. If wanting to take slices without averaging, use 
    plotting_utils.dos_slicing() after loading the data.

    Parameters
    ----------
    tinit : cftime.datetime object
        tinit : cftime.datetime object
        an initial time from which to integrate the wind and tracer time-tendencies 
        present in the data. If this time is not an exact match to a time coordinate 
        in the data, the initial condition used will be the next nearest time in the 
        dataset.
    mass : int
        the SO2 mass of the ensemble, in Tg. Options are 10 or 15
    freq : str
        which frequency dataset to use. Options are:
        'daily'   : the daily data
        '10daily' : the 10-daily averaged data
    qi : int, optional
        Tracer index for the TEM and TEM budget data. Options are:
        0 : no tracer
        1 : AOA
        2 : E90
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
    average_lat : bool, optional
        whether or not to average over the latitude slice specified by latmin,latmax
        If True, slicing is done BEFORE integrating, and otherwise, is done AFTER
        Defaults to True
    average_pres : bool, optional
        whether or not to average over the vertical slice specified by pmin,pmax
        If True, slicing is done BEFORE integrating, and otherwise, is done AFTER
        Defaults to False
    average_time : bool, optional
        whether or not to average over the temporal slice specified by tmin,tmax
        If True, slicing is done BEFORE integrating, and otherwise, is done AFTER.
        Defaults to False
    '''

    # --- check args
    assert freq in ['daily', '10daily'], 'freq must be \'daily\' or \'10daily\' not {}'.format(freq)
    assert qi in [0, 1, 2], 'qi must be 0, 1, or 2, not {}'.format(qi)
    assert tmin <= tinit < tmax, 'tinit ({}) must be between tmin ({}) and tmax ({})!'.format(tinit, tmin, tmax)

    # --- get tracer string
    qstr  = ['','_TRACER-AOA','_TRACER-E90j'][int(qi)]
    qname = [qs.splot('-')[-1] for qs in qstr]
    # turn this option off for tracers, since the non-source-tagged data doesn't exist for tracers
    if(qi != 0): assert not skip_nosrctag, 'if qi>0, then skip_nosrctag must not be False! '\
                                           'Tracer data not defined for the non-sourced-tagged ensemble'

    # --- configure integration settings
    if(qi > 0): 
        # initial condition will be a tracer mixing ratio
        initcond_type = 'q'
    else:
        # initial condition will be the zonal wind
        initcond_type = 'u'

    # --- build string to append to the end of the filenames to specify the slicing
    #slice_args = {'pmin':pmin, 'pmax':pmax, 'latmin':latmin, 'latmax':latmax, 'average':True}
    #slice_strs = {'plev':['_plev{}-{}'.format(pmin,pmax), ''][pmin is None], 
    #              'lat' :['_lat{}-{}'.format(latmin,latmax), ''][latmin is None]}
    #initTime_str = 'IC{}'.format(ttos(tinit))
    #skip_nosrctag_str = ['_srctagonly', ''][skip_nosrctag is False]
    #sfx = '{}{}{}_{}{}'.format(initTime_str, slice_strs['plev'], slice_strs['lat'], freq, skip_nosrctag_str)
    
    # --- build string to append to the end of the filenames to specify the slicing
    tavg_str = ['','Avg'][int(average_time)]
    lavg_str = ['','Avg'][int(average_lat)]
    pavg_str = ['','Avg'][int(average_pres)]
    time_slice_args = {'tmin':tmin, 'tmax':tmax, 'average':average_time}
    lat_slice_args  = {'latmin':latmin, 'latmax':latmax, 'average':average_lat}
    pres_slice_args = {'pmin':pmin, 'pmax':pmax, 'average':average_pres}
    slice_strs = {'plev':['_{}plev{}-{}'.format(pavg_str,pmin,pmax), ''][pmin is None], 
                  'lat' :['_{}lat{}-{}'.format(lavg_str,latmin,latmax), ''][latmin is None], 
                  'time' :['_{}{}--{}'.format(tavg_str,ttos(tmin),ttos(tmax)), ''][tmin is None]}
    initTime_str = 'IC{}'.format(ttos(tinit))
    skip_nosrctag_str = ['_srctagonly', ''][skip_nosrctag is False]
    sfx = '{}{}{}{}_{}Tg_{}{}'.format(initTime_str, slice_strs['plev'], slice_strs['lat'], slice_strs['time'], 
                                        mass, freq, skip_nosrctag_str)

    # --- get data directories
    data, cf, tem, tem_cf, budget, budget_cf = get_data_dirs(freq, mass, qstr, skip_nosrctag)
    N= len(data)

    # --- output files
    impact_str         = '{}/integTend_impact{}.nc'.format(outdir, sfx)
    ensmean_str        = '{}/integTend_ensmean{}.nc'.format(outdir, sfx)
    cf_ensmean_str     = '{}/integTend_cf_ensmean{}.nc'.format(outdir, sfx)
    impact_ensmean_str = '{}/integTend_impact_ensmean{}.nc'.format(outdir, sfx)
    pval_str           = '{}/integTend_pval{}.nc'.format(outdir, sfx)
    coherence_str      = '{}/integTend_coherence{}.nc'.format(outdir, sfx)

    # --- read files if exist
    try:
        ensmean_read, cf_ensmean_read, impact_ensmean_read = 0,0,0
        impact_read, tstat_read, coherence_read = 0,0,0
        if(overwrite): 
            xr.backends.file_manager.FILE_CACHE.clear()
            raise FileNotFoundError    
        impact              = xr.load_dataset(impact_str)
        impact_read         = 1
        ensmean             = xr.load_dataset(ensmean_str)
        ensmean_read        = 1
        cf_ensmean          = xr.load_dataset(cf_ensmean_str)
        cf_ensmean_read     = 1
        impact_ensmean      = xr.load_dataset(impact_ensmean_str)
        impact_ensmean_read = 1
        pval                = xr.load_dataset(pval_str)
        tstat_read          = 1
        coherence           = xr.load_dataset(coherence_str)
        coherence_read      = 1
        return ensmean, cf_ensmean, impact_ensmean, pval, coherence

    # --- else compute
    except FileNotFoundError:

        print('\n-------- integrating tendencies from initial condition at {}...'.format(ttos(tinit)))

        data_ic      = []
        cf_ic        = []

        # loop over ensemble members
        for i in range(N):
            print('working on ensemble member {}...'.format(i+1))
            
            # --- read data
            data[i]      = xr.open_dataset(data[i])
            cf[i]        = xr.open_dataset(cf[i])
            tem[i]       = xr.open_dataset(tem[i])
            tem_cf[i]    = xr.open_dataset(tem_cf[i])
            budget[i]    = xr.open_dataset(budget[i])
            budget_cf[i] = xr.open_dataset(budget_cf[i])
            
            # --- do slice if requested (should be done before tendency integartion)
            if(not average_time):
                #print('slicing in time...')
                tem[i]    = putil.do_slicing(tem[i], **time_slice_args)
                tem_cf[i] = putil.do_slicing(tem_cf[i], **time_slice_args)
                budget[i]    = putil.do_slicing(budget[i], **time_slice_args)
                budget_cf[i] = putil.do_slicing(budget_cf[i], **time_slice_args)
            if(not average_lat):
                #print('slicing in lat...')
                tem[i]    = putil.do_slicing(tem[i], **lat_slice_args)
                tem_cf[i] = putil.do_slicing(tem_cf[i], **lat_slice_args)
                budget[i]    = putil.do_slicing(budget[i], **lat_slice_args)
                budget_cf[i] = putil.do_slicing(budget_cf[i], **lat_slice_args)
            if(not average_pres):
                #print('slicing in plev...')
                tem[i]    = putil.do_slicing(tem[i], **pres_slice_args)
                tem_cf[i] = putil.do_slicing(tem_cf[i], **pres_slice_args)
                budget[i]    = putil.do_slicing(budget[i], **pres_slice_args)
                budget_cf[i] = putil.do_slicing(budget_cf[i], **pres_slice_args)
            debug=False
            if(debug):
                print('tem data shape after slicing: {}'.format(tem[i]['utendepfd'].shape))
                print('tem cf data shape after slicing: {}'.format(tem_cf[i]['utendepfd'].shape))
                print('budget data shape after slicing: {}'.format(budget[i]['utendresvel'].shape))
                print('budget cf data shape after slicing: {}'.format(budget_cf[i]['utendresvel'].shape))

            # --- get integration initial conditions
            if(initcond_type == 'u'):
                initcond_var = 'U'
                data_ic.append(data[i][initcond_var].sel(time=tinit, method='backfill'))
                cf_ic.append(cf[i][initcond_var].sel(time=tinit, method='backfill'))
            elif(initcond_type == 'q'):
                initcond_var = '{}'.format(qname)
                data_ic.append(data[i][initcond_var].sel(time=tinit, method='backfill'))
                cf_ic.append(data_cf[i][initcond_var].sel(time=tinit, method='backfill'))

            # --- integrate tendencies from TEM files
            tend_vars = list(tem[i].data_vars)
            for tv in tend_vars: 
                if '{}tend'.format(initcond_type) not in tv: 
                    tem[i] = tem[i].drop_vars(tv)
            tem[i]    = putil.integrate_tendency(tem[i], data_ic[i])
            tem_cf[i] = putil.integrate_tendency(tem_cf[i], cf_ic[i])

             # --- integrate tendencies from budget files
            tend_vars = list(budget[i].data_vars)
            for tv in tend_vars: 
                if '{}tend'.format(initcond_type) not in tv and '{}TEND'.format(initcond_var.strip('j')) not in tv: 
                        budget[i] = budget[i].drop_vars(tv)
            budget[i]    = putil.integrate_tendency(budget[i], data_ic[i])
            budget_cf[i] = putil.integrate_tendency(budget_cf[i], cf_ic[i])
            
            # --- do slice averaging if requested (must be done after tendency integration)
            if(average_time):
                print('averaging on time slice...')
                tem[i]    = putil.do_slicing(tem[i], **time_slice_args)
                tem_cf[i] = putil.do_slicing(tem_cf[i], **time_slice_args)
                budget[i]    = putil.do_slicing(budget[i], **time_slice_args)
                budget_cf[i] = putil.do_slicing(budget_cf[i], **time_slice_args)
            if(average_lat):
                print('averaging on lat slice...')
                tem[i]    = putil.do_slicing(tem[i], **lat_slice_args)
                tem_cf[i] = putil.do_slicing(tem_cf[i], **lat_slice_args)
                budget[i]    = putil.do_slicing(budget[i], **lat_slice_args)
                budget_cf[i] = putil.do_slicing(budget_cf[i], **lat_slice_args)
            if(average_pres):
                print('averaging on plev slice...')
                tem[i]    = putil.do_slicing(tem[i], **pres_slice_args)
                tem_cf[i] = putil.do_slicing(tem_cf[i], **pres_slice_args)
                budget[i]    = putil.do_slicing(budget[i], **pres_slice_args)
                budget_cf[i] = putil.do_slicing(budget_cf[i], **pres_slice_args)
            if(debug):
                print('tem data shape after slice averaging: {}'.format(tem[i]['utendepfd'].shape))
                print('cf data shape after slice averaging: {}'.format(tem_cf[i]['utendepfd'].shape))
                print('budget data shape after slice averaging: {}'.format(budget[i]['utendresvel'].shape))
                print('budget cf data shape after slice averaging: {}'.format(budget_cf[i]['utendresvel'].shape))

            # --- combine variables (reuse tem list)
            tem[i]    = xr.merge([tem[i], budget[i]])
            tem_cf[i] = xr.merge([tem_cf[i], budget_cf[i]])

        # merge, compute stats
        print('merging data...')
        print('ensemble members found: {}'.format(N))
        tem    = xr.concat(tem, dim='ens')
        tem_cf = xr.concat(tem_cf, dim='ens')

        # ---------- impact
        if(impact_read == 0):
            print('getting impact')
            impact = tem - tem_cf
            impact.to_netcdf(impact_str)

        # ---------- means
        if(ensmean_read == 0):
            print('getting data ensemble mean')
            ensmean = tem.mean('ens') 
            ensmean.to_netcdf(ensmean_str)
        if(cf_ensmean_read == 0):
            print('getting cf ensemble mean')
            cf_ensmean = tem_cf.mean('ens') 
            cf_ensmean.to_netcdf(cf_ensmean_str)
        if(impact_ensmean_read == 0):
            print('getting impact ensemble mean')
            impact_ensmean = impact.mean('ens')
            impact_ensmean.to_netcdf(impact_ensmean_str)

        # ---------- significance
        if(tstat_read == 0):
            tstat = xr.zeros_like(ensmean)
            pval  = xr.zeros_like(ensmean)
            axis = list(tem.dims).index('ens')
            for i, var in enumerate(list(tem.data_vars)):
                print('getting ttest for var {}/{}: {}...'.format(i, len(tem.data_vars), var))
                tstat_var, pval_var = scipy.stats.ttest_rel(tem[var], tem_cf[var], axis=axis)
                tstat[var].values   = tstat_var
                pval[var].values    = pval_var
            pval.to_netcdf(pval_str)

        # ---------- coherence
        if(coherence_read == 0):
            print('getting coherence')
            coherence = np.sign(impact) == np.sign(impact_ensmean)
            coherence = coherence.sum(dim='ens') / N
            coherence.to_netcdf(coherence_str)

        # done, return
        return ensmean, cf_ensmean, impact_ensmean, pval, coherence


