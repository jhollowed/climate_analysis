import numpy as np
import xarray as xr
import 

DATADIR = '/pscratch/sd/j/jhollo/E3SM/historical_data/limvar/analysis'


def get_ensmean_var(var, loc=DATADIR, debug=False):
    '''
    returns a dictionary of the ensemble data, TEM data, TEM budget data, 
    as well as those dataset's counterfactuals, impacts, p-values, and 
    coherence for a specfied variable. This variable can be stored in any
    of the above mentioned datasets. They will all be searched for a variable
    with the matching name.

    Parameters
    ----------
    var : string
        the variable name
    loc : string, optional
        location of the processed ensemble mean data
    debug ;: bool, optional
        if True, print out available variables for each dataset. Defaults to False.
    '''

    print('reading data for variable {}...'.format(var))
    data      = xr.open_dataset('{}/data_ensmean.nc'.format(loc))
    cf        = xr.open_dataset('{}/cf_ensmean.nc'.format(loc))
    impact    = xr.open_dataset('{}/impact_ensmean.nc'.format(loc))
    pval      = xr.open_dataset('{}/pval.nc'.format(loc))
    coherence = xr.open_dataset('{}/impact_coherence.nc'.format(loc))
    if(debug):
        print('available data vars: {}'.format(list(data.data_...vars)))

    tem_data         = xr.open_dataset('{}/tem_data_ensmean.nc'.format(loc))
    tem_cf           = xr.open_dataset('{}/tem_cf_ensmean.nc'.format(loc))
    tem_impact       = xr.open_dataset('{}/tem_impact_ensmean.nc'.format(loc))
    tem_pval         = xr.open_dataset('{}/tem_pval.nc'.format(loc))
    tem_coherence    = xr.open_dataset('{}/tem_impact_coherence.nc'.format(loc))
    budget_data      = xr.open_dataset('{}/budget_data_ensmean.nc'.format(loc))
    budget_cf        = xr.open_dataset('{}/budget_cf_ensmean.nc'.format(loc))
    budget_impact    = xr.open_dataset('{}/budget_impact_ensmean.nc'.format(loc))
    budget_pval      = xr.open_dataset('{}/budget_pval.nc'.format(loc))
    budget_coherence = xr.open_dataset('{}/budget_impact_coherence.nc'.format(loc))
    if(debug):
        print('available TEM vars: {}'.format(list(tem_data.data_vars)))
        print('available TEM budget vars: {}'.format(list(budget_data.data_vars)))

    aoa_tem_data         = xr.open_dataset('{}/tem_data_ensmean_TRACER-AOA.nc'.format(loc))
    aoa_tem_cf           = xr.open_dataset('{}/tem_cf_ensmean_TRACER-AOA.nc'.format(loc))
    aoa_tem_impact       = xr.open_dataset('{}/tem_impact_ensmean_TRACER-AOA.nc'.format(loc))
    aoa_tem_pval         = xr.open_dataset('{}/tem_pval_TRACER-AOA.nc'.format(loc))
    aoa_tem_coherence    = xr.open_dataset('{}/tem_impact_coherence_TRACER-AOA.nc'.format(loc))
    aoa_budget_data      = xr.open_dataset('{}/budget_data_ensmean_TRACER-AOA.nc'.format(loc))
    aoa_budget_cf        = xr.open_dataset('{}/budget_cf_ensmean_TRACER-AOA.nc'.format(loc))
    aoa_budget_impact    = xr.open_dataset('{}/budget_impact_ensmean_TRACER-AOA.nc'.format(loc))
    aoa_budget_pval      = xr.open_dataset('{}/budget_pval_TRACER-AOA.nc'.format(loc))
    aoa_budget_coherence = xr.open_dataset('{}/budget_impact_coherence_TRACER-AOA.nc'.format(loc))
    if(debug):
        print('available AOA TEM vars: {}'.format(list(aoa_tem_data.data_vars)))
        print('available AOA budget vars: {}'.format(list(aoa_budget_data.data_vars)))

    e90_tem_data         = xr.open_dataset('{}/tem_data_ensmean_TRACER-E90j.nc'.format(loc))
    e90_tem_cf           = xr.open_dataset('{}/tem_cf_ensmean_TRACER-E90j.nc'.format(loc))
    e90_tem_impact       = xr.open_dataset('{}/tem_impact_ensmean_TRACER-E90j.nc'.format(loc))
    e90_tem_pval         = xr.open_dataset('{}/tem_pval_TRACER-E90j.nc'.format(loc))
    e90_tem_coherence    = xr.open_dataset('{}/tem_impact_coherence_TRACER-E90j.nc'.format(loc))
    e90_budget_data      = xr.open_dataset('{}/budget_data_ensmean_TRACER-E90j.nc'.format(loc))
    e90_budget_cf        = xr.open_dataset('{}/budget_cf_ensmean_TRACER-E90j.nc'.format(loc))
    e90_budget_impact    = xr.open_dataset('{}/budget_impact_ensmean_TRACER-E90j.nc'.format(loc))
    e90_budget_pval      = xr.open_dataset('{}/budget_pval_TRACER-E90j.nc'.format(loc))
    e90_budget_coherence = xr.open_dataset('{}/budget_impact_coherence_TRACER-E90j.nc'.format(loc))
    if(debug):
        print('available E90 TEM vars: {}'.format(list(e90_tem_data.data_vars)))
        print('available E90 budget vars: {}'.format(list(e90_budget_data.data_vars)))

    # get tropopause data
    data_tropp, cf_tropp, impact_tropp = [data['TROP_P']/100, cf['TROP_P']/100, impact['TROP_P']/100]

    # find which dataset the variable belongs to, read
    print('extracting variable...')
    try:
        data, cf, impact, pval, coherence = [data[var], cf[var], impact[var],
                                             pval[var], coherence[var]]
    except KeyError: pass
    try:
        data, cf, impact, pval, coherence = [tem_data[var], tem_cf[var], tem_impact[var],
                                             tem_pval[var], tem_coherence[var]]
    except KeyError: pass
    try:
        data, cf, impact, pval, coherence = [budget_data[var], budget_cf[var], budget_impact[var],
                                             budget_pval[var], budget_coherence[var]]
    except KeyError: pass
    try:
        if(q == 'aoa'):
            data, cf, impact, pval, coherence = [aoa_tem_data[var], aoa_tem_cf[var], aoa_tem_impact[var],
                                                 aoa_tem_pval[var], aoa_coherence[var]]
            var = '{}_aoa'.format(var)
        elif(q == 'e90'):
            data, cf, impact, pval, coherence = [e90_tem_data[var], e90_tem_cf[var], e90_tem_impact[var],
                                                 e90_tem_pval[var], e90_coherence[var]]
            var = '{}_e90'.format(var)
    except KeyError: pass
    try:
        if(q == 'aoa'):
            data, cf, impact, pval, coherence = [aoa_budget_data[var], aoa_budget_cf[var],
                                                 aoa_budget_impact[var], aoa_budget_pval[var],
                                                 aoa_coherence[var]]
            var = '{}_aoa'.format(var)
        elif(q == 'e90'):
            data, cf, impact, pval, coherence = [e90_budget_data[var], e90_budget_cf[var],
                                                 e90_budget_impact[var], e90_budget_pval[var],
                                                 e90_coherence[var]]
            var = '{}_e90'.format(var)
    except KeyError: pass

    # if variable hasn't been found by here, it doesn't exist
    assert isinstance(data, xr.core.dataarray.DataArray), 'variable {} not found!'.format(var)

    return [data, cf, impact, pval, coherence]
