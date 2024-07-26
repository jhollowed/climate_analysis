'''
Joe Hollowed
University of Michigan 2024

This sciprt does "extra" processing on the limvar ensemble data that was not done
during the original data-processing at Sandia, such as defining new member-level
derived quantities
'''

import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import glob
import pdb
import scipy
import sys


# -------------------------------------------------------


monthly_data = '/pscratch/sd/j/jhollo/E3SM/historical_data/limvar/limvar_monthly'
tembudget_data = '/pscratch/sd/j/jhollo/E3SM/historical_data/limvar/limvar_monthly_tembudget'


print('reading data...')

data_files = sorted(glob.glob('{}/*10Tg*zonalmeans*'.format(monthly_data)))
cf_files   = sorted(glob.glob('{}/*cf*zonalmeans*'.format(monthly_data)))

tem_files    = sorted(glob.glob('{}/*10Tg*TEM*L45_monthlymean.nc'.format(monthly_data)))
tem_cf_files = sorted(glob.glob('{}/*cf*TEM*L45_monthlymean.nc'.format(monthly_data)))
tem_budget_files = sorted(glob.glob('{}/*10Tg*_TEMBudget.nc'.format(tembudget_data)))
tem_budget_cf_files = sorted(glob.glob('{}/*cf*_TEMBudget.nc'.format(tembudget_data)))

q = '_TRACER-AOA'
aoa_tem_files    = sorted(glob.glob('{}/*10Tg*TEM*L45{}_monthlymean.nc'.format(monthly_data, q)))
aoa_tem_cf_files = sorted(glob.glob('{}/*cf*TEM*L45{}_monthlymean.nc'.format(monthly_data, q)))
aoa_tem_budget_files = sorted(glob.glob('{}/*10Tg*_TEMBudget*AOA.nc'.format(tembudget_data)))
aoa_tem_budget_cf_files = sorted(glob.glob('{}/*cf*_TEMBudget*AOA.nc'.format(tembudget_data)))

q = '_TRACER-E90j'
e90_tem_files    = sorted(glob.glob('{}/*10Tg*TEM*L45{}_monthlymean.nc'.format(monthly_data, q)))
e90_tem_cf_files = sorted(glob.glob('{}/*cf*TEM*L45{}_monthlymean.nc'.format(monthly_data, q)))
e90_tem_budget_files = sorted(glob.glob('{}/*10Tg*_TEMBudget*E90j.nc'.format(tembudget_data)))
e90_tem_budget_cf_files = sorted(glob.glob('{}/*cf*_TEMBudget*E90j.nc'.format(tembudget_data)))

# ---------------------------------------------------------------

skip = 0
if(not skip):

    print('computing total residual advection u-tendency...')
    for files in [tem_files, tem_cf_files]: 
        print('working on {} files...'.format(['counterfactual', 'data'][files[0]==tem_files[0]]))
        
        ncheck = 0
        for i in range(len(files)):
           
            # need to use load_dataset rather than open_dataset to 
            # avoid permission issue with the writeout. See:
            # https://github.com/pydata/xarray/issues/2887
            tem = xr.load_dataset(files[i])
            vtem_orig = tem['utendvtem']
            wtem_orig = tem['utendwtem']

            if('utendresid' in tem.data_vars):
                print('====== utendresid already exists in {}; skipping'.format(
                                                        files[i].split('/')[-1]))
                continue
            
            x = tem['utendvtem'] + tem['utendwtem']
            x.attrs['long name'] = 'total utend by residual velocity advection '\
                                   '(utendvtem + utendwtem)'
            x.attrs['units'] = 'm/s'
            x.attrs['history'] = 'computed in extra_processing.py on nersc'
            tem['utendresid'] = x

            # write back to original file
            tem.to_netcdf(files[i])
            
            # ensure we didn't mess up
            tem = xr.load_dataset(files[i])
            check = tem['utendvtem'] != vtem_orig
            if(check.sum() != 0):
                raise RuntimeError('corrupt values in original data file!')
            else:
                ncheck += 1
            # and that what we did was correct
            check = tem['utendresid'] != (vtem_orig + wtem_orig)
            if(check.sum() != 0):
                raise RuntimeError('calculation was not correct!')
            else:
                ncheck += 1
        print('{}/{} checks passed'.format(ncheck, len(tem_files)*2))
        print('done')


# ---------------------------------------------------------------


skip = 0
if(not skip):
    
    print('computing total residual advection aoa-tendency...')
    for files in [aoa_tem_files, aoa_tem_cf_files]:
        print('working on {} files...'.format(['counterfactual', 'data'][files[0]==aoa_tem_files[0]]))

        ncheck = 0
        for i in range(len(aoa_tem_files)):
            
            aoa_tem = xr.load_dataset(files[i])
            vtem_orig = aoa_tem['qtendvtem']
            wtem_orig = aoa_tem['qtendwtem']
            
            if('qtendresid' in aoa_tem.data_vars):
                print('====== qtendresid already exists in {}; skipping'.format(
                                                        files[i].split('/')[-1]))
                continue
            
            x = aoa_tem['qtendvtem'] + aoa_tem['qtendwtem']
            x.attrs['long name'] = 'total qtend by residual velocity advection '\
                                   '(qtendvtem + qtendwtem)'
            x.attrs['units'] = 'm/s'
            x.attrs['history'] = 'computed in extra_processing.py on nersc'
            aoa_tem['qtendresid'] = x

            # write back to original file
            aoa_tem.to_netcdf(files[i])

            # ensure we didn't mess up
            aoa_tem = xr.load_dataset(files[i])
            check = aoa_tem['qtendvtem'] != vtem_orig
            if(check.sum() != 0):
                raise RuntimeError('corrupt values in original data file!')
            else:
                ncheck += 1
            # and that what we did was correct
            check = aoa_tem['qtendresid'] != vtem_orig + wtem_orig
            if(check.sum() != 0):
                raise RuntimeError('calculation was not correct!')
            else:
                ncheck += 1

        print('{}/{} checks passed'.format(ncheck, len(tem_files)*2))
        print('done')


# ---------------------------------------------------------------


skip = 0
if(not skip):

    print('computing total residual advection e90-tendency...')
    for files in [e90_tem_files, e90_tem_cf_files]:
        print('working on {} files...'.format(['counterfactual', 'data'][files[0]==e90_tem_files[0]]))
        
        ncheck = 0
        for i in range(len(tem_files)):
            
            e90_tem = xr.load_dataset(files[i])
            vtem_orig = e90_tem['qtendvtem']
            wtem_orig = e90_tem['qtendwtem']
            
            if('qtendresid' in e90_tem.data_vars):
                print('====== qtendresid already exists in {}; skipping'.format(
                                                        files[i].split('/')[-1]))
                continue
            
            x = e90_tem['qtendvtem'] + e90_tem['qtendwtem']
            x.attrs['long name'] = 'total qtend by residual velocity advection '\
                                   '(qtendvtem + qtendwtem)'
            x.attrs['units'] = 'm/s'
            x.attrs['history'] = 'computed in extra_processing.py on nersc'
            e90_tem['qtendresid'] = x

            # write back to original file
            e90_tem.to_netcdf(files[i])

            # ensure we didn't mess up
            e90_tem = xr.load_dataset(files[i])
            check = e90_tem['qtendvtem'] != vtem_orig
            if(check.sum() != 0):
                raise RuntimeError('corrupt values in original data file!')
            else:
                ncheck += 1
            # and that what we did was correct
            check = e90_tem['qtendresid'] != vtem_orig + wtem_orig
            if(check.sum() != 0):
                raise RuntimeError('calculation was not correct!')
            else:
                ncheck += 1

        print('{}/{} checks passed'.format(ncheck, len(tem_files)*2))
        print('done')


# ---------------------------------------------------------------


skip = 0
if(not skip):

    print('computing e90 loss tendency...')
    for files in [data_files, cf_files]:
        print('working on {} files...'.format(['counterfactual', 'data'][files[0]==data_files[0]]))
        
        ncheck = 0
        for i in range(len(data_files)):
            
            data = xr.load_dataset(files[i])
            e90_orig = data['E90j']
            
            if('E90j_LOSS' in data.data_vars):
                print('====== E90j_LOSS already exists in {}; skipping'.format(
                                                        files[i].split('/')[-1]))
                continue
            
            tau = 90 * 24 * 60 * 60
            x = -1/tau * data['E90j']
            x.attrs['long name'] = 'E90 sink'
            x.attrs['units'] = 'kg/kg/s'
            x.attrs['history'] = 'computed in extra_processing.py on nersc'
            data['E90j_LOSS'] = x

            # write back to original file
            data.to_netcdf(files[i])

            # ensure we didn't mess up
            data = xr.load_dataset(files[i])
            check = data['E90j'] != e90_orig
            if(check.sum() != 0):
                raise RuntimeError('corrupt values in original data file!')
            else:
                ncheck += 1
            # and that what we did was correct
            check = data['E90j_LOSS'] != -1/tau * data['E90j']
            if(check.sum() != 0):
                raise RuntimeError('calculation was not correct!')
            else:
                ncheck += 1

        print('{}/{} checks passed'.format(ncheck, len(tem_files)*2))
        print('done')


# ---------------------------------------------------------------


skip = 0
if(not skip):

    print('computing modified e90 total tendency, tendency diff...')
    for files in [data_files, cf_files]:
        print('working on {} files...'.format(['counterfactual', 'data'][files[0]==data_files[0]]))

        budget_files = [e90_tem_budget_files, e90_tem_budget_cf_files][files[0]==data_files[0]] 
        
        ncheck = 0
        ncheck_tot = 0
        for i in range(len(data_files)):
            
            data = xr.load_dataset(files[i])
            budget_data = xr.load_dataset(budget_files[i])
            e90_orig = data['E90j']
            qttot_orig = budget_data['QTTOTAL']
            qtdiff_orig = budget_data['QTDIFF']
           
            if('QTTOTAL_WITH_LOSS' in data.data_vars):
                print('====== QTOTAL_WITH_LOSS already exists in {}; skipping'.format(
                                                        files[i].split('/')[-1]))
                continue
            
            # update QTTOTAL
            x = budget_data['QTTOTAL'] + data['E90j_LOSS']
            x.attrs = budget_data['QTTOTAL'].attrs
            x.attrs['history'] = 'computed in extra_processing.py on nersc'
            budget_data['QTTOTAL_WITH_LOSS'] = x
            
            if('QTDIFF_WITH_LOSS' in budget_data.data_vars):
                print('====== QTDIFF_WITH_LOSS already exists in {}; skipping'.format(
                                                        files[i].split('/')[-1]))
                continue
            ncheck_tot += 1
            
            # update QTDIFF
            x = data['E90TEND'] - budget_data['QTTOTAL_WITH_LOSS']
            x.attrs = budget_data['QTDIFF'].attrs
            x.attrs['history'] = 'computed in extra_processing.py on nersc'
            budget_data['QTDIFF_WITH_LOSS'] = x

            # write back to original file
            budget_data.to_netcdf(budget_files[i])

            # ensure we didn't mess up
            data = xr.load_dataset(files[i])
            check = data['E90j'] != e90_orig
            if(check.sum() != 0):
                raise RuntimeError('corrupt values in original data file!')
            else:
                ncheck += 1
            budget_data = xr.load_dataset(budget_files[i])
            check = budget_data['QTTOTAL'] != qttot_orig
            if(check.sum() != 0):
                raise RuntimeError('corrupt values in original data file!')
            else:
                ncheck += 1
            budget_data = xr.load_dataset(budget_files[i])
            check = budget_data['QTDIFF'] != qtdiff_orig
            if(check.sum() != 0):
                raise RuntimeError('corrupt values in original data file!')
            else:
                ncheck += 1

            pdb.set_trace()

        print('{}/{} checks passed'.format(ncheck, ncheck_tot))
        print('done')
