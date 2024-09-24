import pdb
import sys
import glob
import numpy as np
import xarray as xr
import PyTEMDiags as pt

# output location
outdir = '/ascldap/users/jphollo/data/limvar/limvar_zonalMeans_monthly_impact'
# tem results location
dataloc = '/ascldap/users/jphollo/data/limvar/limvar_zonalMeans_monthly_timeConcat'
# data location
interploc     = '/gpfs/cldera/data/e3sm/e3sm-cldera/CLDERA-historical/limvar_src_tag/tem_processed/interp_to_pressure_3D'
origloc = '/gpfs/cldera/data/e3sm/e3sm-cldera/CLDERA-historical/limvar_src_tag'

# cmd args
Nens  = int(sys.argv[1])
dry   = bool(int(sys.argv[2]))
print('args: {}'.format(sys.argv))

# ----------------------------------------------------------------

# get ensemble data
print('locating data for ens{}...'.format(Nens))
enshist = sorted(glob.glob('{}/*ens{}{}.eam.h0*.nc'.format(dataloc, Nens, '')))[0]
enshistcf = sorted(glob.glob('{}/*ens{}{}.eam.h0*.nc'.format(dataloc, Nens, '.cf')))[0]
fname = enshist
print('file: {}'.format(enshist))
print('cf file: {}'.format(enshistcf))
if(dry): exit(0)

print('reading data...')
enshistdat = xr.open_dataset(enshist)
enshistcfdat = xr.open_dataset(enshistcf)

print('diffing data...')
ensdiff = enshistdat - enshistcfdat
ensdiff = ensdiff.drop_vars('time_bnds')

print('writing out...')
name = fname.split('/')[-1].split('.nc')[0]
ensdiff.to_netcdf('{}/{}_impact.nc'.format(outdir, name))
