import pdb
import sys
import glob
import numpy as np
import xarray as xr
import PyTEMDiags as pt

# output location
outdir = '/ascldap/users/jphollo/data/limvar/limvar_zonalMeans_monthly_impact_ensmean'
# tem results location
dataloc = '/ascldap/users/jphollo/data/limvar/limvar_zonalMeans_monthly_impact'

# cmd args
dry   = bool(int(sys.argv[1]))
print('args: {}'.format(sys.argv))

# ----------------------------------------------------------------

# get ensemble data
print('locating data...')
enshist = sorted(glob.glob('{}/*ens*.eam.h0*_impact.nc'.format(dataloc)))
print('num files: {}'.format(len(enshist)))
print('first file: {}'.format(enshist[0]))
print('last file: {}'.format(enshist[-1]))
if(dry): exit(0)

print('reading data...')
enshistdat = [xr.open_dataset(dat) for dat in enshist]
enshistdat = xr.concat(enshistdat, dim='ens')

print('averaging data...')
ensmean = enshistdat.mean('ens')

print('finding data deviation...')
ensstd = enshistdat.std('ens')

print('writing out...')
name = enshist[0].split('/')[-1].split('.ens')[0]
ensmean.to_netcdf('{}/{}_impact_ensmean.nc'.format(outdir, name))
ensstd.to_netcdf('{}/{}_impact_ensstd.nc'.format(outdir, name))
