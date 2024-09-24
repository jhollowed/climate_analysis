import pdb
import sys
import glob
import numpy as np
import xarray as xr
import PyTEMDiags as pt

# output location
outdir = '/ascldap/users/jphollo/data/limvar/limvar_zonalMeans_monthly_ensmean'
# data location
temloc = '/ascldap/users/jphollo/data/limvar/limvar_zonalMeans_monthly_timeConcat'

# cmd args
cf    = bool(int(sys.argv[1]))
dry   = bool(int(sys.argv[2]))
print('args: {}'.format(sys.argv))

# ----------------------------------------------------------------

# get ensemble data
print('locating data...')
cfstr = ['', '.cf'][int(cf)]
enshist = sorted(glob.glob('{}/*ens[0-9]{}.eam.h0*.nc'.format(temloc, cfstr)))
print('num files: {}'.format(len(enshist)))
print('first file: {}'.format(enshist[0]))
print('last file: {}'.format(enshist[-1]))
if(dry): exit(0)

print('reading data...')
enshistdat = [xr.open_dataset(dat) for dat in enshist]
enshistdat = xr.concat(enshistdat, dim='ens')

print('averaging data...')
ensmean = enshistdat.mean('ens')
#N = len(enshistdat)
#ensmean = xr.zeros_like(enshistdat[0])
#for i in range(N):
#    print('ens {}...'.format(i+1))
#    ensmean = ensmean + enshistdat[i]
#ensmean = ensmean / N

print('finding data deviation...')
ensstd = enshistdat.std('ens')
#N = len(enshistdat)
#ensstd = xr.zeros_like(enshistdat[0])
#for i in range(N):
#    print('ens {}...'.format(i+1))
#    ensstd = ensstd + (enshistdat[i] - ensmean)**2
#ensstd = np.sqrt(ensstd / N)

print('writing out...')
name = enshist[0].split('/')[-1].split('.ens')[0]
ensmean.to_netcdf('{}/{}{}_ensmean.nc'.format(outdir, name, cfstr))
ensstd.to_netcdf('{}/{}{}_ensstd.nc'.format(outdir, name, cfstr))
