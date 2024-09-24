import pdb
import sys
import glob
import numpy as np
import xarray as xr
import PyTEMDiags as pt

# output location
outdir = '/ascldap/users/jphollo/data/limvar/limvar_zonalMeans_daily_timeConcat'
# zonal mean monthly data location
temloc = '/gpfs/cldera/data/e3sm/e3sm-cldera/CLDERA-historical/limvar_src_tag/tem_processed/interp_to_pressure_zonalMeans'

# cmd args
Nens  = int(sys.argv[1])
cfb   = bool(int(sys.argv[2]))
dry   = bool(int(sys.argv[3]))
print('args: {}'.format(sys.argv))

# ----------------------------------------------------------------

# get ensemble data
print('locating data for ens{}, cf {}...'.format(Nens, cfb))
cfstr = ['', '.cf'][int(cfb)]
enshist = sorted(glob.glob('{}/*ens{}{}.eam.h1*'.format(temloc, Nens, cfstr)))
print('first file: {}'.format(enshist[0]))
print('last file: {}'.format(enshist[-1]))
if(dry): exit(0)

print('reading data...')
enshistdat = [xr.open_dataset(hist) for hist in enshist]

print('combining data...')
ensconcat = xr.concat(enshistdat, dim='time')

print('writing out...')
name = enshist[0].split('/')[-1]
ensconcat.to_netcdf('{}/{}'.format(outdir, name))
