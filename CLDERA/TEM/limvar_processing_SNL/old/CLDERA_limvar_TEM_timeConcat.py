import pdb
import sys
import glob
import numpy as np
import xarray as xr
import PyTEMDiags as pt

# output location
outdir = '/ascldap/users/jphollo/data/limvar/limvar_TEM_timeConcat'
# tem results location
temloc = '/gpfs/cldera/data/e3sm/e3sm-cldera/CLDERA-historical/limvar_src_tag/tem_processed/tem_output'
# data location
interploc     = '/gpfs/cldera/data/e3sm/e3sm-cldera/CLDERA-historical/limvar_src_tag/tem_processed/interp_to_pressure_3D'
origloc = '/gpfs/cldera/data/e3sm/e3sm-cldera/CLDERA-historical/limvar_src_tag'

# cmd args
Nens  = int(sys.argv[1])
cfb   = bool(int(sys.argv[2]))
qi    = int(sys.argv[3])
dry   = bool(int(sys.argv[4]))
print('args: {}'.format(sys.argv))

# ----------------------------------------------------------------

# get ensemble data
print('locating data for ens{}, cf {}, qi {}...'.format(Nens, cfb, qi))
cfstr = ['', '.cf'][int(cfb)]
qstr  = ['','_TRACER-AOA','_TRACER-E90j'][int(qi)]
enshist = sorted(glob.glob('{}/*ens{}{}.eam.h1*TEM*L45{}.nc'.format(temloc, Nens, cfstr, qstr)))
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
