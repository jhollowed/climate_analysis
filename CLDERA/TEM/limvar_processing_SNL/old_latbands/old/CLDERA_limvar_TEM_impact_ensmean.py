import pdb
import sys
import glob
import numpy as np
import xarray as xr
import PyTEMDiags as pt

# output location
outdir = '/ascldap/users/jphollo/data/limvar/limvar_TEM_impact_ensmean'
# tem results location
temloc = '/ascldap/users/jphollo/data/limvar/limvar_TEM_impact'

# cmd args
qi    = int(sys.argv[1])
dry   = bool(int(sys.argv[2]))
print('args: {}'.format(sys.argv))

# ----------------------------------------------------------------

# get ensemble data
print('locating data for qi {}...'.format(qi))
qstr  = ['','_TRACER-AOA','_TRACER-E90j'][int(qi)]
enshist = sorted(glob.glob('{}/*ens*.eam.h1*TEM*L45{}_impact.nc'.format(temloc, qstr)))
print('num files: {}'.format(len(enshist)))
print('first file: {}'.format(enshist[0]))
print('last file: {}'.format(enshist[-1]))
if(dry): exit(0)

print('reading data...')
enshistdat = [xr.open_dataset(dat) for dat in enshist]

print('averaging data...')
N = len(enshistdat)
ensmean = xr.zeros_like(enshistdat[0])
for i in range(N):
    print('ens {}...'.format(i+1))
    ensmean = ensmean + enshistdat[i]
ensmean = ensmean / N

print('finding data deviation...')
N = len(enshistdat)
ensstd = xr.zeros_like(enshistdat[0])
for i in range(N):
    print('ens {}...'.format(i+1))
    ensstd = ensstd + (enshistdat[i] - ensmean)**2
ensstd = np.sqrt(ensstd / N)

print('writing out...')
name = enshist[0].split('/')[-1].split('ens')[0]
ensmean.to_netcdf('{}/{}_impact_ensmean{}.nc'.format(outdir, name, qstr))
ensstd.to_netcdf('{}/{}_impact_ensstd{}.nc'.format(outdir, name, qstr))
