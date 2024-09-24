import pdb
import sys
import glob
import numpy as np
import xarray as xr
import PyTEMDiags as pt
import scipy.stats as stats

# output location
outdir = '/ascldap/users/jphollo/data/limvar/limvar_TEM_ttest'
# tem results location
temloc = '/ascldap/users/jphollo/data/limvar/limvar_TEM_impact_ensmean'

# cmd args
qi    = int(sys.argv[1])
dry   = bool(int(sys.argv[2]))
print('args: {}'.format(sys.argv))

# ----------------------------------------------------------------

# get ensemble data
print('locating data for qi {}...'.format(qi))
qstr  = ['','_TRACER-AOA','_TRACER-E90j'][int(qi)]
impactMeanHist = sorted(glob.glob('{}/*ensmean{}.nc'.format(temloc, qstr)))[0]
impactStdHist  = sorted(glob.glob('{}/*ensstd{}.nc'.format(temloc, qstr)))[0]
print('mean file: {}'.format(impactMeanHist))
print('std file: {}'.format(impactStdHist))
if(dry): exit(0)

print('reading data...')
impactMean = xr.open_dataset(impactMeanHist)
impactStd = xr.open_dataset(impactStdHist)

impactdir = '/ascldap/users/jphollo/data/limvar/limvar_TEM_impact'
N = len(glob.glob('{}/*AOA*'.format(impactdir)))

print('doing ttest...')
ttest = impactMean / (impactStd * np.sqrt(N))

print('getting pvalue...')
pval = xr.zeros_like(ttest)
pval = pval.drop_vars(['time', 'lat', 'plev'])
for v in list(pval.variables):
    print('{}...'.format(v))
    pval[v].values = stats.t.sf(np.abs(ttest[v]), N-1) * 2

print('writing out...')
name = impactMeanHist.split('/')[-1].split('.nc')[0]
pval.to_netcdf('{}/{}_pval{}.nc'.format(outdir, name, qstr))
ttest.to_netcdf('{}/{}_ttest{}.nc'.format(outdir, name, qstr))
