import pdb
import sys
import glob
import numpy as np
import xarray as xr
import PyTEMDiags as pt
import scipy.stats as stats

# output location
outdir = '/ascldap/users/jphollo/data/limvar/limvar_zonalMeans_daily_ttest'
# tem results location
dataloc = '/ascldap/users/jphollo/data/limvar/limvar_zonalMeans_daily_impact_ensmean'

# cmd args
dry   = bool(int(sys.argv[1]))
print('args: {}'.format(sys.argv))

# ----------------------------------------------------------------

# get ensemble data
print('locating data...')
impactMeanHist = sorted(glob.glob('{}/*ensmean.nc'.format(dataloc)))[0]
impactStdHist  = sorted(glob.glob('{}/*ensstd.nc'.format(dataloc)))[0]
print('mean file: {}'.format(impactMeanHist))
print('std file: {}'.format(impactStdHist))
if(dry): exit(0)

print('reading data...')
impactMean = xr.open_dataset(impactMeanHist)
impactStd = xr.open_dataset(impactStdHist)

impactdir = '/ascldap/users/jphollo/data/limvar/limvar_zonalMeans_daily_impact'
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
pval.to_netcdf('{}/{}_pval.nc'.format(outdir, name))
ttest.to_netcdf('{}/{}_ttest.nc'.format(outdir, name))
