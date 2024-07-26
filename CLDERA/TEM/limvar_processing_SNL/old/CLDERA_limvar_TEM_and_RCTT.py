import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import PyTEMDiags as pt
import pdb
import glob
import climate_toolbox as ctb

# directory for writing out processed data, TEM results, spectral mapping files
outdir = '/ascldap/users/jphollo/data/limvar_TEM'

def run_TEM(data, name=None, maxlev=500):

    run_name = name
    L = 45
    data=data.sel({'lev':slice(0, maxlev)})
    p0 = float(data['P0'].values)
    
    print('\n\n ----------  computing TEM for file {}...'.format(run_name))
    ua, va, ta, wap, ps, lat, lev = data['U'], data['V'], data['T'], data['OMEGA'],\
                                    data['PS'], data['lat'], data['lev']
    q = [data['AOA'], data['E90j']]
    #p = ctb.compute_hybrid_pressure(ps, data['hyam'], data['hybm'], dims_like=ta, p0=p0)
    p = lev * 100
    
    tem = pt.TEMDiagnostics(ua, va, ta, wap, p, lat, q=q, p0=p0, L=L,
                            overwrite_map=False, debug_level=2, grid_name='ne30pg2')
    return [tem.to_netcdf(loc=outdir, prefix=run_name, include_attrs=False), 
            tem.q_to_netcdf(loc=outdir, prefix=run_name, include_attrs=False)]

# ----------------------------------------------------------------

# set max time index. The hist files are monthly, so this is a number of months 
# to include in the analysis
tmax = 38
tmax = 12
# set max ensemble members
Nmax = 5

# get ensemble data
print('locating data...')
loc = '/gpfs/cldera/data/e3sm/e3sm-cldera/CLDERA-historical/limvar_src_tag'
ensdirs = sorted(glob.glob('{}/*10Tg.ens*/'.format(loc)))[:Nmax]
enshist = [sorted(glob.glob('{}/archive/atm/hist/*h1*'.format(ensdir)))[:tmax] for ensdir in ensdirs]
cfdirs = sorted(glob.glob('{}/*10Tg.ens*cf/'.format(loc)))[:Nmax]
cfhist = [sorted(glob.glob('{}/archive/atm/hist/*h1*'.format(cfdir)))[:tmax] for cfdir in cfdirs]
N = len(ensdirs)

# compute TEM per-month, per-ensemble member
ensfile = '{}/limvar_src_tag_ens_TEM_tmax{}_Nmax{}.nc'.format(outdir, tmax, Nmax)
overwrite = False
try:
    if(overwrite): raise FileNotFoundError
    ens = xr.open_dataset(ensfile)
except FileNotFoundError:
    print('computing ensemble TEM...')
    enstem = enshist
    for i in range(N):
        print('ens {}/{}...'.format(i+1, N))
        for j in range(len(enshist[i])):
            print('...hist {}/{}...'.format(j+1, len(enshist[i])))
            hist = xr.open_dataset(enshist[i][j])
            name = enshist[i][j].split('/')[-1].split('.nc')[0]
            tem = run_TEM(hist, name)
            enstem[i][j] = xr.merge([xr.open_dataset(tem[0]), xr.open_dataset(tem[1])])
        print('concatenating in time ({} gb)...'.format(np.sum([h.nbytes for h in enstem[i]])/1e9))
        enstem[i] = xr.concat(enstem[i], dim='time')
    print('concatenating ensemble ({} gb)...'.format(np.sum([d.nbytes for d in dat])/1e9))
    enstem = xr.concat(enstem, dim='ens')
    print('writing data...')
    ens.to_netcdf(ensfile)

pdb.set_trace()

print('reading counterfactual data...')
dat = [None]*N
for i in range(N):
    print('cf {}/{}...'.format(i+1, N))
    hist = [None]*len(cfhist[i])
    for j in range(len(hist)):
        print('...hist {}/{}...'.format(j+1, len(hist)))
        hist[j] = xr.open_dataset(cfhist[i][j])[varlist]
    print('concatenating in time ({} gb)...'.format(np.sum([h.nbytes for h in hist])/1e9))
    dat[i] = xr.concat(hist, dim='time')
print('concatenating ensemble ({} gb)...'.format(np.sum([d.nbytes for d in dat])/1e9))
cf = xr.concat(dat, dim='ens')
print('writing data...')
cf.to_netcdf('{}/limvar_src_tag_cf_tmax{}_Nmax{}.nc'.format(out_dir, tmax, Nmax))

pdb.set_trace()

#  combine data...

exit(0)

#  combine data...
outdir = '/ascldap/users/jphollo/data/limvar_TEM'
for i in range(len(daily)):

    overwrite_tem = False
    run_name = daily[i].split('/')[-1].split('.nc')[0]
    try:
        if(overwrite_tem): raise IndexError
        tem_file = glob.glob('{}/{}_TEM*.nc'.format(outdir, run_name))[0]
        tem.append(tem_file)
        print('\n\n ---------- read {} TEM from file {}...'.format(run_name, tem_file.split('/')[-1]))
    except IndexError:
        print('\n\n ----------  computing TEM for file {}...'.format(run_name))
        data = xr.open_dataset(daily[i])
        overwrite_map = False
        L = 45

        data=data.sel({'lev':slice(0, 170)})
        p0 = float(data['P0'].values)
        
        ua, va, ta, wap, ps, lat = data['U'], data['V'], data['T'],\
                                   data['OMEGA'], data['PS'], data['lat']
        p = ctb.compute_hybrid_pressure(ps, data['hyam'], data['hybm'], dims_like=ta, p0=p0)
        weights = np.ones(len(lat))
        
        tem = pt.TEMDiagnostics(ua, va, ta, wap, p, lat, p0=p0, L=L,
                                overwrite_map=overwrite_map, debug=True,
                                grid_name='ne30pg2')
        tem.to_netcdf(loc=outdir, prefix=run_name, include_attrs=True)


