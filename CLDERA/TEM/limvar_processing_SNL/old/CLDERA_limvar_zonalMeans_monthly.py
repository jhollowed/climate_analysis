'''
This file computes and writes out zonal averages on isobars for the 
following monthly-averaged variables. First the vertical remapping is done, 
and then the zonal average is taken via a spectral method.

- BUTGWSPEC
- UTGWORO
- UTGWSPEC
- U
- T
- AOA
- E90
'''

import os.path
import sys
import pdb
import glob
import dask
import numpy as np
import xarray as xr
import PyTEMDiags as pt
from geocat.comp import interpolation as gcinterp

# output directory
outdir = '/gpfs/cldera/data/e3sm/e3sm-cldera/CLDERA-historical/limvar_src_tag/tem_processed/interp_to_pressure_zonalMeans'
# data location
loc    = '/gpfs/cldera/data/e3sm/e3sm-cldera/CLDERA-historical/limvar_src_tag/' 
# TEM data
temloc = '/ascldap/users/jphollo/data/limvar_TEM' 
# map directory for zonal averager
mapdir = '/ascldap/users/jphollo/repos/PyTEMDiags/maps'

Nens  = int(sys.argv[1])
tmini = int(sys.argv[2])
tmaxi = int(sys.argv[3])
cfb   = bool(int(sys.argv[4]))
dry   = bool(int(sys.argv[5]))
print('args: {}'.format(sys.argv))

overwrite=0
g0 = 9.80665

# get ensemble data
print('locating data...')
ensdirs = np.array(sorted(glob.glob('{}/*10Tg.ens*[0-9]/'.format(loc))))
nohud    = np.array(['hud.' not in ss for ss in ensdirs])
ensdirs = ensdirs[nohud]
nocf    = np.array(['.cf.' not in ss for ss in ensdirs])
ensdirs = ensdirs[nocf]
mask    = np.array(['ens{}'.format(Nens) in ss for ss in ensdirs])
ensdir  = ensdirs[mask][0]
enshist = sorted(glob.glob('{}/archive/atm/hist/*h0*'.format(ensdir)))[tmini:tmaxi+1]

# get counterfactual data
cfdirs  = np.array(sorted(glob.glob('{}/*ens*cf/'.format(loc))))
nohud    = np.array(['hud.' not in ss for ss in cfdirs])
cfdirs = cfdirs[nohud]
mask    = np.array(['ens{}'.format(Nens) in ss for ss in cfdirs])
cfdir   = cfdirs[mask][0]
cfhist  = sorted(glob.glob('{}/archive/atm/hist/*h0*'.format(cfdir)))[tmini:tmaxi+1]

# replace ensemble data with counterfactual data if cfb=true
if(cfb):
    ensdirs=cfdirs
    ensdir=cfdir
    enshist=cfhist

# do interpolation
print('-------- working on ENS {}, times {}-{}, cf {} --------'.format(Nens, tmini, tmaxi, cfb))
print('first file is: {}'.format(enshist[0].split('/')[-1]))
print('last file is: {}'.format(enshist[-1].split('/')[-1]))

for i in range(len(enshist)):
   
    # get filename
    fname = enshist[i].split('/')[-1]
    print('\n ----- working on {}'.format(fname))
 
    # get TEM file (so that we know the output latitudes for the zonal average)
    temfile = sorted(glob.glob('{}/{}*'.format(temloc, fname.split('.nc')[0].replace('h0', 'h1'))))[0]
    temfname = temfile.split('/')[-1]
    print('TEM reference file is {}'.format(temfname))
    if(not os.path.isfile(temfile)):
        raise RuntimeError('tem file does not exist! Aborting...')

    # define output file
    outfile = '{}/{}_pinterp_zonalmeans.nc'.format(outdir, fname.split('.nc')[0])
    if(os.path.isfile(outfile) and not overwrite):
        print('output exists; skipping...')
        continue

    # ------ open data
    data    = xr.open_dataset(enshist[i])
    temdata = xr.open_dataset(temfile)

    p0      = float(data['P0'].values)
    lat     = data['lat']
    lev     = data['lev']
    lat_out = temdata['lat']

    # ------ interpolate to isobars
    # the gravity wave variables were not previously vertically interpolated, so do that now
    print('Generating interpolation...')
    if(not dry):
        interp_args = {'ps':data['PS'], 'hyam':data['hyam'], 'hybm':data['hybm'], 'p0':p0,
                       'new_levels':lev.values*100, 'method':'log', 'extrapolate':True, 
                       'variable':'other'}
        gw1  = gcinterp.interp_hybrid_to_pressure(data['BUTGWSPEC'], **interp_args)
        gw2  = gcinterp.interp_hybrid_to_pressure(data['UTGWORO'], **interp_args)
        gw3  = gcinterp.interp_hybrid_to_pressure(data['UTGWSPEC'], **interp_args)
        gw1['plev'] = gw1.plev/100
        gw2['plev'] = gw2.plev/100
        gw3['plev'] = gw3.plev/100

        ua  = gcinterp.interp_hybrid_to_pressure(data['U'], **interp_args)
        aoa = gcinterp.interp_hybrid_to_pressure(data['AOA'], **interp_args)
        e90 = gcinterp.interp_hybrid_to_pressure(data['E90j'], **interp_args)
        ua['plev'] = ua.plev/100
        aoa['plev'] = aoa.plev/100
        e90['plev'] = e90.plev/100

        interp_args['variable'] = 'temperature'
        interp_args['t_bot']    = data['T'].isel(lev=-1)
        interp_args['phi_sfc']  = data['Z3'].isel(lev=-1) * g0
        ta = gcinterp.interp_hybrid_to_pressure(data['T'], **interp_args)
        ta['plev'] = ta.plev/100

    print('Computing interpolation...')
    if(not dry):
        gw1, gw2, gw3, ua, ta, aoa, e90 = dask.compute(gw1, gw2, gw3, ua, ta, aoa, e90)
        plev = gw1['plev']
    
    # ------- do zonal averaging
    print('---- doing zonal averaging...')
    
    # create zonal averager
    L  = 90
    ZM = pt.sph_zonal_averager(lat, lat_out, L, grid_name='180x360', 
                               save_dest=mapdir, debug=True)
    ZM.sph_compute_matrices()
    
    # zonally average
    if(not dry):
        zm_gw1 = ZM.sph_zonal_mean(gw1.T)
        zm_gw2 = ZM.sph_zonal_mean(gw2.T)
        zm_gw3 = ZM.sph_zonal_mean(gw3.T)
        zm_ua  = ZM.sph_zonal_mean(ua.T)
        zm_ta  = ZM.sph_zonal_mean(ta.T)
        zm_aoa = ZM.sph_zonal_mean(e90.T)
        zm_e90 = ZM.sph_zonal_mean(aoa.T)

    # write out
    print('writing to {}'.format(outfile.split('/')[-1]))
    if(not dry):
        data = xr.merge([zm_gw1, zm_gw2, zm_gw3, 
                         zm_ua, zm_ta, zm_aoa, zm_e90, 
                         data['time_bnds']])
        data.to_netcdf(outfile) 
