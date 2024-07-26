import os.path
import sys
import pdb
import glob
import dask
import numpy as np
import xarray as xr
import PyTEMDiags as pt
from geocat.comp import interpolation as gcinterp

outdir = '/gpfs/cldera/data/e3sm/e3sm-cldera/CLDERA-historical/limvar_src_tag/tem_processed/interp_to_pressure_zonalMeans'
ploc   = '/gpfs/cldera/data/e3sm/e3sm-cldera/CLDERA-historical/limvar_src_tag/tem_processed/interp_to_pressure_3D'
loc    = '/gpfs/cldera/data/e3sm/e3sm-cldera/CLDERA-historical/limvar_src_tag/'
mapdir = '/ascldap/users/jphollo/repos/PyTEMDiags/maps'
temloc = '/ascldap/users/jphollo/data/limvar/limvar_TEM_timeConcat'

Nens  = int(sys.argv[1])
tmini = int(sys.argv[2])
tmaxi = int(sys.argv[3])
cfb   = bool(int(sys.argv[4]))
dry   = bool(int(sys.argv[5]))
print('args: {}'.format(sys.argv))

# get ensemble data
print('locating data...')
ensdirs = np.array(sorted(glob.glob('{}/*10Tg.ens*[0-9]/'.format(loc))))
mask    = np.array(['ens{}'.format(Nens) in ss for ss in ensdirs])
ensdir  = ensdirs[mask][0]
enshist = sorted(glob.glob('{}/archive/atm/hist/*h1*'.format(ensdir)))[tmini:tmaxi+1]

cfdirs  = np.array(sorted(glob.glob('{}/*ens*cf/'.format(loc))))
mask    = np.array(['ens{}'.format(Nens) in ss for ss in cfdirs])
cfdir   = cfdirs[mask][0]
cfhist  = sorted(glob.glob('{}/archive/atm/hist/*h1*'.format(cfdir)))[tmini:tmaxi+1]
if(cfb):
    ensdirs=cfdirs
    ensdir=cfdir
    enshist=cfhist

# do interpolation
print('-------- working on ENS {}, times {}-{}, cf {} --------'.format(Nens, tmini, tmaxi, cfb))
print('first file is: {}'.format(enshist[0].split('/')[-1]))
print('last file is: {}'.format(enshist[-1].split('/')[-1]))

for i in range(len(enshist)):
    
    fname = enshist[i].split('/')[-1]
    print('\n ----- working on {}'.format(fname))
    
    pfile = '{}/{}_pinterp.nc'.format(ploc, fname.split('.nc')[0])
    pfname = pfile.split('/')[-1]
    print('pressure-interpolated file is {}'.format(pfname))
    if(not os.path.isfile(pfile)):
        raise RuntimeError('pressure-interpolated file does not exist! Aborting...')

    temfile = sorted(glob.glob('{}/{}*'.format(temloc, pfname.split('.nc')[0])))[0]
    temfname = temfile.split('/')[-1]
    print('TEM reference file is {}'.format(temfname))
    if(not os.path.isfile(temfile)):
        raise RuntimeError('tem file does not exist! Aborting...')

    outfile = '{}/{}_zonalmeans.nc'.format(outdir, pfname.split('.nc')[0])
    if(os.path.isfile(outfile)):
        print('output exists; skipping...')
        continue

    # ------ open data
    temdata = xr.open_dataset(temfile)
    pdata   = xr.open_dataset(pfile)
    data    = xr.open_dataset(enshist[i])

    p0      = float(data['P0'].values)
    lat     = data['lat']
    lev     = pdata['plev']
    lat_out = temdata['lat']

    # create zonal averager
    L  = 90
    ZM = pt.sph_zonal_averager(lat, lat_out, L, grid_name='180x360', 
                               save_dest=mapdir, debug=True)
    ZM.sph_compute_matrices()
    
    # zonally average
    print('---- doing zonal averaging...')
    if(not dry):
        print('u...')
        zm_ua    = ZM.sph_zonal_mean(pdata['U'].T)
        print('v...')
        zm_va    = ZM.sph_zonal_mean(pdata['V'].T)
        print('t...')
        zm_ta    = ZM.sph_zonal_mean(pdata['T'].T)
        print('wap...')
        zm_wap   = ZM.sph_zonal_mean(pdata['OMEGA'].T)
        print('aoa...')
        zm_aoa   = ZM.sph_zonal_mean(pdata['AOA'].T)
        print('e90...')
        zm_e90   = ZM.sph_zonal_mean(pdata['E90j'].T)

    # write out
    print('writing to {}'.format(outfile.split('/')[-1]))
    if(not dry):
        data = xr.merge([zm_ua, zm_va, zm_ta, zm_wap, zm_aoa, zm_e90])
        data.to_netcdf(outfile) 
