'''
Joe Hollowed 2024

This computes the Transformed Eulerian Mean (TEM) and asssocaited quantities for a limvar ensemble
member. Command line arguments are:

Usage
-----
python ./limvar_TEM.py [Nens] [massMag] [tmini] [tmaxi] [overwrite] [dry]

Parameters
----------
Nens : int
    the ensemble member specified as an integer
massMag : int
    the mass magnitude of the ensemble. Valid values are:
    0, 1, 3, 5, 7, 10, 13, 15
    if massMag = 0, this sepcifies the counterfactual ensemble
tmini : int
    index of the first history file to include in the calculation
tmaxi : int
    index of the last history file to include in the calculation
overwrite : int
    whether or not to overwrite the TEM data if it already exists
    at the specified outpit file. 0 = False, 1 = True, 
    If False, then the script skips the interpolation for any file that 
    already exists
dry : int
    whether or not to do a "dry run" of the function, printing information
    about the execution but skipping the actual interpolation. 0 = False, 1 = True

Note that tmini, tmaxi specify only the *index* of the identified time files, and do
not represent a time in physical units 
'''

import pdb
import sys
import glob
import numpy as np
import xarray as xr
import PyTEMDiags as pt
from combine_zinterp_native import combine_interp_native_data

# directory for writing out processed data, TEM results
outdir = '/gpfs/cldera/data/e3sm/e3sm-cldera/CLDERA-historical/limvar_src_tag/tem_processed/tem_output'
# vertically interpolated data location (all data below hybm=0 transition)
interploc = '/gpfs/cldera/data/e3sm/e3sm-cldera/CLDERA-historical/limvar_src_tag/tem_processed/interp_to_pressure_3D'

# cmd args
Nens      = int(sys.argv[1])
massMag   = int(sys.argv[2])
tmini     = int(sys.argv[3])
tmaxi     = int(sys.argv[4])
overwrite = bool(int(sys.argv[5]))
dry       = bool(int(sys.argv[6]))
print('args: {}'.format(sys.argv))

cfb = massMag == 0 # counterfactual flag 
L   = 45           # spherical harmonic max order for TEM zonal averaging

# ----------------------------------------------------------------

# get ensemble member data
print('locating data...')
if(not cfb):
    enshist = sorted(glob.glob('{}/*{}Tg*ens{}.eam*h1*'.format(
                               interploc, massMag, Nens)))[tmini:tmaxi+1]
else:
    enshist = sorted(glob.glob('{}/*ens{}.cf.*h1*'.format(
                               interploc, Nens)))[tmini:tmaxi+1]

# compute TEM for each history file in this time range
print('computing ensemble TEM...')
for i in range(len(enshist)):
    print('\n ========== working on time {} ({}/{})'.format(tmini+i, i, tmaxi-tmini))
    
    name = enshist[i].split('/')[-1].split('.nc')[0]
    
    # check if file exists
    existing = glob.glob('{}/{}*TEM*'.format(outdir, name))
    if(len(existing) > 0 and not overwrite):
        print('TEM file exists; skipping...')
        continue
    
    # call helper function to combine vertically interpolated and native data
    print('merging interpolated and native data...')
    if(not dry): hist = combine_interp_native_data(enshist[i])
    else: hist = None
    
    # run TEM
    print('---- computing TEM for file {}...'.format(name))
    if(dry): exit(0)
     
    ua, va, ta, wap, lat, p = hist['U'], hist['V'], hist['T'], hist['OMEGA'],\
                              hist['lat'], hist['plev']*100
    q                       = [hist['AOA'], hist['E90j']]
    p0                      = float(hist['P0'].values)
  
    tem = pt.TEMDiagnostics(ua, va, ta, wap, p, lat, q=q, p0=p0, L=L,
                            overwrite_map=False, debug_level=1, grid_name='ne30pg2')

    # run TEM
    tem.to_netcdf(loc=outdir, prefix=name, include_attrs=False)
    # run TEM for tracers
    tem.q_to_netcdf(loc=outdir, prefix=name, include_attrs=False)
