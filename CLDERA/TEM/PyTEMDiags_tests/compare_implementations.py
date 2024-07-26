# Joe Hollowed
# University of Michigan 2023
#
# Code to perform test runs of the PyTEM package, as well as compare results
# to implementations by Tom Ehrmann and Christiane Jablonowski

import pdb
import glob
import numpy as np
import PyTEMDiags
import xarray as xr
import climate_toolbox as ctb
import matplotlib.pyplot as plt
import climate_artist as cla
from matplotlib.ticker import ScalarFormatter

top = '/pscratch/sd/j/jhollo/PyTEMDiags/output'
jhf = sorted(glob.glob('{}/TEM*180x360*.nc'.format(top)))
cj = xr.open_dataset(glob.glob('{}/*p_converted_z*'.format(top))[0])
ln = xr.open_dataset(glob.glob('{}/*lev1-33*'.format(top))[0])

# get PyTEM data
p_jh_var_or_coord = 'coord'
#p_jh_var_or_coord = 'var'
if(p_jh_var_or_coord == 'coord'): jh = xr.open_dataset(jhf[0])
else: jh = xr.open_dataset(jhf[1])

jh_nt, cj_nt, ln_nt = [len(jh['time']), len(cj['time']), len(ln['time'])]
jh_nlev, cj_nlev, ln_nlev = [len(jh['lev']), len(cj['lev']), len(ln['lev'])]
jh_nlat, cj_nlat, ln_nlat = [len(jh['lat']), len(cj['lat']), len(ln['lat'])]
print('nt: {}, {}, {}'.format(jh_nt, cj_nt, ln_nt))
print('nlev: {}, {}, {}'.format(jh_nlev, cj_nlev, ln_nlev))
print('nlat: {}, {}, {}'.format(jh_nlat, cj_nlat, ln_nlat))

# take monthly mean
jh = jh.mean('time')
ln = ln.mean('time')
cj = cj.mean(['time', 'lon'])

#limit to cj lat
ln['lat'] = jh.lat
jh = jh.sel(lat=cj.lat)
ln = ln.sel(lat=cj.lat)

# match lev to cj
ln = ln.sel(lev=ln.lev[1:-1])
jh = jh.sel(lev=cj.lev)

# verify
jh_nlev, cj_nlev, ln_nlev = [len(jh['lev']), len(cj['lev']), len(ln['lev'])]
jh_nlat, cj_nlat, ln_nlat = [len(jh['lat']), len(cj['lat']), len(ln['lat'])]
print('mod:')
print('nlev: {}, {}, {}'.format(jh_nlev, cj_nlev, ln_nlev))
print('nlat: {}, {}, {}'.format(jh_nlat, cj_nlat, ln_nlat))

LAT, LEV = np.meshgrid(jh['lat'], jh['lev'])
LAT, LEV = LAT.T, LEV.T

allvars = ['utendepfd', 'vtem', 'wtem']
titles = ['utendepfd [m/s2 * 10^-5]', 'vtem [m/s]', 'wtem [m/s]']
scaling = [1e5, 1, 1]
cmaps = [plt.cm.rainbow, plt.cm.rainbow, plt.cm.rainbow]
levels = [np.arange(-10, 6.25, 1.25), 
          np.linspace(-0.4, 1.6, 11), 
          np.linspace(-0.005, 0.005, 11)] 
difflevels_cj = [11, 11, 11]
difflevels_ln = [11, 11, 11]

for i in range(len(allvars)):
    fig = plt.figure(figsize=(15, 10))
    axjh = fig.add_subplot(2, 3, 4)
    axcj = fig.add_subplot(2, 3, 2)
    axln = fig.add_subplot(2, 3, 3)
    axjh_cj = fig.add_subplot(2, 3, 5)
    axjh_ln = fig.add_subplot(2, 3, 6)
    allax = [axjh, axcj, axln, axjh_cj, axjh_ln]

    var_jh = jh[allvars[i]] * scaling[i] 
    var_cj = cj[allvars[i]] * scaling[i]
    var_ln = ln[allvars[i]] * scaling[i]
    var_jh_cj = var_jh - var_cj
    var_jh_ln = var_jh - var_ln

    varjh = var_jh.transpose('lat', 'lev')
    var_cj = var_cj.transpose('lat', 'lev')
    var_ln = var_ln.transpose('lat', 'lev')
    var_jh_cj = var_jh_cj.transpose('lat', 'lev')
    var_jh_ln = var_jh_ln.transpose('lat', 'lev')

    # --- vars
    im=axjh.contourf(LAT, LEV, var_jh, cmap=cmaps[i], levels=levels[i], extend='both')
    axjh.contour(LAT, LEV, var_jh, colors='k', levels=levels[i], linewidths=0.5)
    plt.colorbar(im, ax=axjh, location='bottom', extend='both')
   
    im = axcj.contourf(LAT, LEV, var_cj, cmap=cmaps[i], levels=levels[i], extend='both')
    axcj.contour(LAT, LEV, var_cj, colors='k', levels=levels[i], linewidths=0.5)
    plt.colorbar(im, ax=axcj, location='top', extend='both')
    
    im=axln.contourf(LAT, LEV, var_ln, cmap=cmaps[i], levels=levels[i], extend='both')
    axln.contour(LAT, LEV, var_ln, colors='k', levels=levels[i], linewidths=0.5)
    plt.colorbar(im, ax=axln, location='top', extend='both')
     
    # --- vars diff
    im = axjh_cj.contourf(LAT, LEV, var_jh_cj, cmap=cmaps[i], 
                          levels=difflevels_cj[i], extend='both')
    axjh_cj.contour(LAT, LEV, var_jh_cj, colors='k', levels=difflevels_cj[i], 
                    linewidths=0.5)
    plt.colorbar(im, ax=axjh_cj, location='bottom', extend='both')
    
    im = axjh_ln.contourf(LAT, LEV, var_jh_ln, cmap=cmaps[i], 
                          levels=difflevels_ln[i], extend='both')
    axjh_ln.contour(LAT, LEV, var_jh_ln, colors='k', levels=difflevels_ln[i], 
                    linewidths=0.5)
    plt.colorbar(im, ax=axjh_ln, location='bottom', extend='both')
    
    # --- format
    axjh.set_xlabel('lat [deg]')
    axjh_cj.set_xlabel('lat [deg]')
    axjh_ln.set_xlabel('lat [deg]')
    axcj.set_xlabel('lev [hPa]')
    axln.set_xlabel('lev [hPa]')
    axln.yaxis.tick_right()
    axln.yaxis.set_label_position("right")
    axjh_ln.set_xlabel('lev [hPa]')
    axjh_ln.yaxis.tick_right()
    axjh_ln.yaxis.set_label_position("right")

    axjh.set_title('PyTEMDiags')
    axcj.set_title('Fortran')
    axln.set_title('Matlab')
    axjh_cj.set_title('(Python - Fortran)')
    axjh_ln.set_title('(Python - Matlab)')

    for ax in allax: 
        ax.grid(alpha=0.33)
        ax.set_ylim([1, np.max(LEV)])
        ax.invert_yaxis()
        ax.set_yscale('log')
        ax.yaxis.set_major_formatter(ScalarFormatter())
    plt.suptitle(titles[i], fontsize=14)
       
    plt.tight_layout()
    plt.savefig('figs/{}_{}_compare_weighted.png'.format(allvars[i], p_jh_var_or_coord))
