# eruption_movies.py
# Joe Hollowed
# CLIMATE 589, 4/15/2022
#
# renders video frames of the eruption in horizontal cross, vertical cross, 
# and AzimuthalEquidistant projection

import pdb
import glob
import numpy as np
import xarray as xr
import matplotlib as mpl
import cartopy.crs as ccrs
import artist_utils as claut
import climate_toolbox as ctb
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
from climate_artist import vertical_slice as pltvert
from climate_artist import horizontal_slice as plthor

# ============================================================


def circulation_snapshot(runf, tsnap, title, savedest, inclTracers=True):
   
    if(isinstance(tsnap, list)):
        ti, tf = tsnap[0], tsnap[1]
        tsnap = slice(ti, tf)
        tismean=True
    else:
        tismean=False

    # params
    lat0 = 15.15
    lon0 = 120.35
    dlat = 0.5
    psel = 45
    minc = -8
    clipc = -15
    maxc = -4
    clevels = 9

    # read data, transform time coordinate to ndays
    dat = xr.open_dataset(runf)
    td = ctb.time2day(dat['time'])
    dat = dat.assign_coords(time=td)
     
    # get vars at time snapshot/slice, take zonal means
    if(tismean):
        method = None
        mean = ['lon', 'time']
        ttitle = 'mean of days {}-{}'.format(ti, tf)
    else: 
        method = 'nearest'
        mean = 'lon'
        ttitle = 'day {}'.format(tsnap)
   
    U   = dat['U'].sel({'time':tsnap}, method=method).mean(mean)
    V   = dat['V'].sel({'time':tsnap}, method=method).mean(mean)
    T   = dat['T'].sel({'time':tsnap}, method=method).mean(mean)
    #OM  = dat['OMEGA'].sel({'time':tsnap}, method=method).mean(mean)
    if inclTracers:
        C   = dat['SAI_SO2'].sel({'time':tsnap}, method=method).sel({'lon':120.35}, method='nearest')
    lat = dat['lat']
    lev = dat['lev']

    # if tsnap was a slice, take time mean
    if len(U.shape) == 3:
        U = U.mean('time')
        V = V.mean('time')
        #OM = OM.mean('time')

    # ---------- plot ----------
    data_crs = ccrs.PlateCarree()
    #levels = np.linspace(minc, maxc, clevels)
    levels=10
    uleva = np.arange(0, np.max(U), 15)
    ulev = np.hstack([-uleva[1:][::-1], uleva])
    
    Tlev = np.arange(170, 310, 15)

    #cmap = claut.ncar_rgb_to_cmap(gmt)
    cmap = mpl.cm.rainbow
    #cmap = mpl.cm.Spectral_r
        
    print('\n\n=============== {}'.format(title)) 

    #fig = plt.figure(figsize=(8,10))
    
    figu = plt.figure()
    axu = figu.add_subplot(111)
    #axu = fig.add_subplot(221)
    #axv = fig.add_subplot(222)
    #axom = fig.add_subplot(223)
    
    figT = plt.figure()
    axT = figT.add_subplot(111)

    # ----- vertical
    
    # ============= U

    pltargs = {'levels':ulev, 'cmap':cmap, 'zorder':0}
    pltargs_c = {'levels':ulev, 'colors':'k', 'linewidths':0.6, 'zorder':1}
    cArgs = {'orientation':'horizontal', 'location':'top', 'label':'U [m/s]'}
    var_dict = [{'var':U, 'plotType':'contourf', 'plotArgs':pltargs, 'colorArgs':cArgs}]
    var_dict_c = [{'var':U, 'plotType':'contour', 'plotArgs':pltargs_c, 'colorFormatter':None}]    
    cf = pltvert(lat, lev, var_dict, ax=axu, plot_zscale=True, slice_at='zonal mean', xlabel='')
    pltvert(lat, lev, var_dict_c, ax=axu, plot_zscale=False, inverty=False, slice_at='', xlabel='')
    cf[0].set_ticks(ulev)
    axu.set_ylabel('p  [hPa]')
    axu.set_xlabel('latitude')
    figu.suptitle('{}'.format(title), fontsize=14)
   
    if inclTracers:
        #pltargs_c = {'levels':[-5, -4], 'colors':'m', 'linewidths':1.6, "linestyles":'-'}
        pltargs_c = {'levels':[-7, -5.5], 'colors':'m', 'linewidths':1.6, "linestyles":'-'}
        cargs = {'fmt':'%.0f', 'manual':[(18.8, 201), (14.89, 26.2)]} 
        var_dict_c = [{'var':np.log10(C), 'plotType':'contour', 'plotArgs':pltargs_c, 'colorArgs':cargs}]
        pltvert(lat, lev, var_dict_c, ax=axu, plot_zscale=False,inverty=False,slice_at='',xlabel='')
        axu.plot([0,0],[100,100],'-m', label='log10(concentration)')
        axu.legend(loc='lower right', fancybox=False, framealpha=1)
    
    # ============= T
    
    pltargs = {'levels':Tlev, 'cmap':cmap, 'zorder':0}
    pltargs_c = {'levels':Tlev, 'colors':'k', 'linewidths':0.6, 'zorder':1}
    cArgs = {'orientation':'horizontal', 'location':'top', 'label':'T [K]'}
    var_dict = [{'var':T, 'plotType':'contourf', 'plotArgs':pltargs, 'colorArgs':cArgs}]
    var_dict_c = [{'var':T, 'plotType':'contour', 'plotArgs':pltargs_c, 'colorFormatter':None}]    
    cf = pltvert(lat, lev, var_dict, ax=axT, plot_zscale=True, slice_at='zonal mean', xlabel='')
    pltvert(lat, lev, var_dict_c, ax=axT, plot_zscale=False, inverty=False, slice_at='', xlabel='')
    cf[0].set_ticks(Tlev)
    axT.set_ylabel('p  [hPa]')
    axT.set_xlabel('latitude')
    figT.suptitle('{}'.format(title), fontsize=14)
   
    if inclTracers:
        #pltargs_c = {'levels':[-5, -4], 'colors':'m', 'linewidths':1.6, "linestyles":'-'}
        pltargs_c = {'levels':[-7, -5.5], 'colors':'m', 'linewidths':1.6, "linestyles":'-'}
        cargs = {'fmt':'%.0f', 'manual':[(18.8, 201), (14.89, 26.2)]} 
        var_dict_c = [{'var':np.log10(C), 'plotType':'contour', 'plotArgs':pltargs_c, 'colorArgs':cargs}]
        pltvert(lat, lev, var_dict_c, ax=axT, plot_zscale=False,inverty=False,slice_at='',xlabel='')
        axu.plot([0,0],[100,100],'-m', label='log10(concentration)')
        axu.legend(loc='lower right', fancybox=False, framealpha=1)

    
    if(0):
        pltargs = {'levels':vlev, 'cmap':cmap, 'zorder':0}
        pltargs_c = {'levels':vlev, 'colors':'k', 'linewidths':0.6, 'linestyles':'-', 'zorder':1}
        cArgs = {'orientation':'horizontal', 'location':'top', 'label':'V [m/s]'}
        var_dict = [{'var':V, 'plotType':'contourf', 'plotArgs':pltargs, 'colorArgs':cArgs}]
        var_dict_c = [{'var':V, 'plotType':'contour', 'plotArgs':pltargs_c, 'colorFormatter':None}]    
        pltvert(lat, lev, var_dict, ax=axv, plot_zscale=True, slice_at='zonal mean', xlabel='')
        pltvert(lat, lev, var_dict_c, ax=axv, plot_zscale=False, inverty=False, slice_at='', xlabel='')
        axv.set_ylabel('p  [hPa]')
        axv.set_xlabel('lat  [deg]')
        
        pltargs = {'levels':levels, 'cmap':cmap, 'zorder':0}
        pltargs_c = {'levels':levels, 'colors':'k', 'linewidths':0.6, 'linestyles':'-', 'zorder':1}
        cArgs = {'orientation':'horizontal', 'location':'top', 'label':'OM [Pa/s]'}
        var_dict = [{'var':OM, 'plotType':'contourf', 'plotArgs':pltargs, 'colorArgs':cArgs}]
        var_dict_c = [{'var':OM, 'plotType':'contour', 'plotArgs':pltargs_c, 'colorFormatter':None}]    
        pltvert(lat, lev, var_dict, ax=axom, plot_zscale=True, slice_at='zonal mean', xlabel='')
        pltvert(lat, lev, var_dict_c, ax=axom, plot_zscale=False, inverty=False, slice_at='', xlabel='')
        axom.set_ylabel('p  [hPa]')
        axom.set_xlabel('lat  [deg]')

    #fig.suptitle('{}, {}'.format(title, ttitle), fontsize=14)
    #fig.suptitle('{}, day 150'.format(title, ttitle), fontsize=14)
    plt.tight_layout()
    #plt.savefig('{}/{}_t{}.png'.format(savedest, title, tsnap), dpi=150)
    plt.figure(figu.number)
    plt.savefig('{}/{}_U.png'.format(savedest, title), dpi=200)
    plt.figure(figT.number)
    plt.savefig('{}/{}_T.png'.format(savedest, title), dpi=200)
    #plt.show()
        





if(__name__ == '__main__'):


    gmt = '/glade/u/home/jhollowed/repos/climate_analysis/CLDERA/SAI/analysis/cmaps/cmp_flux.rgb'

    eam_hsw = '/global/cscratch1/sd/jhollo/E3SM/E3SMv2_cases/hsw_validate_cases/'\
              'E3SM_ne16_L72_FIDEAL_30year/run' 
    hsw_6momeans = '{}/E3SM_ne16_L72_FIDEAL_30year.eam.h0.0001-01-01-00000.regrid.2x2.nc'.format(eam_hsw)
    dest = './figs/eam_hsw_circulation'
    circulation_snapshot(hsw_6momeans, [1825, 3650], 'EAM HSW ne16L72, 5 year mean', 
                         dest, inclTracers=False)

    ic = '{}/E3SM_ne16_L72_FIDEAL_30year.eam.i.0006-01-01-00000.nc.regrid.2x2.nc'.format(eam_hsw)
    circulation_snapshot(ic, 0, 'EAM HSW ne16L72, IC from spinup year 6', dest, inclTracers=False)
    


    #dest = '../figs/inithist_by_year'
    #ic_files = sorted(glob.glob('{}/*.i.*regrid*'.format(eam_hsw)))

    #for ic in ic_files:
    #    day = int(ic.split('.i.')[-1].split('-01')[0])
    #    circulation_snapshot(ic, 0, 'EAM HSW ne16L72, year {}'.format(day), dest, inclTracers=False)
    
