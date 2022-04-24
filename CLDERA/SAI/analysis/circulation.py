# eruption_movies.py
# Joe Hollowed
# CLIMATE 589, 4/15/2022
#
# renders video frames of the eruption in horizontal cross, vertical cross, 
# and AzimuthalEquidistant projection

import numpy as np
import xarray as xr
import matplotlib as mpl
import cartopy.crs as ccrs
import artist_utils as claut
import climate_toolbox as ctb
import matplotlib.pyplot as plt
from pdb import set_trace as st
from matplotlib.offsetbox import AnchoredText
from climate_artist import vertical_slice as pltvert
from climate_artist import horizontal_slice as plthor

# ============================================================


def circulation_snapshot(runf, tsnap, title, savedest):
   
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
    run = xr.open_dataset(runf)
    td = ctb.time2day(run['time'])
    run = run.assign_coords(time=td)
     
    # get vars at time snapshot/slice, take zonal means
    U   = run['U'].sel({'time':tsnap}, method='nearest').mean('lon')
    V   = run['V'].sel({'time':tsnap}, method='nearest').mean('lon')
    OM  = run['OMEGA'].sel({'time':tsnap}, method='nearest').mean('lon')
    lat = run['lat']
    lev = run['lev']

    # if tsnap was a slice, take time mean
    if len(U.shape) == 3:
        U = U.mean('time')
        V = V.mean('time')
        OM = OM.mean('time')

    # ---------- plot ----------
    data_crs = ccrs.PlateCarree()
    #levels = np.linspace(minc, maxc, clevels)
    levels=10
    uleva = np.arange(0, np.max(U), 15)
    ulev = np.hstack([-uleva[1:][::-1], uleva])
    vlev = np.linspace(-3, 2, 10)
    #cmap = claut.ncar_rgb_to_cmap(gmt)
    cmap = mpl.cm.rainbow
    #cmap = mpl.cm.Spectral_r
        
    print('\n\n=============== {}'.format(title)) 

    #fig = plt.figure(figsize=(8,10))
    fig = plt.figure()
    axu = fig.add_subplot(111)
    #axu = fig.add_subplot(221)
    #axv = fig.add_subplot(222)
    #axom = fig.add_subplot(223)

    # ----- vertical
    
    pltargs = {'levels':ulev, 'cmap':cmap, 'zorder':0}
    #pltargs_c = {'levels':ulev, 'colors':'k', 'linewidths':0.6, 'linestyles':'-', 'zorder':1}
    pltargs_c = {'levels':ulev, 'colors':'k', 'linewidths':0.6, 'zorder':1}
    cArgs = {'orientation':'horizontal', 'location':'top', 'label':'U [m/s]'}
    var_dict = [{'var':U, 'plotType':'contourf', 'plotArgs':pltargs, 'colorArgs':cArgs}]
    var_dict_c = [{'var':U, 'plotType':'contour', 'plotArgs':pltargs_c, 'colorFormatter':None}]    
    pltvert(lat, lev, var_dict, ax=axu, plot_zscale=True, slice_at='zonal mean', xlabel='')
    pltvert(lat, lev, var_dict_c, ax=axu, plot_zscale=False, inverty=False, slice_at='', xlabel='')
    axu.set_ylabel('p  [hPa]')
    axu.set_xlabel('lat  [deg]')
    
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
    
    fig.suptitle('{}, day {}'.format(title, tsnap), fontsize=14)
    plt.tight_layout()
    plt.show()
    plt.savefig('{}/{}_t{}.png'.format(savedest, title, tsnap), dpi=150)
        





if(__name__ == '__main__'):

    #gmt = '/glade/u/home/jhollowed/repos/climate_analysis/CLDERA/SAI/analysis/cmaps/amwg256.rgb'
    #gmt = '/glade/u/home/jhollowed/repos/climate_analysis/CLDERA/SAI/analysis/cmaps/GMT_no_green.rgb'
    #gmt = '/glade/u/home/jhollowed/repos/climate_analysis/CLDERA/SAI/analysis/cmaps/GHRSST_anomaly.rgb'
    gmt = '/glade/u/home/jhollowed/repos/climate_analysis/CLDERA/SAI/analysis/cmaps/cmp_flux.rgb'
    #gmt = '/glade/u/home/jhollowed/repos/climate_analysis/CLDERA/SAI/analysis/cmaps/WhBlGrYeRe.rgb'
    whs = '/glade/scratch/jhollowed/CAM/cases/sai_runs/SE_ne16L72_whs_saiv2_fix0_tau0_qsplit1/'\
          'run/SE_ne16L72_whs_saiv2_fix0_tau0_qsplit1.cam.h0.0001-01-01-00000.nc'
    whsdest = '/glade/u/home/jhollowed/repos/climate_analysis/CLDERA/SAI/analysis/figs2/whs'
    whsgdest = '/glade/u/home/jhollowed/repos/climate_analysis/CLDERA/SAI/analysis/figs2/whsg'
    whs_massfix = '/glade/scratch/jhollowed/CAM/cases/sai_runs/SE_ne16L72_whs_saiv2_fix1_tau0_qsplit1/'\
                  'SE_ne16L72_whs_saiv2_fix1_tau0_qsplit1.cam.h0.0001-01-01-00000.nc'
    amip = '/glade/scratch/jhollowed/CAM/cases/sai_runs/E3SM_AMIP_ne30_L72_SAI/'\
           'AMIPcase_ne30_L72_SAI.eam.h0.0001-01-01-00000.nc.regrid.2x2.nc'
    amipdest = '/glade/u/home/jhollowed/repos/climate_analysis/CLDERA/SAI/analysis/figs2/amip'
    
    #circulation_snapshot(whs, 60, 'SE ne30L72, CAM HSW', whsdest)
    circulation_snapshot(amip, 0, 'SE ne30L72, EAM AMIP', amipdest)
    
