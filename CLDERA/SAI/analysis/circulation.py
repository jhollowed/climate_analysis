# eruption_movies.py
# Joe Hollowed
# CLIMATE 589, 4/15/2022
#
# renders circulation snapshots as zonal means for single times, or means of time windows

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


def circulation_snapshot(runf, tsnap, title, ttitle=None, savedest=None, inclTracers=True):
   
    if(isinstance(tsnap, list)):
        ti, tf = tsnap[0], tsnap[1]
        tsnap = slice(tsnap[0], tsnap[1])
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
        if(ttitle is None):
            ttitle = 'mean of days {}-{}'.format(ti, tf)
    else: 
        method = 'nearest'
        mean = 'lon'
        if(ttitle is None):
            ttitle = 'day {}'.format(tsnap)
   
    U   = dat['U'].sel({'time':tsnap}, method=method).mean(mean)
    T   = dat['T'].sel({'time':tsnap}, method=method).mean(mean)
    #V   = dat['V'].sel({'time':tsnap}, method=method).mean(mean)
    #OM  = dat['OMEGA'].sel({'time':tsnap}, method=method).mean(mean)
    if inclTracers:
        C   = dat['SAI_SO2'].sel({'time':tsnap}, method=method).sel({'lon':120.35}, method='nearest')
    lat = dat['lat']
    lev = dat['lev']

    # if tsnap was a slice, take time mean
    if len(U.shape) == 3:
        U = U.mean('time')
        T = T.mean('time')
        #V = V.mean('time')
        #OM = OM.mean('time')

    # ---------- plot ----------
    data_crs = ccrs.PlateCarree()
    Tlev = np.linspace(150, 300, 11)
    uleva = np.arange(0, np.max(U), 15)
    ulev = np.hstack([-uleva[1:][::-1], uleva])
    
    #cmap = claut.ncar_rgb_to_cmap(gmt)
    cmap_u = mpl.cm.rainbow
    cmap_T = mpl.cm.YlOrRd
    #cmap = mpl.cm.Spectral_r
        
    
    print('\n\n=============== {}'.format(title)) 

    fig = plt.figure(figsize=(6, 5))
    axu = fig.add_subplot(111)
    #axu = fig.add_subplot(221)
    #axv = fig.add_subplot(222)
    #axom = fig.add_subplot(223)
    
    figT = plt.figure(figsize=(6,5))
    axT = figT.add_subplot(111)

    # ----- vertical
    
    pltargs = {'levels':ulev, 'cmap':cmap_u, 'zorder':0}
    pltargs_c = {'levels':ulev, 'colors':'k', 'linewidths':0.6, 'zorder':1}
    cArgs = {'orientation':'horizontal', 'location':'top', 'label':'U [m/s]'}
    cArgs_c = {'fmt':'%d'}
    var_dict = [{'var':U, 'plotType':'contourf', 'plotArgs':pltargs, 'colorArgs':cArgs}]
    var_dict_c = [{'var':U, 'plotType':'contour', 'plotArgs':pltargs_c, 'colorArgs':cArgs_c}]    
    cf = pltvert(lat, lev, var_dict, ax=axu, plot_zscale=True, slice_at='zonal mean', xlabel='')
    pltvert(lat, lev, var_dict_c, ax=axu, plot_zscale=False, inverty=False, slice_at='', xlabel='')
    cf[0].set_ticks(ulev)
    
    pltargs = {'levels':Tlev, 'cmap':cmap_T, 'zorder':0}
    pltargs_c = {'levels':Tlev, 'colors':'k', 'linewidths':0.6, 'zorder':1}
    cArgs = {'orientation':'horizontal', 'location':'top', 'label':'T [K]'}
    cArgs_c = {'fmt':'%d'}
    var_dict = [{'var':T, 'plotType':'contourf', 'plotArgs':pltargs, 'colorArgs':cArgs}]
    var_dict_c = [{'var':T, 'plotType':'contour', 'plotArgs':pltargs_c, 'colorArgs':cArgs_c}]    
    cf = pltvert(lat, lev, var_dict, ax=axT, plot_zscale=True, slice_at='zonal mean', xlabel='')
    pltvert(lat, lev, var_dict_c, ax=axT, plot_zscale=False, inverty=False, slice_at='', xlabel='')
    cf[0].set_ticks(Tlev)

    if inclTracers:
        #pltargs_c = {'levels':[-5, -4], 'colors':'m', 'linewidths':1.6, "linestyles":'-'}
        pltargs_c = {'levels':[-7, -5.5], 'colors':'m', 'linewidths':1.6, "linestyles":'-'}
        cargs = {'fmt':'%.0f', 'manual':[(18.8, 201), (14.89, 26.2)]} 
        var_dict_c = [{'var':np.log10(C), 'plotType':'contour', 'plotArgs':pltargs_c, 'colorArgs':cargs}]
        pltvert(lat, lev, var_dict_c, ax=axu, plot_zscale=False,inverty=False,slice_at='',xlabel='')
        axu.plot([0,0],[100,100],'-m', label='log10(concentration)')
        axu.legend(loc='lower right', fancybox=False, framealpha=1)

    axu.set_ylabel('p  [hPa]')
    axu.set_xlabel('latitude')
    axT.set_ylabel('p  [hPa]')
    axT.set_xlabel('latitude')
   
    if(0): # currentlty disable V, OMEGA plots
        pltargs = {'levels':vlev, 'cmap':cmap, 'zorder':0}
        pltargs_c = {'levels':vlev, 'colors':'k', 'linewidths':0.6, 'linestyles':'-', 'zorder':1}
        cArgs = {'orientation':'horizontal', 'location':'top', 'label':'V [m/s]'}
        cArgs_c = {'fmt':'%.0f'}
        var_dict = [{'var':V, 'plotType':'contourf', 'plotArgs':pltargs, 'colorArgs':cArgs}]
        var_dict_c = [{'var':V, 'plotType':'contour', 'plotArgs':pltargs_c, 'colorArgs':cArgs_c}]    
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
    plt.tight_layout()
    if(savedest is not None):
        fig.savefig('{}/{}_t{}_U.png'.format(savedest, title, tsnap), dpi=150)
        figT.savefig('{}/{}_t{}_T.png'.format(savedest, title, tsnap), dpi=150)
    else: 
        plt.show()
        





if(__name__ == '__main__'):

    #gmt = '/glade/u/home/jhollowed/repos/climate_analysis/CLDERA/SAI/analysis/cmaps/amwg256.rgb'
    #gmt = '/glade/u/home/jhollowed/repos/climate_analysis/CLDERA/SAI/analysis/cmaps/GMT_no_green.rgb'
    #gmt = '/glade/u/home/jhollowed/repos/climate_analysis/CLDERA/SAI/analysis/cmaps/GHRSST_anomaly.rgb'
    gmt = '/glade/u/home/jhollowed/repos/climate_analysis/CLDERA/SAI/analysis/cmaps/cmp_flux.rgb'
    #gmt = '/glade/u/home/jhollowed/repos/climate_analysis/CLDERA/SAI/analysis/cmaps/WhBlGrYeRe.rgb'
    #whs = '/glade/scratch/jhollowed/CAM/cases/sai_runs/SE_ne16L72_whs_saiv2_fix0_tau0_qsplit1/'\
    #      'run/SE_ne16L72_whs_saiv2_fix0_tau0_qsplit1.cam.h0.0001-01-01-00000.nc'
   
    #whs_spinup = '/glade/scratch/jhollowed/CAM/cases/sai_runs/SE_ne16L72_whs_sai_spinup/run/SE_ne16L72_whs_sai_spinup.cam.h0.0001-10-28-00000.nc'
    #whs = '/glade/scratch/jhollowed/CAM/cases/sai_runs/SE_ne16L72_whs_sai_fix0_tau0_nsplit1_nodiff0/'\
    #      'run/SE_ne16L72_whs_sai_fix0_tau0_nsplit1_nodiff0.cam.h0.0001-01-01-00000.nc'
    #whs_massfix = '/glade/scratch/jhollowed/CAM/cases/sai_runs/SE_ne16L72_whs_saiv2_fix1_tau0_qsplit1/'\
    #              'SE_ne16L72_whs_saiv2_fix1_tau0_qsplit1.cam.h0.0001-01-01-00000.nc'
    #whsdest = '/glade/u/home/jhollowed/repos/climate_analysis/CLDERA/SAI/analysis/figs/'
    #whsgdest = '/glade/u/home/jhollowed/repos/climate_analysis/CLDERA/SAI/analysis/figs2/whsg'
    
    
    #amip = '/glade/scratch/jhollowed/CAM/cases/sai_runs/E3SM_AMIP_ne30_L72_SAI_juneclimo/'\
    #       'E3SM_case_ne30_L72_SAI_amip_juneclimo.eam.h0.0001-01-01-00000.regird.2x2.nc'
    #amipdest = '/glade/u/home/jhollowed/repos/climate_analysis/CLDERA/SAI/analysis/figs/amip'
    
    #circulation_snapshot(whs, 1, 'SE ne30L72, CAM HSW', whsdest)
    #circulation_snapshot(whs, [0,60], 'SE ne30L72, CAM HSW', whsdest)
    #circulation_snapshot(amip, 0.5, 'SE ne30L72, EAM AMIP', amipdest)
    #circulation_snapshot(amipJune, 0, 'SE ne30L72, EAM AMIP', amipdest, inclTracers=False)
    #circulation_snapshot(whs_spinup, [0,60], 'SE ne30L72, CAM HSW', whsdest, inclTracers=False)

    # SPINUP FOR IMPROVED WHS TOP RUNS
    whs_improved_initfiles = sorted(glob.glob('/global/cscratch1/sd/jhollo/E3SM/E3SMv2_cases/'\
                                      'hsw_validate_cases/E3SM_ne16_L72_FIDEAL_10year_spinup_withNewTeq/'\
                                      'run/*eam.i*regrid.2x2.nc'))[::-1]
    dest = '/global/homes/j/jhollo/repos/climate_analysis/CLDERA/SAI/analysis/figs/updated_whs_apole'
    for i in range(len(whs_improved_initfiles)):
        initfile = whs_improved_initfiles[i]
        year = int(initfile.split('-')[0].split('.')[-1])
        if(year != 5): continue
        circulation_snapshot(initfile, 0, 'SE ne16L72, EAM HSW', 
                             ttitle='year {}'.format(year), savedest=dest, inclTracers=False)
    
