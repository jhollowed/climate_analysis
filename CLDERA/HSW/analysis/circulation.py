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
from matplotlib import colors
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
from climate_artist import vertical_slice as pltvert
from climate_artist import horizontal_slice as plthor

# ============================================================


def circulation_snapshot(runf, tsnap, title, ttitle=None, maketitle=True, 
                         savedest=None, inclTracers=True, auto_levels=False):
   
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
    if(auto_levels):
        Tlev = 9
    else:
        Tlev = np.linspace(180, 300, 9).astype(int)
    
    #uleva = np.arange(0, np.max(U), 15)
    #ulev = np.hstack([-uleva[1:][::-1], uleva])
    if(auto_levels):
        ulev = 9
    else:
        ulev = np.hstack([np.linspace(-30, 0, 6), np.linspace(15, 75, 5)]).astype(int)
    divnorm=colors.TwoSlopeNorm(vcenter=0.)
    
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

    if(tismean):
        annot = 'mean from year 5-10'
    else:
        annot = 'zonal mean at year 5'

    # ----- vertical
    pltargs = {'levels':ulev, 'cmap':cmap_u, 'zorder':0, 'norm':divnorm}
    pltargs_c = {'levels':ulev, 'colors':'k', 'linewidths':0.6, 'zorder':1, 'norm':divnorm}
    cArgs = {'orientation':'horizontal', 'location':'top', 'label':'U [m/s]','aspect':30,'format':'%d'}
    cArgs_c = {'fmt':'%d'}
    var_dict = [{'var':U, 'plotType':'contourf', 'plotArgs':pltargs, 'colorArgs':cArgs}]
    var_dict_c = [{'var':U, 'plotType':'contour', 'plotArgs':pltargs_c, 'colorArgs':cArgs_c}]    
    cf = pltvert(lat, lev, var_dict, ax=axu, plot_zscale=True,
                 annotation=annot, annotation_alpha=0.85, annotation_loc='lower right')
    pltvert(lat, lev, var_dict_c, ax=axu, plot_zscale=False, inverty=False, annotation='', xlabel='')
    if(not(auto_levels)): cf[0].set_ticks(ulev)
    
    pltargs = {'levels':Tlev, 'cmap':cmap_T, 'zorder':0}
    pltargs_c = {'levels':Tlev, 'colors':'k', 'linewidths':0.6, 'zorder':1}
    cArgs = {'orientation':'horizontal', 'location':'top', 'label':'T [K]', 'aspect':30, 'format':'%d'}
    cArgs_c = {'fmt':'%d'}
    var_dict = [{'var':T, 'plotType':'contourf', 'plotArgs':pltargs, 'colorArgs':cArgs}]
    var_dict_c = [{'var':T, 'plotType':'contour', 'plotArgs':pltargs_c, 'colorArgs':cArgs_c}]    
    cf = pltvert(lat, lev, var_dict, ax=axT, plot_zscale=True, 
                 annotation=annot, annotation_alpha=0.85, annotation_loc='lower right')
    pltvert(lat, lev, var_dict_c, ax=axT, plot_zscale=False, inverty=False, annotation='', xlabel='')
    if(not(auto_levels)): cf[0].set_ticks(Tlev)

    if inclTracers:
        #pltargs_c = {'levels':[-5, -4], 'colors':'m', 'linewidths':1.6, "linestyles":'-'}
        pltargs_c = {'levels':[-7, -5.5], 'colors':'m', 'linewidths':1.6, "linestyles":'-'}
        cargs = {'fmt':'%.0f', 'manual':[(18.8, 201), (14.89, 26.2)]} 
        var_dict_c = [{'var':np.log10(C), 'plotType':'contour', 'plotArgs':pltargs_c, 'colorArgs':cargs}]
        pltvert(lat, lev, var_dict_c, ax=axu, plot_zscale=False,inverty=False,annotation='',xlabel='')
        axu.plot([0,0],[100,100],'-m', label='log10(concentration)')
        axu.legend(loc='lower right', fancybox=False, framealpha=1)

    axu.set_ylabel('p  [hPa]')
    axu.set_xlabel('latitude')
    axT.set_ylabel('p  [hPa]')
    axT.set_xlabel('latitude')
   
    if(maketitle):
        fig.suptitle('{}, {}'.format(title, ttitle), fontsize=14)
        figT.suptitle('{}, {}'.format(title, ttitle), fontsize=14)
    fig.tight_layout()
    figT.tight_layout()
    if(savedest is not None):
        fig.savefig('{}/{}_t{}_U.png'.format(savedest, title, tsnap), dpi=150)
        figT.savefig('{}/{}_t{}_T.png'.format(savedest, title, tsnap), dpi=150)
        plt.show()
    else: 
        plt.show()
        





if(__name__ == '__main__'):

    #gmt = '/glade/u/home/jhollowed/repos/climate_analysis/CLDERA/SAI/analysis/cmaps/amwg256.rgb'
    #gmt = '/glade/u/home/jhollowed/repos/climate_analysis/CLDERA/SAI/analysis/cmaps/GMT_no_green.rgb'
    #gmt = '/glade/u/home/jhollowed/repos/climate_analysis/CLDERA/SAI/analysis/cmaps/GHRSST_anomaly.rgb'
    gmt = '/glade/u/home/jhollowed/repos/climate_analysis/CLDERA/SAI/analysis/cmaps/cmp_flux.rgb'
    #gmt = '/glade/u/home/jhollowed/repos/climate_analysis/CLDERA/SAI/analysis/cmaps/WhBlGrYeRe.rgb'

    # SPINUP FOR REBASED RUNS
    if(0):
        initfiles = sorted(glob.glob('/global/cscratch1/sd/jhollo/E3SM/E3SMv2_cases/'\
                                     'hsw_cases/E3SM_ne16_L72_FIDEAL_10year_'\
                                     'spinup/run/*eam.i*regrid.2x2.nc'))[::-1]
        dest='./figs/inithist_by_year_rebase'
        
        for i in range(len(initfiles)):
            initfile = initfiles[i]
            year = int(initfile.split('-')[0].split('.')[-1])
            circulation_snapshot(initfile, 0, 'SE_ne16pg2_L72_EAM_HSW_year{}'.format(year), 
                                 maketitle=False, 
                                 ttitle='year{}'.format(year), savedest=dest, inclTracers=False)
    
    # TIME MEAN OF SPINUP FOR REBASED FIDEAL
    if(0):
        avg_file = '/global/cscratch1/sd/jhollo/E3SM/E3SMv2_cases/hsw_cases/'\
                   'E3SM_ne16_L72_FIDEAL_10year_spinup/run/'\
                   'E3SM_ne16_L72_FIDEAL_10year_spinup.eam.h0.0001-01-01-00000'\
                   '.nc.regrid.2x2.nc'
        dest='./figs/eam_hsw_ciruclation_rebase'
        
        circulation_snapshot(avg_file, [int(360*3), int(360*7)], 
                             'SE_ne16pg2_L72_EAM_HSW_5yearAvg', maketitle=False, 
                             ttitle='', inclTracers=False, savedest=dest)
    
    # DIFFERENCE IN CIRUCLATION FROM PRE-REBASE STATE
    if(1):
        initfiles = sorted(glob.glob('/global/cscratch1/sd/jhollo/E3SM/E3SMv2_cases/'\
                                     'hsw_cases/E3SM_ne16_L72_FIDEAL_10year_'\
                                     'spinup/run/*eam.i*PREREBASE_DIFF'))[::-1]
        dest='./figs/inithist_by_year_rebase'
        
        for i in range(len(initfiles)):
            initfile = initfiles[i]
            year = int(initfile.split('-')[0].split('.')[-1])
            circulation_snapshot(initfile,0,'SE_ne16pg2_L72_EAM_HSW_year{}_PREREBASE_DIFF'.format(year), 
                                 maketitle=False, 
                                 ttitle='year{}'.format(year), savedest=dest, inclTracers=False,
                                 auto_levels=True)
        
