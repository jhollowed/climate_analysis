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


def animate_eruption(runf, title, savedest, tracer='SO2', globe=True, demo=False):
   
    # params
    lat0 = 15.15
    lon0 = 120.35
    lon_center = 0
    dlat = 0.5
    psel = 45
    minc = -12
    clipc = -15
    maxc = -4
    clevels = 9

    # read data
    run = xr.open_dataset(runf)
    c = run['SAI_{}'.format(tracer)]
   
    # for horizontal slices
    lat = c['lat']
    lon = c['lon']
    chor = c.sel({'lev':psel}, method='nearest')           # replace zeros with tiny value, for log
    mask = np.logical_or(chor == 0, chor < 10**clipc)
    chor = np.ma.masked_array(chor, mask).filled(10**clipc) 
    chor = np.log10(chor)
    
    # for vertical slice
    lev = c['lev']
    latsel = slice(lat0-dlat, lat0+dlat)
    #cvert = np.log10(c.sel({'lat':latsel}).mean('lat'))
    #mask = np.logical_or(cvert == -float('inf'), cvert < minc)
    #cvert = np.ma.masked_array(cvert, mask).filled(minc)
    cvert = c.sel({'lat':latsel}).mean('lat')              # replace zeros with tiny value, for log
    mask = np.logical_or(cvert == 0, cvert < 10**clipc)
    cvert = np.ma.masked_array(cvert, mask).filled(10**clipc) 
    cvert = np.log10(cvert)

    # get time in number of days
    td = ctb.time2day(run['time'])
    
    # ---------- plot ----------
    levels = np.linspace(minc, maxc, clevels)
    cmap = claut.ncar_rgb_to_cmap(gmt)
    #cmap = claut.ncar_rgb_to_cmap(gmt, norm=True)
    #cmap = mpl.cm.YlGnBu
    data_crs = ccrs.PlateCarree()
        
    print('\n\n=============== {}'.format(title)) 

    for k in range(len(td)):

        if(demo): k+= 10
        #if(k%2 != 0): continue     # every other time sample
        if(td[k] > 30): continue    # stop after day 40

        print('--------- {}'.format(k)) 
        if(globe):
            fig = plt.figure(figsize=(10,6))
            spec = fig.add_gridspec(2, 2)
            ax1 = fig.add_subplot(spec[1,0], projection=ccrs.PlateCarree(lon_center)) # horizontal
            ax2 = fig.add_subplot(spec[0,0]) # vertical
            ax3 = fig.add_subplot(spec[:,1], projection=ccrs.AzimuthalEquidistant(lon0, 90)) # pole
        else:
            fig = plt.figure(figsize=(8,10))
            ax1 = fig.add_subplot(212, projection=ccrs.PlateCarree(lon_center)) # horizontal
            ax2 = fig.add_subplot(211) # vertical

        # ----- horizontal
        pltargs = {'levels':levels, 'cmap':cmap}
        pltargs_c = {'levels':levels, 'colors':'k', 'linewidths':0.6, 'linestyles':'-', }
        var_dict = [{'var':chor[k], 'plotType':'contourf', 'plotArgs':pltargs, 
                     'colorFormatter':None}]
        var_dict_c = [{'var':chor[k], 'plotType':'contour', 'plotArgs':pltargs_c, 'colorFormatter':None}]
        plthor(lon, lat, var_dict, ax=ax1, slice_at='p=45 hPa')
        plthor(lon, lat, var_dict_c, ax=ax1, slice_at='')
        
        
        # ----- vertical
        
        pltargs = {'levels':levels, 'cmap':cmap, 'zorder':0}
        pltargs_c = {'levels':levels, 'colors':'k', 'linewidths':0.6, 'linestyles':'-', 'zorder':1}
        cArgs = {'orientation':'horizontal', 'location':'top', 'label':'log10(concentration)', 'format':'%.1f'}
        var_dict = [{'var':cvert[k], 'plotType':'contourf', 'plotArgs':pltargs, 'colorArgs':cArgs}]
        var_dict_c = [{'var':cvert[k], 'plotType':'contour', 'plotArgs':pltargs_c, \
                       'colorFormatter':None}]    
        pltvert(lon, lev, var_dict, ax=ax2, plot_zscale=True, center_x=lon_center, slice_at='lat=15.15 deg', xlabel='', slice_at_loc='upper left')
        pltvert(lon, lev, var_dict_c, ax=ax2, plot_zscale=False, inverty=False, center_x=lon_center, slice_at='', xlabel='')
        ax2.set_ylabel('p  [hPa]')
        
        #ax2.set_xticks([0, 60, 120, 180, 240])
        #ax2.set_xlim(np.array(ax1.get_xlim()) - 120)
        ax2.set_xticks(np.array([-180, -120, -60, 0, 60, 120, 180]))
        #ax2.set_xlim(np.array(ax1.get_xlim())+120)
        #ax2.set_xticklabels([])
        #ax2.xaxis.set_tick_params(direction='in', which='both')
        #ax2.xaxis.set_ticks_position('both')


        if(globe):
            # ----- horizontal
            pltargs = {'levels':levels, 'cmap':cmap}
            pltargs_c = {'levels':levels, 'colors':'k', 'linewidths':0.6}
            var_dict = [{'var':chor[k], 'plotType':'contourf', 'plotArgs':pltargs, 'colorFormatter':None}]
            var_dict_c = [{'var':chor[k], 'plotType':'contour', 'plotArgs':pltargs_c}]
            gridlinesArgs = {'draw_labels':False}
            plthor(lon, lat, var_dict, ax=ax3, gridlinesArgs=gridlinesArgs, slice_at='p=45 hPa')
            #plthor(lat, lon, var_dict_c, ax=ax1)
        
        fig.suptitle('{}, day {}'.format(title, round(td[k]+0.)), fontsize=14)
        ax2.set_aspect('auto')
        plt.subplots_adjust(hspace=0)

        if(demo):
            plt.show()
        else:
            plt.savefig('{}/{:03d}.png'.format(savedest, k), dpi=150)
        





if(__name__ == '__main__'):

    gmt = '/glade/u/home/jhollowed/repos/climate_analysis/CLDERA/SAI/analysis/cmaps/GMT_no_green.rgb'
    #gmt = '/glade/u/home/jhollowed/repos/climate_analysis/CLDERA/SAI/analysis/cmaps/WhBlGrYeRe.rgb'
    whs = '/glade/scratch/jhollowed/CAM/cases/sai_runs/SE_ne16L72_whs_sai_fix0_tau0_nsplit1_nodiff0/'\
          'run/SE_ne16L72_whs_sai_fix0_tau0_nsplit1_nodiff0.cam.h0.0001-01-01-00000.nc'
    whsdest = '/glade/u/home/jhollowed/repos/climate_analysis/CLDERA/SAI/analysis/figs/whsg_symm_init_ash'
    #whsgdest = '/glade/u/home/jhollowed/repos/climate_analysis/CLDERA/SAI/analysis/figs3/whsg'
    #whs_massfix = '/glade/scratch/jhollowed/CAM/cases/sai_runs/SE_ne16L72_whs_saiv2_fix1_tau0_qsplit1/'\
    #              'SE_ne16L72_whs_saiv2_fix1_tau0_qsplit1.cam.h0.0001-01-01-00000.nc'
    amip = '/glade/scratch/jhollowed/CAM/cases/sai_runs/E3SM_AMIP_ne30_L72_SAI_juneclimo/'\
           'E3SM_case_ne30_L72_SAI_amip_juneclimo.eam.h0.0001-01-01-00000.regird.2x2.nc'
    amipdest = '/glade/u/home/jhollowed/repos/climate_analysis/CLDERA/SAI/analysis/figs/amip_ash'
    
    #animate_eruption(whs, 'SE ne30L72, CAM HSW', whsgdest, globe=True)
    #animate_eruption(whs, 'SE ne30L72, CAM WHS, SO2', whsdest, globe=True, demo=False)
    #animate_eruption(amip, 'SE ne30L72, EAM AMIP, SO2', amipdest, globe=True, demo=False, tracer='SO2')
    animate_eruption(whs, 'SE ne30L72, CaM HSW, ash', whsdest, globe=True, demo=False, tracer='ASH')
    
    
