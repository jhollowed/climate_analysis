# eruption_movies.py
# Joe Hollowed
# CLIMATE 589, 4/15/2022
#
# renders video frames of the eruption in horizontal cross, vertical cross, 
# and AzimuthalEquidistant projection

import pdb
import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import climate_toolbox as ctb
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
from climate_artist import horizontal_slice as plthor
from climate_artist import vertical_slice as pltvert

# ============================================================


def animate_eruption(runf, title, savedest, tracer='SO2', globe=True):
   
    # params
    lat0 = 15.15
    lon0 = 120.35
    dlat = 0.5
    psel = 45
    minc = -10
    maxc = -4
    clevels = 8

    # read data
    run = xr.open_dataset(runf)
    c = run['SAI_{}'.format(tracer)]
    
    # for horizontal slices
    lat = c['lat']
    lon = c['lon']
    chor = np.log10(c.sel({'lev':psel}, method='nearest'))
    chor = np.ma.masked_array(chor, np.logical_or(chor == -float('inf'), 
                                                  chor < minc))
    
    # for vertical slice
    lev = c['lev']
    latsel = slice(lat0-dlat, lat0+dlat)
    cvert = np.log10(c.sel({'lat':latsel}).mean('lat'))
    cvert = np.ma.masked_array(cvert, np.logical_or(cvert == -float('inf'), 
                                                    cvert < minc))

    # get time in number of days
    td = ctb.time2day(run['time'])
    
    # ---------- plot ----------
    levels = np.linspace(minc, maxc, clevels)
    cmap = ctb.ncar_rgb_to_cmap(gmt)
    data_crs = ccrs.PlateCarree()
        
    print('\n\n=============== {}'.format(title)) 

    for k in range(len(td)):

        if(k%2 != 0): continue

        print('--------- {}'.format(k)) 
        if(globe):
            fig = plt.figure(figsize=(10,6))
            spec = fig.add_gridspec(2, 2)
            ax1 = fig.add_subplot(spec[0,0], projection=ccrs.PlateCarree(lon0)) # horizontal
            ax2 = fig.add_subplot(spec[1,0]) # vertical
            ax3 = fig.add_subplot(spec[:,1], projection=ccrs.AzimuthalEquidistant(lon0, 90)) # pole
        else:
            fig = plt.figure(figsize=(8,10))
            ax1 = fig.add_subplot(211, projection=ccrs.PlateCarree(lon0)) # horizontal
            ax2 = fig.add_subplot(212) # vertical

        # ----- horizontal
        pltargs = {'levels':levels, 'cmap':cmap}
        pltargs_c = {'levels':levels, 'colors':'k', 'linewidths':0.6, 'linestyles':'-', }
        var_dict = [{'var':chor[k], 'plotType':'contourf', 'plotArgs':pltargs, 
                     'colorFormatter':None}]
        var_dict_c = [{'var':chor[k], 'plotType':'contour', 'plotArgs':pltargs_c, 'colorFormatter':None}]
        plthor(lon, lat, var_dict, ax=ax1)
        plthor(lon, lat, var_dict_c, ax=ax1)
        text_box = AnchoredText('p=45 hPa', frameon=True, loc='lower left', pad=0.5)
        text_box.set_zorder(100)
        plt.setp(text_box.patch, facecolor='white', alpha=1)
        ax1.add_artist(text_box)
        
        
        # ----- vertical
        
        # recenter on lon0
        lonmax = lon0+180
        lonroll = np.searchsorted(lon, lonmax)
        
        lonv = lon - lon0
        lonv[lonv > 180] = lonv[lonv > 180] - 360
        lonv = np.roll(lonv, -lonroll)
        lonlabs = np.roll(lon, -lonroll)

        
        pltargs = {'levels':levels, 'cmap':cmap, 'zorder':0}
        pltargs_c = {'levels':levels, 'colors':'k', 'linewidths':0.6, 'linestyles':'-', 'zorder':1}
        cArgs = {'orientation':'horizontal', 'location':'bottom', 'label':'log10(c_SO2 [kg/kg])'}
        var_dict = [{'var':cvert[k], 'plotType':'contourf', 'plotArgs':pltargs, 'colorArgs':cArgs}]
        var_dict_c = [{'var':cvert[k], 'plotType':'contour', 'plotArgs':pltargs_c, 'colorFormatter':None}]    


        pltvert(lonv, lev, var_dict, ax=ax2, plot_zscale=True)
        pltvert(lonv, lev, var_dict_c, ax=ax2, plot_zscale=False, inverty=False)
        text_box = AnchoredText('lat=15.15 deg', frameon=True, loc='lower left', pad=0.5)
        text_box.set_zorder(100)
        plt.setp(text_box.patch, facecolor='white', alpha=1)
        ax2.add_artist(text_box)
        ax2.set_ylabel('p  [hPa]')
        #ax2.set_xticks(ax1.get_xticks())
        ax2.set_xticks([-120, -60, 0, 60, 120])
        ax2.set_xlim(ax1.get_xlim())
        #pdb.set_trace()
        ax2.set_xticklabels([])
        ax2.xaxis.set_tick_params(direction='in', which='both')
        ax2.xaxis.set_ticks_position('both')


        if(globe):
            # ----- horizontal
            pltargs = {'levels':levels, 'cmap':cmap}
            pltargs_c = {'levels':levels, 'colors':'k', 'linewidths':0.6}
            var_dict = [{'var':chor[k], 'plotType':'contourf', 'plotArgs':pltargs, 'colorFormatter':None}]
            var_dict_c = [{'var':chor[k], 'plotType':'contour', 'plotArgs':pltargs_c}]
            gridlinesArgs = {'draw_labels':False}
            plthor(lon, lat, var_dict, ax=ax3, gridlinesArgs=gridlinesArgs)
            #plthor(lat, lon, var_dict_c, ax=ax1)
        
        fig.suptitle('{}, day {}'.format(title, round(td[k]+0.)), fontsize=14)
        plt.tight_layout()
        plt.savefig('{}/{:03d}.png'.format(savedest, k), dpi=150)
        





if(__name__ == '__main__'):

    gmt = '/glade/u/home/jhollowed/repos/climate_analysis/CLDERA/SAI/analysis/GMT_no_green.rgb'
    whs = '/glade/scratch/jhollowed/CAM/cases/sai_runs/SE_ne16L72_whs_saiv2_fix0_tau0_qsplit1/'\
          'run/SE_ne16L72_whs_saiv2_fix0_tau0_qsplit1.cam.h0.0001-01-01-00000.nc'
    whsdest = '/glade/u/home/jhollowed/repos/climate_analysis/CLDERA/SAI/analysis/figs/whs'
    whsgdest = '/glade/u/home/jhollowed/repos/climate_analysis/CLDERA/SAI/analysis/figs2/whsg'
    whs_massfix = '/glade/scratch/jhollowed/CAM/cases/sai_runs/SE_ne16L72_whs_saiv2_fix1_tau0_qsplit1/'\
                  'SE_ne16L72_whs_saiv2_fix1_tau0_qsplit1.cam.h0.0001-01-01-00000.nc'
    amip = '/glade/scratch/jhollowed/CAM/cases/sai_runs/E3SM_AMIP_ne30_L72_SAI/'\
           'AMIPcase_ne30_L72_SAI.eam.h0.0001-01-01-00000.nc.regrid.2x2.nc'
    amipdest = '/glade/u/home/jhollowed/repos/climate_analysis/CLDERA/SAI/analysis/figs/amip'
    
    animate_eruption(whs, 'SE ne30L72, CAM WHS', whsgdest, globe=True)
    #animate_eruption(whs, 'SE ne30L72, CAM WHS', whsdest, globe=False)
    #animate_eruption(amip, 'SE ne30L72, EAM AMIP', amipdest, globe=False)
    
