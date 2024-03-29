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
from climate_artist import horizontal_slice as plthor

# ============================================================


def circulation_snapshot(runf, title, savedest, inclTracers=True, vertarr=False):
   

    # read data, transform time coordinate to ndays
    dat = xr.open_dataset(runf)
    th = dat['nsteph']*1800/60/60
    dat = dat.assign_coords(time=th)
     
    T050   = dat['T050']
    T850   = dat['T850']
    lat = dat['lat']
    lon = dat['lon']

    if(not vertarr):
        # horizontally arrange plots
        orientation='vertical'
        location='right'
        ax1n = 121
        ax2n = 122
        figsize=(10, 5)
    else:
        # vertically arrange plots
        orientation='horizontal'
        location='top'
        ax1n = 211
        ax2n = 212
        figsize=(7, 7)
    
    # ---------- plot at each time ----------

    data_crs = ccrs.PlateCarree()
    Tlev050 = np.arange(185, 216, 3)
    Tlev850 = np.arange(260, 306, 5)
    #cmap = claut.ncar_rgb_to_cmap(gmt)
    cmap = mpl.cm.rainbow
    #cmap = mpl.cm.Spectral_r
    pltargs050 = {'levels':Tlev050, 'cmap':cmap, 'zorder':0}
    pltargs850 = {'levels':Tlev850, 'cmap':cmap, 'zorder':0}
    cArgs = {'orientation':orientation, 'location':location, 'label':'T [K]', 'aspect':30}
        
    figT = plt.figure(figsize=figsize)
        
    print('\n\n=============== {}'.format(title)) 

    for t in range(len(th)):
        print('------- t: {}'.format(th.values[t]))
        if(th.values[t] % 6 != 0):
            continue
        
        plt.clf()
        ax050 = figT.add_subplot(ax1n, projection=data_crs)
        ax850 = figT.add_subplot(ax2n, projection=data_crs)
        
        
        var_dict = [{'var':T050[t], 'plotType':'contourf', 'plotArgs':pltargs050, 'colorArgs':cArgs}]
        cf = plthor(lon, lat, var_dict, ax=ax050, annotation='50 hPa', xlabel='', coastlines=False)
        cf[0].set_ticks(Tlev050)
        cf[0].ax.set_xticklabels(Tlev050, rotation=90)
        #ax050.set_ylabel('lat [deg]')
        #ax050.set_xlabel('lon [deg]')
        
        #if(vertarr):
        #    var_dict = [{'var':T850[t], 'plotType':'contourf', 'plotArgs':pltargs850, 'colorFormatter':None}]
        #else:
        var_dict = [{'var':T850[t], 'plotType':'contourf', 'plotArgs':pltargs850, 'colorArgs':cArgs}]
        cf = plthor(lon, lat, var_dict, ax=ax850, annotation='850 hPa', xlabel='', coastlines=False)
        #if(not vertarr):
        cf[0].set_ticks(Tlev850)
        cf[0].ax.set_xticklabels(Tlev850,rotation=90)
        #ax850.set_ylabel('lat [deg]')
        #ax850.set_xlabel('lon [deg]')
       
        #plt.subplots_adjust(wspace=0.2)
        #figT.suptitle('{}, hour {:.0f}'.format(title, th.values[t]), fontsize=14)
        plt.tight_layout()
        plt.savefig('{}/{}_hour{:3.0f}.png'.format(savedest, title, th.values[t]), dpi=150)
        #plt.show()
        





if(__name__ == '__main__'):

    #gmt = '../../SAI/analysis/cmaps/GMT_no_green.rgb'
    eam_hsw = '/global/cscratch1/sd/jhollo/E3SM/E3SMv2_cases/hsw_cases/pre_rebase/'\
          'E3SM_ne16_L72_FIDEAL_30day/run' 
    hsw_6hrmeans = '{}/E3SM_ne16_L72_FIDEAL_30day.eam.h0.0001-01-01-00000.regrid.2x2.nc'.format(eam_hsw)
    #eam_hsw = '/global/cscratch1/sd/jhollo/E3SM/E3SMv2_cases/sai_cases/'\
    #         'HSW_SAI_ne16pg2_L72_passive/run' 
    #hsw_6hrmeans = '{}/HSW_SAI_ne16pg2_L72_passive.eam.h1.0001-01-01-00000.'\
    #               'regrid.91x180_aave.nc'.format(eam_hsw)
    dest = './figs/pre_rebase/eam_hsw_30day_animation'
    circulation_snapshot(hsw_6hrmeans, 'EAM HSW ne16L72, 30-day evolution', dest, vertarr=True)
