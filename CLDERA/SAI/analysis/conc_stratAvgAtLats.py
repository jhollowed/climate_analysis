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
from wrappers import dpres_hybrid_ccm as pwgt

# ============================================================


def conc_strat(runf, title, savedest, tracer='SO2', globe=True):

    run = xr.open_dataset(runf)
    td = ctb.time2day(run['time'])

    c = (run['SAI_{}'.format(tracer)]).mean('lon')
    
    lat = run['lat']
    lev = run['lev']
    lat0 = 15.15
    dlat = 0.5
    cmap = ctb.ncar_rgb_to_cmap(gmt)
    
    # average over the steratosphere w/ weighting
    #hyaivert = hyai.sel({'ilev':slice(1, 70)})
    #hybivert = hybi.sel({'ilev':slice(1, 70)})
    #pdb.set_trace()
    #pp = pwgt(PS, P0, hyaivert, hybivert)
    
    latsel = [-60, 0, 15, 60]
    cvert = c.sel({'lev':45}, method='nearest')
    cvert = np.log10(cvert.sel({'lat':latsel}, method='nearest')).T
    #cvert = np.ma.masked_array(cvert, np.logical_or(cvert == -float('inf'), 
    #                                                cvert < -20))
    latt = cvert['lat']

    colors = cmap(np.linspace(0.33, 0.95, len(latt)))
    for i in range(len(cvert)):
        plt.plot(td, cvert[i], label = 'lat = {:.0f} deg'.format(latt.values[i]), color=colors[i])
    plt.xlabel('t  [days]', fontsize=14)
    plt.ylabel('log10(c_SO2  [kg/kg])', fontsize=14)
    plt.legend()
    plt.title('{}, 45 hPa'.format(title))
    plt.ylim([-10, -4])
    plt.savefig('{}/{}.png'.format(savedest, title.split()[-1]), dpi=300)
    return
    
    clog = np.log10(c)
    clog = np.ma.masked_array(clog, np.logical_or(clog == -float('inf'), 
                                                    clog < -20))
     

    # get time in number of days
    
    # ---------- plot ----------
    
    lon0 = 120.35
    minc = -12
    maxc = -4
    clevels = 8
    levels = np.linspace(minc, maxc, clevels)
    data_crs = ccrs.PlateCarree()
        
    print('\n\n=============== {}'.format(title)) 

    for k in range(len(td)):
        
        if(k%10) != 0: continue

        fig = plt.figure(figsize=(7,6))
        ax1 = fig.add_subplot(111)

        # ----- vertical
        pltargs = {'levels':levels, 'cmap':cmap}
        pltargs_c = {'levels':levels, 'colors':'k', 'linewidths':0.6, 'linestyles':'-'}
        cArgs = {'orientation':'horizontal', 'location':'bottom', 'label':'log10(c_SO2 [kg/kg])'}
        var_dict = [{'var':clog[k], 'plotType':'contourf', 'plotArgs':pltargs, 'colorArgs':cArgs}]
        var_dict_c = [{'var':clog[k], 'plotType':'contour', 'plotArgs':pltargs_c, 'colorFormatter':None}]
        pltvert(lat, lev, var_dict, ax=ax1)
        pltvert(lat, lev, var_dict_c, ax=ax1, plot_zscale=False)
        text_box = AnchoredText('lon=120 deg', frameon=True, loc=4, pad=0.5)
        text_box.set_zorder(9)
        plt.setp(text_box.patch, facecolor='white', alpha=1)
        ax1.add_artist(text_box)
        ax1.set_ylabel('p  [hPa]')
        ax1.set_xlabel('lat  [deg]')

        fig.suptitle('{}, day {}'.format(title, round(td[k]+0.)), fontsize=16)
        plt.tight_layout()
        plt.savefig('{}/{}{:03d}.png'.format(savedest, title.split()[-1], k), dpi=300)
        





if(__name__ == '__main__'):

    gmt = '/glade/u/home/jhollowed/repos/climate_analysis/CLDERA/SAI/analysis/GMT_no_green.rgb'
    whs = '/glade/scratch/jhollowed/CAM/cases/sai_runs/SE_ne16L72_whs_saiv2_fix0_tau0_qsplit1/'\
          'run/SE_ne16L72_whs_saiv2_fix0_tau0_qsplit1.cam.h0.0001-01-01-00000.nc'
    whsdest = '/glade/u/home/jhollowed/repos/climate_analysis/CLDERA/SAI/analysis/figs/'
    amip = '/glade/scratch/jhollowed/CAM/cases/sai_runs/E3SM_AMIP_ne30_L72_SAI/'\
           'AMIPcase_ne30_L72_SAI.eam.h0.0001-01-01-00000.nc.regrid.2x2.nc'
    amipdest = '/glade/u/home/jhollowed/repos/climate_analysis/CLDERA/SAI/analysis/figs/'
    
    conc_strat(whs, 'SE ne30L72, CAM WHS', whsdest)
    conc_strat(amip, 'SE ne30L72, EAM AMIP', amipdest)
    #animate_eruption(whs, 'SE ne30L72, CAM WHS', whsdest, globe=False)
    #animate_eruption(amip, 'SE ne30L72, EAM AMIP', amipdest, globe=False)
    
