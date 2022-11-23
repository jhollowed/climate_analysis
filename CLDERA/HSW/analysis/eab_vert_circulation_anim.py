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
from matplotlib import colors

# ============================================================


def circulation_snapshot(runf, title, savedest):
   

    # read data, transform time coordinate to ndays
    try:
        dat = xr.open_dataset('./tmpdata/U_zonalMean_30day.nc')
        th = dat['time']
        U = dat['U']
        print('read tmp U data')
    except FileNotFoundError:
        print('reading data')
        dat = xr.open_dataset(runf)
        th = dat['nsteph']*1800/60/60
        dat = dat.assign_coords(time=th)
        
        print('taking mean of U data')
        U   = dat['U'].sel({'time':slice(0, 30*24)}).mean('lon')
        print('writing U data')
        U.to_netcdf('./tmpdata/U_zonalMean_30day.nc')
    
    lat = dat['lat']
    lev = dat['lev']

    ulev = np.hstack([np.linspace(-35, 0, 8), np.linspace(15, 90, 6)]).astype(int)
    divnorm=colors.TwoSlopeNorm(vcenter=0.)


    # ---------- plot at each time ----------

    print('Plotting')
    fig = plt.figure(figsize=(6, 5))
    cmap = mpl.cm.rainbow
    pltargs = {'levels':ulev, 'cmap':cmap, 'zorder':0, 'norm':divnorm}
    cArgs = {'orientation':'horizontal', 'location':'top', 'label':'U [m/s]','aspect':30,'format':'%d'}
        
    print('\n\n=============== {}'.format(title)) 

    for t in range(len(th)):
        print('------- t: {}'.format(th.values[t]))
        if(th.values[t] % 6 != 0):
            continue
        plt.clf()
        ax = fig.add_subplot(111)
        ttitle = 'day {}'.format(int(th.values[t]/24))

        var_dict = [{'var':U[t], 'plotType':'contourf', 'plotArgs':pltargs, 'colorArgs':cArgs}]
        cf = pltvert(lat, lev, var_dict, ax=ax, plot_zscale=True, annotation=None)
        cf[0].set_ticks(ulev)
        ax.set_ylabel('p  [hPa]')
        ax.set_xlabel('latitude')
        
        fig.suptitle('{}\n{}'.format(title, ttitle), fontsize=14)
        fig.tight_layout()
        fig.savefig('{}/{}_t{:03d}_U.png'.format(savedest, title, int(th.values[t])), dpi=150)
        





if(__name__ == '__main__'):

    eam_hsw = '/global/cscratch1/sd/jhollo/E3SM/E3SMv2
    cases/sai_cases/'\
              'HSW_SAI_ne16pg2_L72_passive/run' 
    hsw_6hrmeans = '{}/HSW_SAI_ne16pg2_L72_passive.eam.h1.0001-01-01-00000.'\
                   'regrid.91x180_aave.nc'.format(eam_hsw)
    dest = './figs/eab_vert_anims'
    circulation_snapshot(hsw_6hrmeans, 'EAM HSW ne16L72, 30-day U evolution', dest)
