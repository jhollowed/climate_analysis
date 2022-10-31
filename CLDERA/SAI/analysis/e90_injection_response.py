import pdb
import glob
import numpy as np
import xarray as xr
import matplotlib as mpl
import artist_utils as aut
import cartopy.crs as ccrs
import climate_toolbox as ctb
from matplotlib import colors
import matplotlib.pyplot as plt
from metpy.units import units as u
import matplotlib.ticker as ticker
from matplotlib.offsetbox import AnchoredText
from climate_artist import vertical_slice as pltvert
from climate_artist import horizontal_slice as plthor

# ============================================================


def e90_injection_response(run1, run2, sfx, overwrite=False, 
                           figsavedest=None, datsavedest='.', inj_delay=0):

    # windows in days over which to compoute means for horizontal slices
    # second month, fourth month
    # offset if injection_dealy > 0
    t_windows = [slice(30+inj_delay, 60+inj_delay), slice(120+inj_delay, 150+inj_delay)]
    t_labels = ['month 2\n', 'month 5\n']
    t_labels = ['', '']

    # params
    kso2 = 1/((25 * u.day).to(u.s))
    ksulf = 1/((360*u.day).to(u.s))
    tf = (9 * u.hr).to(u.s)
    Mso2 = 1.7e10 * u.kg
    w = 2.04

    # read data, transform time coordinate to ndays
    print('reading data 1')
    dat1 = xr.open_dataset(run1)
    td1 = ctb.time2day(dat1['time'])
    dat1 = dat1.assign_coords(time=td1)
    print('found {} timesteps from day {} to {}'.format(len(td1), min(td1), max(td1)))
    
    print('reading data 2')
    dat2 = xr.open_dataset(run2)
    td2 = ctb.time2day(dat2['time'])
    dat2 = dat2.assign_coords(time=td2)
    print('found {} timesteps from day {} to {}'.format(len(td2), min(td2), max(td2)))

    # assume that these are the same for both datasets
    time = td1
    lat = dat1['lat']
    lev = dat1['lev']
    
    # get "anomalies", read from file if previosly computed, or compute and then write out
    try:
        if(overwrite == True): raise FileNotFoundError
        E90diff = xr.open_dataset('{}/E90diff_{}.nc'.format(datsavedest, sfx))['E90diff']
        print('read E90diff')
    except FileNotFoundError:
        print('computing E90 diff...')
        # ---- diff of zonal means
        E901   = dat1['E90'].mean('lon') * 1e9
        E902   = dat2['E90'].mean('lon') * 1e9
        E90diff = E901 - E902
        E90diff.to_dataset(name='E90diff').to_netcdf('{}/E90diff_{}.nc'.format(datsavedest, sfx))
        print('wrote udiff')


    # ==========================  plotting  ==================================

    print('plotting...')
    # define colormaps, colormap levels, format, and contour labels per-dataset
    
    #clev_E90_divnorm=colors.TwoSlopeNorm(vmin=-10, vcenter=0., vmax=30)
    clev_E90 = np.linspace(-10, 10, 11)
    clev_E90_fmt = '%.0f'
    clev_E90_labels = []
    cmap_E90 = mpl.cm.rainbow
    
    label_fs = 14
    tick_fs = 12
 
    # first 2 rows:
    # heating, T anom, U anom, for 2 30-day windows
    # last 2 rows:
    # T anom vs t, T anom vs t at select levels
    fig   = plt.figure(figsize=(10, 5))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    
    print('PLOTTING E90')
    
    pltargs = {'levels':clev_E90, 'cmap':cmap_E90, 'zorder':0}#, 'norm':clev_E90_divnorm}
    pltargs_c = {'levels':clev_E90, 'colors':'k', 'linewidths':0.6, 'zorder':1}
    cArgs = {'orientation':'horizontal', 'location':'top', 
             'label':'{}E90 anomaly [ppb]'.format(t_labels[0]), 'aspect':20, 'format':clev_E90_fmt}
    cArgs_c = {'fmt':clev_E90_fmt, 'levels':clev_E90_labels}
    var_dict = [{'var':E90diff.sel({'time':t_windows[0]}).mean('time'), 
                 'plotType':'contourf', 'plotArgs':pltargs, 'colorArgs':cArgs}]
    var_dict_c = [{'var':E90diff.sel({'time':t_windows[0]}).mean('time'), 
                   'plotType':'contour', 'plotArgs':pltargs_c, 'colorArgs':cArgs_c}]    
    cf = pltvert(lat, lev, var_dict, ax=ax1, plot_zscale=False,annotation='', gridlines=True)
    pltvert(lat, lev, var_dict_c, ax=ax1, plot_zscale=False, inverty=False, annotation='')
    cf[0].set_ticks(clev_E90)
    cf[0].ax.tick_params(rotation=90)
    ax1.set_ylim([10, 900])
    ax1.invert_yaxis()
    
    # lower plot (second time window) will have no colobar; levels agree with above
    var_dict = [{'var':E90diff.sel({'time':t_windows[1]}).mean('time'), 
        'plotType':'contourf', 'plotArgs':pltargs, 'colorArgs':cArgs}]#, 'colorFormatter':None}]
    var_dict_c = [{'var':E90diff.sel({'time':t_windows[1]}).mean('time'), 
                   'plotType':'contour', 'plotArgs':pltargs_c, 'colorArgs':cArgs_c}]    
    cArgs['label'] = '{}E90 anomaly [ppb]'.format(t_labels[1])
    cf = pltvert(lat, lev, var_dict, ax=ax2, plot_zscale=False,annotation='', gridlines=True)
    pltvert(lat, lev, var_dict_c, ax=ax2, plot_zscale=False, inverty=False, annotation='')
    cf[0].set_ticks(clev_E90)
    cf[0].ax.tick_params(rotation=90) 
    ax2.set_ylim([10, 900])
    ax2.invert_yaxis()
    
    aut.format_ticks(ax1)
    aut.format_ticks(ax2) 
    ax1.set_ylabel('pressure  [hPa]')
    ax2.set_ylabel('pressure  [hPa]')
    ax2.yaxis.set_label_position("right")
    ax2.yaxis.tick_right() 
    ax1.set_xlabel('latitude')
    ax2.set_xlabel('latitude')
    fig.suptitle(sfx)
     
    plt.tight_layout()
    if(figsavedest is not None):
        fig.savefig('{}/e90_injection_response{}.png'.format(figsavedest, sfx), dpi=300)
    else: 
        plt.show()
    

    print('PLOTTING E90 LINES')
    lineax.plot(time, E90diff.sel(lat=slice(-15,15)).mean('lat').sel(lev=80, method='nearest'), 
                label=linelabel)


if(__name__ == '__main__'):
    
    linefig = plt.figure()
    lineax =linefig.add_subplot(111)
    configs = ['17', '34', '85', '170']
    linelabels = ['1x SO2', '2x SO2', '5x SO2', '10x SO2']

    fig_dest = '/global/homes/j/jhollo/repos/climate_analysis/CLDERA/SAI/analysis/figs/'\
               'heating_response'
    data_dest = '/global/cscratch1/sd/jhollo/E3SM/E3SMv2_cases/sai_cases/processes_pathways/e90'
    data_source = '/global/cscratch1/sd/jhollo/E3SM/E3SMv2_cases/sai_cases/HSW_SAI_ne16pg2_L72_heatingEns_e90st80_CLONES'

    for j in range(len(configs)):
        linelabel = linelabels[j]
        runs1_config = 'heatingEns_e90st80__cldera_sai_MSO2_{}'.format(configs[j])
        runs2_config = 'heatingEns_e90st80_passive'
        sfx = '_{}-{}'.format(runs1_config, runs2_config)
        
        rundir1 = '{}/HSW_SAI_ne16pg2_L72_{}/run/'.format(data_source, runs1_config)
        rundir2 = '{}/HSW_SAI_ne16pg2_L72_{}/run/'.format(data_source, runs2_config)
        runs1 = glob.glob('{}/*regrid*aave*nc'.format(rundir1))
        runs2 = glob.glob('{}/*regrid*aave*.nc'.format(rundir2))
        
        if(len(runs1) > 1):
            run1 = '{}/{}_concat_hist.nc'.format(data_dest, runs1_config)
            ctb.concat_run_outputs(runs1, outFile=run1, histnum=0, regridded=True, component='eam')
        else: 
            run1 = runs1[0]
        if(len(runs2) > 1):
            run2 = '{}/{}_concat_hist.nc'.format(data_dest, runs2_config)
            ctb.concat_run_outputs(runs2, outFile=run2, histnum=0, regridded=True, component='eam')
        else: 
            run2 = runs2[0]
        print('\n')
        
        e90_injection_response(run1, run2, overwrite=False, sfx=sfx, 
                               figsavedest=fig_dest, datsavedest=data_dest, inj_delay=0) 

    lineax.grid(alpha=0.5)
    lineax.set_xlabel('time [days]')
    lineax.set_ylabel('E90 anomaly [ppb]')
    lineax.set_title('zonal mean, +-15 deg meridional mean, at 80 hPa')
    lineax.legend()
    plt.show()
