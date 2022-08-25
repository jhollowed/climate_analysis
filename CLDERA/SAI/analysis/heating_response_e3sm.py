import pdb
import glob
import numpy as np
import xarray as xr
import matplotlib as mpl
import artist_utils as aut
import cartopy.crs as ccrs
import artist_utils as claut
import climate_toolbox as ctb
from matplotlib import colors
import matplotlib.pyplot as plt
from metpy.units import units as u
import matplotlib.ticker as ticker
from matplotlib.offsetbox import AnchoredText
from climate_artist import vertical_slice as pltvert
from climate_artist import horizontal_slice as plthor

# ============================================================


def prelim_e3sm_fig(run1, run2, sfx='', savedest=None, datsavedest='.', inj_delay=0):

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

    # toggle to force re-computation
    overwrite=False
    file_ext = 'pathways_{}.nc'.format(sfx)
    
    # get "anomalies", read from file if previosly computed, or compute and then write out
    try:
        if(overwrite == True): raise FileNotFoundError
        Udiff = xr.open_dataset('{}/Udiff_{}'.format(datsavedest, file_ext))['Udiff']
        print('read Udiff')
    except FileNotFoundError:
        print('computing udiff...')
        # ---- diff of zonal means
        U1   = dat1['U'].mean('lon')
        U2   = dat2['U'].mean('lon')
        Udiff = U1 - U2
        Udiff.to_dataset(name='Udiff').to_netcdf('{}/Udiff_{}'.format(datsavedest, file_ext))
        print('wrote udiff')

    try:
        if(overwrite == True): raise FileNotFoundError
        Tdiff = xr.open_dataset('{}/Tdiff_{}'.format(datsavedest, file_ext))['Tdiff']
        print('read Tdiff')
        Tdiff_global = xr.open_dataset('{}/Tdiff_global_{}'.format(datsavedest, file_ext))['Tdiff_global']
        print('read Tdiff_global')
    except FileNotFoundError:
        print('computing Tdiff...')
        # ---- diff of zonal means
        T1   = dat1['T'].mean('lon')
        T2   = dat2['T'].mean('lon')
        Tdiff = T1-T2
        Tdiff.to_dataset(name='Tdiff').to_netcdf('{}/Tdiff_{}'.format(datsavedest, file_ext)) 
        print('wrote Tdiff')
        # ---- get global T anom via weighted avg
        print('computing Tdiff_global...')
        weights = np.cos(np.deg2rad(lat))
        weights.name = 'weights'
        Tdiff_global  = (Tdiff.weighted(weights)).mean('lat')
        Tdiff_global.to_dataset(name='Tdiff_global').to_netcdf(
                                '{}/Tdiff_global_{}'.format(datsavedest, file_ext))
        print('wrote Tdiff_global')
    
    try:
        if(overwrite == True): raise FileNotFoundError
        Hdiff = xr.open_dataset('{}/Hdiff_{}'.format(datsavedest, file_ext))['Hdiff']
        print('read Hdiff')
    except FileNotFoundError:
        print('computing Hdiff...', end='\r')
        H1   = dat1['SAI_HEAT']+dat1['SAI_COOL']
        H1   = H1.mean('lon')
        print('computing Hdiff... H1', end='\r')
        H2   = dat2['SAI_HEAT']+dat2['SAI_COOL']
        H2   = H2.mean('lon')
        print('computing Hdiff... H1 H2')
        # ---- diff of zonal means
        Hdiff = H1-H2
        Hdiff.to_dataset(name='Hdiff').to_netcdf('{}/Hdiff_{}'.format(datsavedest, file_ext))
        print('wrote Hdiff')


    # get total masses of SO2, sulfate for all time, compute analytic solution
    try:
        if(overwrite == True): raise FileNotFoundError
        SO21 = xr.open_dataset('{}/so21_{}'.format(datsavedest, file_ext))['tot_SO2_mass']
        sulf1 = xr.open_dataset('{}/sulf1_{}'.format(datsavedest, file_ext))['tot_SULFATE_mass']
        print('read tracers 1')
    except FileNotFoundError:
        print('computing tracers 1...')
        SO21 = (dat1['SO2'] * dat1['AIR_MASS']).sum(('lat', 'lon', 'lev'))
        sulf1 = (dat1['SULFATE'] * dat1['AIR_MASS']).sum(('lat', 'lon', 'lev'))
        SO21.to_dataset(name='tot_SO2_mass').to_netcdf('{}/so21_{}'.format(datsavedest, file_ext))
        sulf1.to_dataset(name='tot_SULFATE_mass').to_netcdf('{}/sulf1_{}'.format(datsavedest, file_ext))
        print('wrote tracers 1')
    
    try:
        if(overwrite == True): raise FileNotFoundError
        SO22 = xr.open_dataset('{}/so22_{}'.format(datsavedest, file_ext))['tot_SO2_mass']
        sulf2 = xr.open_dataset('{}/sulf2_{}'.format(datsavedest, file_ext))['tot_SULFATE_mass']
        print('read tracers 2')
    except FileNotFoundError:
        print('computing tracers 2...')
        SO22 = (dat2['SO2'] * dat2['AIR_MASS']).sum(('lat', 'lon', 'lev'))
        sulf2 = (dat2['SULFATE'] * dat2['AIR_MASS']).sum(('lat', 'lon', 'lev'))
        SO22.to_dataset(name='tot_SO2_mass').to_netcdf('{}/so22_{}'.format(datsavedest, file_ext))
        sulf2.to_dataset(name='tot_SULFATE_mass').to_netcdf('{}/sulf2_{}'.format(datsavedest, file_ext))
        print('wrote tracers 2')
    
    # analytic total mass with V(:) = 1, sum Vk = 1
    print('computing analytic tracers...')
    ttime = (time * u.day).to(u.s)
    tmin = np.array([min(t, tf).m for t in ttime]) * ttime[0].u
    V = 1
    A = Mso2 / tf
    SO2_analytic = A/kso2 * V * np.exp(-kso2*ttime) * (np.exp(kso2 * tmin) - 1) 
    sulf_analytic = (w*A/((ksulf - kso2)*ksulf) * V * \
                        np.exp(-ksulf*ttime) * (\
                        kso2*(1 - np.exp(ksulf * tmin)) - ksulf*np.exp((ksulf-kso2)*ttime) * \
                        (1 - np.exp(kso2 * tmin))))

    # ==========================  plotting  ==================================

    print('plotting...')
    # define colormaps, colormap levels, format, and contour labels per-dataset
    clev_H = np.hstack([np.linspace(-0.03, 0, 4), 
                        np.linspace(0.05, 0.3, 6)])
    clev_H_divnorm=colors.TwoSlopeNorm(vmin=-0.03, vcenter=0., vmax=0.3)
    clev_H_fmt = '%.2f'
    clev_H_labels = []
    cmap_H = mpl.cm.RdYlGn_r
    
    clev_T = np.linspace(-8, 8, 9)
    clev_T_fmt = '%.0f'
    clev_T_labels = []
    cmap_T = mpl.cm.RdYlBu_r
    
    clev_U_divnorm=colors.TwoSlopeNorm(vmin=-10, vcenter=0., vmax=30)
    clev_U = np.linspace(-10, 30, 9)
    clev_U_fmt = '%.0f'
    clev_U_labels = []
    cmap_U = mpl.cm.rainbow
    
    clev_Tg = np.linspace(-3, 4, 8)
    clev_Tg_fmt = '%.0f'
    clev_Tg_labels = []
    
    label_fs = 14
    tick_fs = 12

    
    # first 2 rows:
    # heating, T anom, U anom, for 2 30-day windows
    # last 2 rows:
    # T anom vs t, T anom vs t at select levels
    fig   = plt.figure(figsize=(8, 14))
    gs    = fig.add_gridspec(nrows=5, ncols=3)
    axH0  = fig.add_subplot(gs[0, 0])
    axT0  = fig.add_subplot(gs[0, 1])
    axU0  = fig.add_subplot(gs[0, 2])
    axH1  = fig.add_subplot(gs[1, 0])
    axT1  = fig.add_subplot(gs[1, 1])
    axU1  = fig.add_subplot(gs[1, 2])
    axTg  = fig.add_subplot(gs[2, :])
    axTgl = fig.add_subplot(gs[3, :])
    axTr  = fig.add_subplot(gs[4, :])
    top_axes = [axH0, axT0, axU0, axH1, axT1, axU1]
    bottom_axes = [axTg, axTgl, axTr]
    
    # ---------- HDIFF PLOTS ----------
    print('PLOTTING HEATING')

    pltargs = {'levels':clev_H, 'cmap':cmap_H, 'zorder':0, 'norm':clev_H_divnorm}
    pltargs_c = {'levels':clev_H, 'colors':'k', 'linewidths':0.6, 'zorder':1}
    cArgs = {'orientation':'horizontal', 'location':'top', 
             'label':'{}Heating  rate [K/day]'.format(t_labels[0]), 'aspect':20, 'format':clev_H_fmt}
    cArgs_c = {'fmt':clev_H_fmt, 'levels':clev_H_labels}
    var_dict = [{'var':Hdiff.sel({'time':t_windows[0]}).mean('time'), 
                 'plotType':'contourf', 'plotArgs':pltargs, 'colorArgs':cArgs}]
    var_dict_c = [{'var':Hdiff.sel({'time':t_windows[0]}).mean('time'), 
                   'plotType':'contour', 'plotArgs':pltargs_c, 'colorArgs':cArgs_c}]    
    cf = pltvert(lat, lev, var_dict, ax=axH0, plot_zscale=False, annotation='', gridlines=True)
    pltvert(lat, lev, var_dict_c, ax=axH0, plot_zscale=False, inverty=False, annotation='')
    cf[0].set_ticks(clev_H)
    cf[0].ax.tick_params(rotation=90)
    
    # lower plot (second time window) will have no colobar; levels agree with above
    var_dict = [{'var':Hdiff.sel({'time':t_windows[1]}).mean('time'), 
        'plotType':'contourf', 'plotArgs':pltargs, 'colorArgs':cArgs}]#, 'colorFormatter':None}]
    var_dict_c = [{'var':Hdiff.sel({'time':t_windows[1]}).mean('time'), 
                   'plotType':'contour', 'plotArgs':pltargs_c, 'colorArgs':cArgs_c}]
    cArgs['label'] = '{}Heating  rate [K/day]'.format(t_labels[1])
    cf = pltvert(lat, lev, var_dict, ax=axH1, plot_zscale=False,annotation='', gridlines=True)
    pltvert(lat, lev, var_dict_c, ax=axH1, plot_zscale=False, inverty=False, annotation='')
    cf[0].set_ticks(clev_H)
    cf[0].ax.tick_params(rotation=90)
    
    # ---------- TDIFF PLOTS ----------
    print('PLOTTING TEMP')

    pltargs = {'levels':clev_T, 'cmap':cmap_T, 'zorder':0}
    pltargs_c = {'levels':clev_T, 'colors':'k', 'linewidths':0.6, 'zorder':1}
    cArgs = {'orientation':'horizontal', 'location':'top', 
             'label':'{}T anomaly [K]'.format(t_labels[0]), 'aspect':20, 'format':clev_T_fmt}
    cArgs_c = {'fmt':clev_T_fmt, 'levels':clev_T_labels}
    var_dict = [{'var':Tdiff.sel({'time':t_windows[0]}).mean('time'), 
                 'plotType':'contourf', 'plotArgs':pltargs, 'colorArgs':cArgs}]
    var_dict_c = [{'var':Tdiff.sel({'time':t_windows[0]}).mean('time'), 
                   'plotType':'contour', 'plotArgs':pltargs_c, 'colorArgs':cArgs_c}]    
    cf = pltvert(lat, lev, var_dict, ax=axT0, plot_zscale=False,annotation='', gridlines=True)
    pltvert(lat, lev, var_dict_c, ax=axT0, plot_zscale=False, inverty=False, annotation='')
    cf[0].set_ticks(clev_T)
    cf[0].ax.tick_params(rotation=90)
    
    # lower plot (second time window) will have no colobar; levels agree with above
    var_dict = [{'var':Tdiff.sel({'time':t_windows[1]}).mean('time'), 
        'plotType':'contourf', 'plotArgs':pltargs, 'colorArgs':cArgs}]#, 'colorFormatter':None}]
    var_dict_c = [{'var':Tdiff.sel({'time':t_windows[1]}).mean('time'), 
                   'plotType':'contour', 'plotArgs':pltargs_c, 'colorArgs':cArgs_c}]    
    cArgs['label'] = '{}Heating  rate [K/day]'.format(t_labels[1])
    cf = pltvert(lat, lev, var_dict, ax=axT1, plot_zscale=False,annotation='', gridlines=True)
    pltvert(lat, lev, var_dict_c, ax=axT1, plot_zscale=False, inverty=False, annotation='')
    cf[0].set_ticks(clev_T)
    cf[0].ax.tick_params(rotation=90)
    
    # ---------- UDIFF PLOTS ----------
    print('PLOTTING WIND')
    
    pltargs = {'levels':clev_U, 'cmap':cmap_U, 'zorder':0, 'norm':clev_U_divnorm}
    pltargs_c = {'levels':clev_U, 'colors':'k', 'linewidths':0.6, 'zorder':1}
    cArgs = {'orientation':'horizontal', 'location':'top', 
             'label':'{}U anomaly [m/s]'.format(t_labels[0]), 'aspect':20, 'format':clev_U_fmt}
    cArgs_c = {'fmt':clev_U_fmt, 'levels':clev_U_labels}
    var_dict = [{'var':Udiff.sel({'time':t_windows[0]}).mean('time'), 
                 'plotType':'contourf', 'plotArgs':pltargs, 'colorArgs':cArgs}]
    var_dict_c = [{'var':Udiff.sel({'time':t_windows[0]}).mean('time'), 
                   'plotType':'contour', 'plotArgs':pltargs_c, 'colorArgs':cArgs_c}]    
    cf = pltvert(lat, lev, var_dict, ax=axU0, plot_zscale=False,annotation='', gridlines=True)
    pltvert(lat, lev, var_dict_c, ax=axU0, plot_zscale=False, inverty=False, annotation='')
    cf[0].set_ticks(clev_U)
    cf[0].ax.tick_params(rotation=90)
    
    # lower plot (second time window) will have no colobar; levels agree with above
    var_dict = [{'var':Udiff.sel({'time':t_windows[1]}).mean('time'), 
        'plotType':'contourf', 'plotArgs':pltargs, 'colorArgs':cArgs}]#, 'colorFormatter':None}]
    var_dict_c = [{'var':Udiff.sel({'time':t_windows[1]}).mean('time'), 
                   'plotType':'contour', 'plotArgs':pltargs_c, 'colorArgs':cArgs_c}]    
    cArgs['label'] = '{}Heating  rate [K/day]'.format(t_labels[1])
    cf = pltvert(lat, lev, var_dict, ax=axU1, plot_zscale=False,annotation='', gridlines=True)
    pltvert(lat, lev, var_dict_c, ax=axU1, plot_zscale=False, inverty=False, annotation='')
    cf[0].set_ticks(clev_U)
    cf[0].ax.tick_params(rotation=90)
    
    
    # ---------- global T diff plot ----------
    print('PLOTTING GLOBAL TEMP')
    
    TIME, LEV = np.meshgrid(time, lev)
    cf = axTg.contourf(TIME, LEV, Tdiff_global.T, cmap=cmap_T, levels=clev_Tg, extend='both')
    cfl = axTg.contour(TIME, LEV, Tdiff_global.T, levels=clev_Tg, colors='k', linewidths=0.6, zorder=1)
    cb = fig.colorbar(cf, ax=axTg, orientation='horizontal', location='top',
                      aspect=40, extend='both', extendrect=False, format=clev_Tg_fmt)
    cb.set_label('Global mean T anomaly [K]', fontsize=tick_fs)
    axTg.clabel(cfl, levels=clev_Tg_labels, inline=True, fmt=clev_Tg_fmt, fontsize=tick_fs)
    cb.set_ticks(clev_Tg)

    axTg.set_yscale('log')
    axTg.invert_yaxis()
    axTg.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y,pos: \
                                   ('{{:.{:1d}f}}'.format(int(np.maximum(-np.log10(y),0)))).format(y)))
    axTg.set_ylabel('pressure [hPa]', fontsize=label_fs)
    
    # ---------- global T diff plot at levels ----------
    print('PLOTTING GLOBAL NLEV TEMP')
   
    cls = plt.cm.viridis(np.linspace(0.2, 0.8, 4))
    Tg_lev = [3, 10, 30, 500]
    for i in range(len(Tg_lev)):
        axTgl.plot(time, Tdiff_global.sel({'lev':Tg_lev[i]}, method='nearest'), 
                   label='{} hPa'.format(Tg_lev[i]), color=cls[i])
    axTgl.set_ylabel('Global Mean\nT anomaly [K]', fontsize=label_fs)
    axTgl.plot(axTgl.get_xlim(), [0, 0], '--k', lw=2)
    axTgl.legend(fontsize=tick_fs, fancybox=False, loc='upper left', ncol=2)
    
    # ---------- total tracer masses ----------
    print('PLOTTING TRACER MASSES')
   
    axTr.plot(time, SO21/1e9, '-r',  label='SO$_2$, active')
    axTr.plot(time, SO22/1e9, '--r', label='SO$_2$, passive')
    axTr.plot(time, sulf2/1e9,'-', color='orange',label='sulfate, active')
    axTr.plot(time, sulf2/1e9,'--', color='orange', label='sulfate, passive')
    axTr.set_ylabel('Total mass [Tg]', fontsize=label_fs)
    axTr.set_xlabel('time [days]', fontsize=label_fs)
    axTr.legend(fontsize=tick_fs, fancybox=False, loc='center right', ncol=2)
    #axTr.set_ylim([0, 30])
    
    
    # -------- formatting --------
    for ax in top_axes:
        ax.set_xlabel('')
        ax.set_ylabel('')
        aut.format_ticks(ax)
    
    for ax in [axH0, axH1, axU0, axU1]:
        ax.set_ylabel('pressure  [hPa]')
    for ax in [axT0, axT1]:
        ax.yaxis.set_ticklabels([])
    for ax in [axU0, axU1]:
         ax.yaxis.set_label_position("right")
         ax.yaxis.tick_right()
    
    for ax in [axH1, axT1, axU1]:
        ax.set_xlabel('latitude')
    for ax in [axH0, axT0, axU0]:
        ax.xaxis.set_ticklabels([])
    
    for ax in [axTg, axTgl, axTr]:
        ax.set_xlim([0, 180])
        aut.format_ticks(ax)
        ax.grid(lw=0.3, alpha=0.75)
     
    plt.tight_layout()
    if(savedest is not None):
        fig.savefig('{}/prelim_e3smFig_pthwy_{}.png'.format(savedest, sfx), dpi=300)
    else: 
        plt.show()
        





if(__name__ == '__main__'):

    dest = '/global/homes/j/jhollo/repos/climate_analysis/CLDERA/SAI/analysis/figs/'\
           'e3sm_pathway_figs'
    data = '/global/cscratch1/sd/jhollo/E3SM/E3SMv2_cases/sai_cases'
    ds = '/global/cscratch1/sd/jhollo/E3SM/E3SMv2_cases/sai_cases/processes_pathways'
    #ds = '/global/cscratch1/sd/jhollo/E3SM/E3SMv2_cases/sai_cases/processes_pathways/delay_15days'

    runs1_config = 'allActive_np4'
    runs2_config = 'passive_np4'
    sfx = '{}-{}'.format(runs1_config, runs2_config)
    
    runs1 = '{}/E3SM_ne16_L72_FIDEAL_SAI_{}/run/'.format(data, runs1_config)
    runs2 = '{}/E3SM_ne16_L72_FIDEAL_SAI_{}/run/'.format(data, runs2_config)
    
    run1 = '{}/{}_concat_hist.nc'.format(ds, runs1_config)
    run2 = '{}/{}_concat_hist.nc'.format(ds, runs2_config)
    
    # first concatenate runs to single dataset
    ctb.concat_run_outputs(runs1, outFile=run1, histnum=0, regridded=True, component='eam')
    ctb.concat_run_outputs(runs2, outFile=run2, histnum=0, regridded=True, component='eam')
    print('\n')

    #prelim_e3sm_fig(run1, run2, sfx=sfx, savedest=dest, datsavedest=ds)
    prelim_e3sm_fig(run1, run2, sfx=sfx, savedest=dest, datsavedest=ds, inj_delay=15)
    
