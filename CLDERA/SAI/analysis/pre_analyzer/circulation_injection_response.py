import pdb
import glob
import warnings
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


def analyze_response(run1, run2, run1_native=None, run2_native=None, overwrite=False,  
                     sfx='', savedest=None, datsavedest='.', inj_delay=0):

    # windows in days over which to compoute means for horizontal slices
    # second month, fourth month
    # offset if injection_dealy > 0
    t_windows = [slice(30+inj_delay, 60+inj_delay), slice(120+inj_delay, 150+inj_delay)]
    t_labels = ['month 2\n', 'month 5\n']
    #t_labels = ['', '']

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

    # if data lengths don't match, but time frequency does, warn user and trim to match
    dat1dt = np.diff(dat1.nsteph)[0]
    dat2dt = np.diff(dat2.nsteph)[0]
    assert dat1dt == dat2dt, "dataset time sampling does not match!"\
                            "\ndata 1: samples every {} timesteps"\
                            "\ndata 2: samples every {} timesteps".format(dat1dt, dat2dt)
    trim_dat1 = False
    trim_dat2 = False
    if max(td1) > max(td2):
        warnings.warn('data lengths do not match in time; trimming data 1 to {} days'.format(max(td2)))
        dat1 = dat1.sel({'time':slice(min(td2), max(td2))})
        trim_dat1 = True
    if max(td1) < max(td2):
        warnings.warn('data lengths do not match in time; trimming data 2 to {} days'.format(max(td2)))
        dat2 = dat2.sel({'time':slice(min(td1), max(td1))})
        trim_dat2 = True

    # read native if passed
    if(run1_native is not None and run2_native is not None):
        print('reading native data 1')
        ndat1 = xr.open_dataset(run1_native)
        ndat1 = ndat1.assign_coords(time=td1)
        if(trim_dat1):
            ndat1 = ndat1.sel({'time':slice(min(td2), max(td2))})
        
        print('reading native data 2')
        ndat2 = xr.open_dataset(run2_native)
        ndat2 = ndat2.assign_coords(time=td2)
        if(trim_dat2):
            ndat2 = ndat2.sel({'time':slice(min(td2), max(td2))})

        print('using native grid data for computing total tracer masses...')
        ttdat1 = ndat1
        ttdat2 = ndat2
        ttdims = ('ncol', 'lev')
        ttstr='native'
        file_ext = 'pathways_nativeMass{}.nc'.format(sfx)
    else:
        print('run1_native AND/OR run2_native = None;')
        print('using interpolated grid data for computing total tracer masses...')
        ttdat1 = dat1
        ttdat2 = dat2
        ttdims = ('lat', 'lon', 'lev')
        ttstr=''
        file_ext = 'pathways{}.nc'.format(sfx)
   
    # update time arrays in case either datasets were trimmed above
    td1 = dat1['time'].values
    td2 = dat1['time'].values
    # assume that these are the same for both datasets
    time = td1
    lat = dat1['lat']
    lev = dat1['lev']
    
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
        Tdiff_global = xr.open_dataset('{}/Tdiff_global_{}'.format(
                                             datsavedest, file_ext))['Tdiff_global']
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
        SO21 = xr.open_dataset('{}/so21_{}'.format(datsavedest,file_ext))['tot_SO2_mass']
        sulf1 = xr.open_dataset('{}/sulf1_{}'.format(datsavedest,file_ext))['tot_SULFATE_mass']
        print('read tracers 1')
    except FileNotFoundError:
        print('computing tracers 1...')
        SO21 = (ttdat1['SO2'] * ttdat1['AIR_MASS']).sum(ttdims)
        sulf1 = (ttdat1['SULFATE'] * ttdat1['AIR_MASS']).sum(ttdims)
        SO21.to_dataset(name='tot_SO2_mass').to_netcdf('{}/so21_{}'.format(datsavedest, file_ext))
        sulf1.to_dataset(name='tot_SULFATE_mass').to_netcdf('{}/sulf1_{}'.format(
                                                            datsavedest, file_ext))
        print('wrote tracers 1')
    
    try:
        if(overwrite == True): raise FileNotFoundError
        SO22 = xr.open_dataset('{}/so22_{}'.format(datsavedest,file_ext))['tot_SO2_mass']
        sulf2 = xr.open_dataset('{}/sulf2_{}'.format(datsavedest,file_ext))['tot_SULFATE_mass']
        print('read tracers 2')
    except FileNotFoundError:
        print('computing tracers 2...')
        SO22 = (ttdat2['SO2'] * ttdat2['AIR_MASS']).sum(ttdims)
        sulf2 = (ttdat2['SULFATE'] * ttdat2['AIR_MASS']).sum(ttdims)
        SO22.to_dataset(name='tot_SO2_mass').to_netcdf('{}/so22_{}'.format(datsavedest, file_ext))
        sulf2.to_dataset(name='tot_SULFATE_mass').to_netcdf('{}/sulf2_{}'.format(
                                                                          datsavedest, file_ext))
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
    clev_H_cticks = clev_H
    clev_H_fmt = '%.2f'
    cmap_H = mpl.cm.RdYlGn_r
    
    clev_T_divnorm=colors.TwoSlopeNorm(vmin=-4, vcenter=0., vmax=8)
    clev_T = np.linspace(-4, 8, 7)
    clev_T_cticks = clev_T[clev_T % 1==0]
    clev_T_fmt = '%.0f'
    cmap_T = mpl.cm.RdYlBu_r
    
    clev_U_divnorm=colors.TwoSlopeNorm(vmin=-10, vcenter=0., vmax=15)
    clev_U = np.linspace(-10, 15, 11)
    clev_U_cticks = clev_U[clev_U % 1==0]
    clev_U_fmt = '%.0f'
    cmap_U = mpl.cm.rainbow
    
    clev_Tg_divnorm=colors.TwoSlopeNorm(vmin=-2, vcenter=0., vmax=4)
    clev_Tg = np.linspace(-2, 4, 13)
    clev_Tg_cticks = clev_Tg[clev_Tg % 1==0]
    clev_Tg_fmt = '%.0f'
    
    label_fs = 14
    tick_fs = 12
    contour_lw = 0.4

    
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
    cArgs = {'orientation':'horizontal', 'location':'top', 
             'label':'{}Heating  rate [K/day]'.format(t_labels[0]), 'aspect':20, 'format':clev_H_fmt}
    var_dict = [{'var':Hdiff.sel({'time':t_windows[0]}).mean('time'), 
                 'plotType':'contourf', 'plotArgs':pltargs, 'colorArgs':cArgs}]
    cf = pltvert(lat, lev, var_dict, ax=axH0, plot_zscale=False, grid=True)
    cf[0].set_ticks(clev_H_cticks)
    cf[0].ax.tick_params(rotation=90)
    
    # lower plot (second time window) will have no colobar; levels agree with above
    var_dict = [{'var':Hdiff.sel({'time':t_windows[1]}).mean('time'), 
        'plotType':'contourf', 'plotArgs':pltargs, 'colorArgs':cArgs}]#, 'colorFormatter':None}]
    cArgs['label'] = '{}Heating  rate [K/day]'.format(t_labels[1])
    cf = pltvert(lat, lev, var_dict, ax=axH1, plot_zscale=False,grid=True)
    cf[0].set_ticks(clev_H_cticks)
    cf[0].ax.tick_params(rotation=90)
    
    # ---------- TDIFF PLOTS ----------
    print('PLOTTING TEMP')

    pltargs = {'levels':clev_T, 'cmap':cmap_T, 'zorder':0, 'norm':clev_T_divnorm}
    cArgs = {'orientation':'horizontal', 'location':'top', 
             'label':'{}T anomaly [K]'.format(t_labels[0]), 'aspect':20, 'format':clev_T_fmt}
    var_dict = [{'var':Tdiff.sel({'time':t_windows[0]}).mean('time'), 
                 'plotType':'contourf', 'plotArgs':pltargs, 'colorArgs':cArgs}]
    cf = pltvert(lat, lev, var_dict, ax=axT0, plot_zscale=False,grid=True)
    cf[0].set_ticks(clev_T_cticks)
    cf[0].ax.tick_params(rotation=90)
    
    # lower plot (second time window) will have no colobar; levels agree with above
    var_dict = [{'var':Tdiff.sel({'time':t_windows[1]}).mean('time'), 
        'plotType':'contourf', 'plotArgs':pltargs, 'colorArgs':cArgs}]#, 'colorFormatter':None}]
    cArgs['label'] = '{}T anomaly [K]'.format(t_labels[1])
    cf = pltvert(lat, lev, var_dict, ax=axT1, plot_zscale=False,grid=True) 
    cf[0].set_ticks(clev_T_cticks)
    cf[0].ax.tick_params(rotation=90)
 
    # ---------- UDIFF PLOTS ----------
    print('PLOTTING WIND')
    
    pltargs = {'levels':clev_U, 'cmap':cmap_U, 'zorder':0, 'norm':clev_U_divnorm}
    cArgs = {'orientation':'horizontal', 'location':'top', 
             'label':'{}U anomaly [m/s]'.format(t_labels[0]), 'aspect':20, 'format':clev_U_fmt}
    var_dict = [{'var':Udiff.sel({'time':t_windows[0]}).mean('time'), 
                 'plotType':'contourf', 'plotArgs':pltargs, 'colorArgs':cArgs}]
    cf = pltvert(lat, lev, var_dict, ax=axU0, plot_zscale=False, grid=True)
    cf[0].set_ticks(clev_U_cticks)
    cf[0].ax.tick_params(rotation=90)
    
    # lower plot (second time window) will have no colobar; levels agree with above
    var_dict = [{'var':Udiff.sel({'time':t_windows[1]}).mean('time'), 
        'plotType':'contourf', 'plotArgs':pltargs, 'colorArgs':cArgs}]#, 'colorFormatter':None}]
    cArgs['label'] = '{}U anomaly [m/s]'.format(t_labels[1])
    cf = pltvert(lat, lev, var_dict, ax=axU1, plot_zscale=False, grid=True)
    cf[0].set_ticks(clev_U_cticks)
    cf[0].ax.tick_params(rotation=90)
    
    
    # ---------- global T diff plot ----------
    print('PLOTTING GLOBAL TEMP')
   
    TIME, LEV = np.meshgrid(time, lev)
    cf = axTg.contourf(TIME, LEV, Tdiff_global.T, cmap=cmap_T, levels=clev_Tg, extend='both', 
                       norm=clev_Tg_divnorm)
    cfl = axTg.contour(TIME, LEV, Tdiff_global.T, levels=clev_Tg, colors='k', linewidths=contour_lw, 
                       zorder=1)
    cb = fig.colorbar(cf, ax=axTg, orientation='horizontal', location='top',
                      aspect=40, extend='both', extendrect=False, format=clev_Tg_fmt)
    cb.set_label('Global mean T anomaly [K]', fontsize=tick_fs)
    #axTg.clabel(cfl, levels=clev_Tg_cticks, inline=True, fmt=clev_Tg_fmt, fontsize=tick_fs)
    cfl.collections[np.where(clev_Tg == 0)[0][0]].set_linewidth(contour_lw*1.66)
    cb.set_ticks(clev_Tg_cticks)

    axTg.set_yscale('log')
    axTg.invert_yaxis()
    axTg.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y,pos: \
                                   ('{{:.{:1d}f}}'.format(int(np.maximum(-np.log10(y),0)))).format(y)))
    axTg.set_ylabel('pressure [hPa]', fontsize=label_fs)

    # TMP VIEW SURF
    if(0):
        axTg.set_yscale('linear')
        axTg.set_ylim([max(lev), 900])
    
    
    # ---------- global T diff plot at levels ----------
    print('PLOTTING GLOBAL NLEV TEMP')
   
    cls = plt.cm.viridis(np.linspace(0.2, 0.8, 4))
    #Tg_lev = [3, 30, 100, 850]
    Tg_lev = [50, 30, 100]
    for i in range(len(Tg_lev)):
        axTgl.plot(time, Tdiff_global.sel({'lev':Tg_lev[i]}, method='nearest'), 
                   label='{} hPa'.format(Tg_lev[i]), color=cls[i])
    axTgl.set_ylabel('Global Mean\nT anomaly [K]', fontsize=label_fs)
    axTgl.plot(axTgl.get_xlim(), [0, 0], '--k', lw=2)
    axTgl.legend(fontsize=tick_fs, fancybox=False, loc='upper right', ncol=2)
    
    # ---------- total tracer masses ----------
    print('PLOTTING TRACER MASSES')
   
    axTr.plot(time, SO21/1e9, '-r',  label='SO$_2$, active')
    #axTr.plot(time, SO22/1e9, '--r', label='SO$_2$, passive')
    axTr.plot(time, sulf2/1e9,'-', color='orange',label='sulfate, active')
    #axTr.plot(time, sulf2/1e9,'--', color='orange', label='sulfate, passive')
    axTr.set_ylabel('Total mass [Tg]', fontsize=label_fs)
    axTr.set_xlabel('time [days]', fontsize=label_fs)
    axTr.legend(fontsize=tick_fs, fancybox=False, loc='upper right', ncol=2)
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
        ax.set_xlim([0, np.max(time)])
        aut.format_ticks(ax)
        ax.grid(lw=0.3, alpha=0.75)
     
    plt.tight_layout()
    if(savedest is not None):
        fig.savefig('{}/heating_response{}.png'.format(savedest, sfx), dpi=300)
        plt.show()
    else: 
        plt.show()
        





if(__name__ == '__main__'):

    fig_dest = '/global/homes/j/jhollo/repos/climate_analysis/CLDERA/SAI/analysis/pre_analyzer/figs/heating_response_1.5year'
    data_dest = '/global/cscratch1/sd/jhollo/E3SM/E3SMv2_cases/sai_cases/processes_pathways_1.5year'
    data_source = '/global/cscratch1/sd/jhollo/E3SM/E3SMv2_cases/sai_cases'

    # for ens
    data_source_ens = '/global/cscratch1/sd/jhollo/E3SM/E3SMv2_cases/sai_cases/gamma_clones'
    histnum=2

    #gamma = [4, 8, 10]
    #gamma = [11, 12, 13]
    gamma = [10]
    for g in gamma:
        print('\n---------- GAMMA = {}'.format(g))
        
        runs1_config = 'allActive__cldera_sai_gammaq_{0}__cldera_sai_gammatau_{0}'.format(g)
        runs2_config = 'passive'
        #sfx = '_{}-{}'.format(runs1_config, runs2_config)
        sfx = '_gamma{}'.format(g)
        
        rundir1 = '{}/HSW_SAI_ne16pg2_L72_{}/run/'.format(data_source_ens, runs1_config)
        rundir2 = '{}/HSW_SAI_ne16pg2_L72_{}/run/'.format(data_source, runs2_config)
        runs1 = glob.glob('{}/*eam.h{}*regrid*aave*nc'.format(rundir1, histnum))
        runs2 = glob.glob('{}/*eam.h{}*regrid*aave*nc'.format(rundir2, histnum))
        
        if(len(runs1) > 1):
            run1 = '{}/{}_concat_hist.nc'.format(data_dest, runs1_config)
            ctb.concat_run_outputs(rundir1, outFile=run1, histnum=histnum, regridded=True, component='eam')
        else: 
            run1 = runs1[0]
        if(len(runs2) > 1):
            run2 = '{}/{}_concat_hist.nc'.format(data_dest, runs2_config)
            ctb.concat_run_outputs(rundir2, outFile=run2, histnum=histnum, regridded=True, component='eam')
        else: 
            run2 = runs2[0]
        print('\n')
        
        use_native_mass = True
        if(use_native_mass):
            runs1_native = glob.glob('{}/*eam.h{}.*0.nc'.format(rundir1, histnum))
            runs2_native = glob.glob('{}/*eam.h{}.*0.nc'.format(rundir2, histnum))
            if(len(runs1_native) > 1):
                run1_native = '{}/{}_concat_hist.nc'.format(data_dest, runs1_config)
                ctb.concat_run_outputs(rundir1, outFile=run1_native, 
                                       histnum=0, regridded=False, component='eam')
            else: 
                run1_native = runs1_native[0]
            if(len(runs2_native) > 1):
                run2_native = '{}/{}_concat_hist.nc'.format(data_dest, runs2_config)
                ctb.concat_run_outputs(rundir2, outFile=run2_native, 
                                       histnum=0, regridded=False, component='eam')
            else: 
                run2_native = runs2_native[0]
        else:
            runs1_native = None
            runs2_native = None
        print('\n')

        analyze_response(run1, run2, run1_native, run2_native, overwrite=False,
                         sfx=sfx, savedest=fig_dest, datsavedest=data_dest, inj_delay=0) 
