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


def heating_movie(run1, run2, run1_native=None, run2_native=None, overwrite=False,  
                  sfx='', savedest=None, datsavedest='.', inj_delay=0):

    # windows in days over which to compoute means for horizontal slices
    # second month, fourth month
    # offset if injection_dealy > 0
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

    # read native if passed
    if(run1_native is not None and run2_native is not None):
        print('reading native data 1')
        ndat1 = xr.open_dataset(run1_native)
        ndat1 = ndat1.assign_coords(time=td1)
        
        print('reading native data 2')
        ndat2 = xr.open_dataset(run2_native)
        ndat2 = ndat2.assign_coords(time=td2)

        print('using native grid data for computing total tracer masses...')
        ttdat1 = ndat1
        ttdat2 = ndat2
        ttdims = ('ncol', 'lev')
        ttstr='native'
        file_ext = 'pathways_nativeMass_day{}.nc'.format(sfx)
    else:
        print('run1_native AND/OR run2_native = None;')
        print('using interpolated grid data for computing total tracer masses...')
        ttdat1 = dat1
        ttdat2 = dat2
        ttdims = ('lat', 'lon', 'lev')
        ttstr=''
        file_ext = 'pathways_day{}.nc'.format(sfx)
    
    # assume that these are the same for both datasets
    time = td1 * u.day
    t0 = time.min()
    tf = time.max()
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
        Vdiff = xr.open_dataset('{}/Vdiff_{}'.format(datsavedest, file_ext))['Vdiff']
        print('read Vdiff')
    except FileNotFoundError:
        print('computing vdiff...')
        # ---- diff of zonal means
        V1   = dat1['V'].mean('lon')
        V2   = dat2['V'].mean('lon')
        Vdiff = V1 - V2
        Vdiff.to_dataset(name='Vdiff').to_netcdf('{}/Vdiff_{}'.format(datsavedest, file_ext))
        print('wrote vdiff')

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
    ttime = time.to(u.s)
    tmin = np.array([min(t, tf).m for t in ttime]) * ttime[0].u
    V = 1
    A = Mso2 / tf
    SO2_analytic = A/kso2 * V * np.exp(-kso2*ttime) * (np.exp(kso2 * tmin) - 1) 
    sulf_analytic = (w*A/((ksulf - kso2)*ksulf) * V * \
                        np.exp(-ksulf*ttime) * (\
                        kso2*(1 - np.exp(ksulf * tmin)) - ksulf*np.exp((ksulf-kso2)*ttime) * \
                        (1 - np.exp(kso2 * tmin))))
        
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
    
    clev_V_divnorm=colors.TwoSlopeNorm(vmin=-10, vcenter=0., vmax=30)
    clev_V = np.linspace(-10, 30, 9)
    clev_V_fmt = '%.0f'
    clev_V_labels = []
    cmap_V = mpl.cm.rainbow
    
    # clev_Tg = np.linspace(-3, 4, 8)  #old
    clev_Tg = np.linspace(-3, 6, 10)
    clev_Tg_fmt = '%.0f'
    clev_Tg_labels = []
    
    label_fs = 14
    tick_fs = 12

    
    # first 2 rows:
    # heating, T anom, U anom, for 2 30-day windows
    # last 2 rows:
    # T anom vs t, T anom vs t at select levels
    fig   = plt.figure(figsize=(14, 9))
    gs    = fig.add_gridspec(nrows=4, ncols=4)

    
    # ===================  loop over time =============================

    for j in range(len(time)):
        if j == 0: continue # skip IC
        ti = time[j]
        timei = time[time <= ti]
        tslice = slice(t0.m, ti.m)

        # ==========================  plotting  ==================================
        
        print('\n==================== TIME {}/{} ===================='.format(ti.m, tf.m))

        # clear figure
        fig.clf()
        axH  = fig.add_subplot(gs[0:2, 0])
        axT  = fig.add_subplot(gs[0:2, 1])
        axU  = fig.add_subplot(gs[2:,  0])
        axV  = fig.add_subplot(gs[2:,  1])
        axTg  = fig.add_subplot(gs[0:2, 2:])
        axTgl = fig.add_subplot(gs[2,   2:])
        axTr  = fig.add_subplot(gs[3,   2:])

        
        # ---------- HDIFF PLOTS ----------
        print('PLOTTING HEATING')

        pltargs = {'levels':clev_H, 'cmap':cmap_H, 'zorder':0, 'norm':clev_H_divnorm}
        cArgs = {'orientation':'horizontal', 'location':'top', 'aspect':20, 
                 'label':'{}Heating  rate [K/day]'.format(t_labels[0]), 'format':clev_H_fmt}
        var_dict = [{'var':Hdiff.sel({'time':ti}, method='nearest'), 
                     'plotType':'contourf', 'plotArgs':pltargs, 'colorArgs':cArgs}]
        cf = pltvert(lat, lev, var_dict, ax=axH, plot_zscale=False, grid=True)
        cf[0].set_ticks(clev_H)
        cf[0].ax.tick_params(rotation=90)
        
        
        # ---------- TDIFF PLOTS ----------
        print('PLOTTING TEMP')

        pltargs = {'levels':clev_T, 'cmap':cmap_T, 'zorder':0}
        cArgs = {'orientation':'horizontal', 'location':'top', 
                 'label':'{}T anomaly [K]'.format(t_labels[0]), 'aspect':20, 'format':clev_T_fmt}
        var_dict = [{'var':Tdiff.sel({'time':ti}, method='nearest'), 
                     'plotType':'contourf', 'plotArgs':pltargs, 'colorArgs':cArgs}]
        cf = pltvert(lat, lev, var_dict, ax=axT, plot_zscale=False,grid=True)
        cf[0].set_ticks(clev_T)
        cf[0].ax.tick_params(rotation=90)    
     
        # ---------- UDIFF PLOTS ----------
        print('PLOTTING U')
        
        pltargs = {'levels':clev_U, 'cmap':cmap_U, 'zorder':0, 'norm':clev_U_divnorm}
        cArgs = {'orientation':'horizontal', 'location':'top', 
                 'label':'{}U anomaly [m/s]'.format(t_labels[0]), 'aspect':20, 'format':clev_U_fmt}
        var_dict = [{'var':Udiff.sel({'time':ti}, method='nearest'), 
                     'plotType':'contourf', 'plotArgs':pltargs, 'colorArgs':cArgs}]
        cf = pltvert(lat, lev, var_dict, ax=axU, plot_zscale=False,grid=True)
        cf[0].set_ticks(clev_U)
        cf[0].ax.tick_params(rotation=90)
        
        # ---------- VDIFF PLOTS ----------
        print('PLOTTING V')
        
        pltargs = {'levels':clev_V, 'cmap':cmap_V, 'zorder':0, 'norm':clev_V_divnorm}
        cArgs = {'orientation':'horizontal', 'location':'top', 
                 'label':'{}V anomaly [m/s]'.format(t_labels[0]), 'aspect':20, 'format':clev_U_fmt}
        var_dict = [{'var':Vdiff.sel({'time':ti}, method='nearest'), 
                     'plotType':'contourf', 'plotArgs':pltargs, 'colorArgs':cArgs}]
        cf = pltvert(lat, lev, var_dict, ax=axV, plot_zscale=False,grid=True)
        cf[0].set_ticks(clev_V)
        cf[0].ax.tick_params(rotation=90)
        
        
        # ---------- global T diff plot ----------
        print('PLOTTING GLOBAL TEMP')
        
        TIME, LEV = np.meshgrid(timei.m, lev)
        Tdiff_dat = Tdiff_global.sel({'time':tslice})
        cf = axTg.contourf(TIME, LEV, Tdiff_dat.T, cmap=cmap_T, levels=clev_Tg, extend='both')
        cfl = axTg.contour(TIME, LEV, Tdiff_dat.T, levels=clev_Tg, colors='k', linewidths=0.6,zorder=1)
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
        Tg_lev = [3, 30, 100, 850]
        for i in range(len(Tg_lev)):
            axTgl.plot(timei, Tdiff_dat.sel({'lev':Tg_lev[i]}, method='nearest'), 
                       label='{} hPa'.format(Tg_lev[i]), color=cls[i])
        axTgl.set_ylabel('Global Mean\nT anomaly [K]', fontsize=label_fs)
        axTgl.plot(axTgl.get_xlim(), [0, 0], '--k', lw=2)
        axTgl.legend(fontsize=tick_fs, fancybox=False, loc='lower right', ncol=4)
        
        # ---------- total tracer masses ----------
        print('PLOTTING TRACER MASSES')
        
        SO21i = SO21.sel({'time':tslice})
        SO22i = SO22.sel({'time':tslice})
        sulf1i = sulf1.sel({'time':tslice})
        sulf2i = sulf2.sel({'time':tslice})

        axTr.plot(timei, SO21i/1e9, '-r',  label='SO$_2$, active')
        #axTr.plot(timei, SO22i/1e9, '--r', label='SO$_2$, passive')
        axTr.plot(timei, sulf1i/1e9,'-', color='orange',label='sulfate, active')
        #axTr.plot(timei, sulf2i/1e9,'--', color='orange', label='sulfate, passive')
        axTr.set_ylabel('Total mass [Tg]', fontsize=label_fs)
        axTr.set_xlabel('time [days]', fontsize=label_fs)
        axTr.legend(fontsize=tick_fs, fancybox=False, loc='lower right', ncol=2)
        #axTr.set_ylim([0, 30])
        
        # -------- formatting --------
        zonal_axes = [axH, axT, axU, axV]
        summ_axes = [axTg, axTgl, axTr]
        top_axes = [axH, axT, axTg, axTgl]
        right_axes = [axT, axV]
        
        for ax in zonal_axes:
            ax.set_ylabel('pressure  [hPa]')
            ax.set_xlabel('lat [deg]')
        
        for ax in zonal_axes:
            aut.format_ticks(ax)
        for ax in top_axes:
            ax.set_xlabel('')
        for ax in right_axes:
            ax.set_ylabel('')
        for ax in summ_axes:
            aut.format_ticks(ax)
            ax.set_xlim([t0, tf])
            aut.format_ticks(ax)
            ax.grid(lw=0.3, alpha=0.75)
        axTr.set_ylim([-2, 33])
        axTgl.set_ylim([-4, 8])
        
        plt.tight_layout()
        if(savedest is not None):
            fig.savefig('{}/heating_response{}_frame{:04d}.png'.format(savedest, sfx, j), dpi=100)
        else: 
            plt.show()
        





if(__name__ == '__main__'):

    fig_dest = '/global/homes/j/jhollo/repos/climate_analysis/CLDERA/SAI/analysis/figs/'\
               'pre_logcool/heating_response_1year_movie'
    data_dest = '/global/cscratch1/sd/jhollo/E3SM/E3SMv2_cases/sai_cases/processes_pathways_1year'
    data_source = '/global/cscratch1/sd/jhollo/E3SM/E3SMv2_cases/sai_cases/pre_logcool'
    histnum=2

    runs1_config = 'allActive'
    runs2_config = 'passive'
    sfx = '_{}-{}'.format(runs1_config, runs2_config)
    
    rundir1 = '{}/HSW_SAI_ne16pg2_L72_{}/run/'.format(data_source, runs1_config)
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

    heating_movie(run1, run2, run1_native, run2_native, overwrite=False,
                  sfx=sfx, savedest=fig_dest, datsavedest=data_dest, inj_delay=0) 
