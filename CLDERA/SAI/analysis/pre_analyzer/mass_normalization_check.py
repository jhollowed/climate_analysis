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


def mass_norm_fig(runs, run_labels, figsavedest=None, datsavedest='.', inj_delay=0, maxtime=None):

    # read data, transform time coordinate to ndays
    print('reading data')
    dat = [xr.open_dataset(run) for run in runs]
    td = [ctb.time2day(d['time']) for d in dat]
    dat = [dat[i].assign_coords(time=td[i]) for i in range(len(dat))]
    for i in range(len(dat)):
        print('found {} timesteps from day {} to {} for run {}'.format(
                                          len(td[i]), min(td[i]), max(td[i]), run_labels[i]))
    # temp
    if(maxtime is not None):
        dat = [d.sel({'time':slice(0, maxtime)}) for d in dat]
    time = [d['time'].values for d in dat]
    
    # toggle to force re-computation
    overwrite=False
    
    # get total masses of SO2, sulfate for all time, compute analytic solution
    SO2_mass = np.zeros(len(dat), dtype=object)
    sulf_mass = np.zeros(len(dat), dtype=object)
    for i in range(len(dat)):
        try:
            if(overwrite == True): raise FileNotFoundError
            SO2_mass[i] = xr.open_dataset('{}/so2_{}.nc'.format(
                                          datsavedest, run_labels[i]))['tot_so2_mass']
            sulf_mass[i] = xr.open_dataset('{}/sulf_{}.nc'.format(
                                           datsavedest, run_labels[i]))['tot_sulf_mass']
            print('read tracers for run {}'.format(run_labels[i]))
        except FileNotFoundError:
            print('computing tracer masses for run {}...'.format(run_labels[i]))
            SO2_mass[i] = (dat[i]['SO2'] * dat[i]['AIR_MASS']).sum(('ncol', 'lev'))
            sulf_mass[i] = (dat[i]['SULFATE'] * dat[i]['AIR_MASS']).sum(('ncol', 'lev'))
            SO2_mass[i].to_dataset(name='tot_so2_mass').to_netcdf('{}/so2_{}.nc'.format(
                                                                  datsavedest, run_labels[i]))
            sulf_mass[i].to_dataset(name='tot_sulf_mass').to_netcdf('{}/sulf_{}.nc'.format(
                                                                    datsavedest, run_labels[i]))
            print('wrote tracer masses')

    # analytic total mass with V(:) = 1, sum Vk = 1
    print('computing analytic tracer masses...')
    # params
    kso2 = 1/((25 * u.day).to(u.s))
    ksulf = 1/((360*u.day).to(u.s))
    tf = (9 * u.hr).to(u.s)
    Mso2 = 1.7e10 * u.kg
    w = 2.04
    # uses time array from first run
    ttime = ((time[0] - inj_delay) * u.day).to(u.s)
    tmin = np.array([min(t, tf).m for t in ttime]) * ttime[0].u
    V = 1
    A = Mso2 / tf
    SO2_analytic = A/kso2 * V * np.exp(-kso2*ttime) * (np.exp(kso2 * tmin) - 1) 
    sulf_analytic = (w*A/((ksulf - kso2)*ksulf) * V * \
                        np.exp(-ksulf*ttime) * (\
                        kso2*(1 - np.exp(ksulf * tmin)) - ksulf*np.exp((ksulf-kso2)*ttime) * \
                        (1 - np.exp(kso2 * tmin))))
    for j in range(len(ttime)):
        if(ttime[j] < 0): SO2_analytic[j], sulf_analytic[j] = 0*u.kg,0*u.kg

    # ==========================  plotting  ==================================

    print('plotting...')
    
    label_fs = 14
    tick_fs = 12 
    fig   = plt.figure()
    ax  = fig.add_subplot(111) 

    so2c = plt.cm.plasma(np.linspace(0.15, 0.85, len(dat)))
    sulfc = plt.cm.viridis(np.linspace(0.15, 0.85, len(dat)))
  
    for i in range(len(dat)):
        ax.plot(time[i], SO2_mass[i]/1e9,  color=so2c[i], lw=3,label='SO$_2$, {}'.format(run_labels[i]))
        ax.plot(time[i], sulf_mass[i]/1e9, color=sulfc[i], lw=3,label='sulfate, {}'.format(run_labels[i]))
        print('SO2 mass {} diffs: {}'.format(run_labels[i], 
        [np.max(np.abs((SO2_mass[i] - SO2_mass[k]).values/1e9)) for k in range(len(dat))]))
        print('SULF mass {} diffs: {}'.format(run_labels[i], 
        [np.max(np.abs((sulf_mass[i] - sulf_mass[k]).values/1e9)) for k in range(len(dat))]))
    ax.plot(time[0], SO2_analytic/1e9, '--k', lw=1, label='SO$_2$, analytic')
    ax.plot(time[0], sulf_analytic/1e9,'--k',  lw=1, label='sulfate, analytic')
    ax.set_ylabel('Total mass [Tg]', fontsize=label_fs)
    ax.set_xlabel('time [days]', fontsize=label_fs)
    #ax.legend(fontsize=tick_fs, fancybox=False, loc='center right', ncol=2)
    ax.legend(fontsize=tick_fs, fancybox=False, ncol=1)
    aut.format_ticks(ax)
    ax.grid(lw=0.3, alpha=0.75)
    

    plt.tight_layout()
    if(figsavedest is not None):
        fig.savefig('{}/hswpp_totMass{}.png'.format(figsavedest,run_labels[i]), dpi=300)
        plt.show()
    else:
        plt.show() 


if(__name__ == '__main__'):

    data = '/global/cscratch1/sd/jhollo/E3SM/E3SMv2_cases/sai_cases'
    ds = '/global/cscratch1/sd/jhollo/E3SM/E3SMv2_cases/sai_cases/processes_pathways/masscheck'
    fs = '/global/homes/j/jhollo/repos/climate_analysis/CLDERA/SAI/analysis/figs/injection'
    
    #run_labels = ['allActive_delay15days', 'allActive_np4_delay15days']
    #delay=15
    #run_labels = ['passive', 'allActive', 'passive_np4', 'allActive_np4']
    #delay=0
    run_labels = ['passive', 'allActive']
    delay=0
    
    runs = [''] * len(run_labels)
    for i in range(len(run_labels)):
        runs[i] = glob.glob('{}/HSW_SAI_ne16pg2_L72_{}/run/*eam*h0*0.nc'.format(
                                                                              data, run_labels[i]))[0]
        #concatf = '{}/{}_concat_hist.nc'.format(ds, run_config)
        #ctb.concat_run_outputs(run, outFile=concatf, histnum=0, regridded=False, component='eam')
        #print('\n')

    mass_norm_fig(runs, run_labels, figsavedest=fs, datsavedest=ds, inj_delay=delay)
    
