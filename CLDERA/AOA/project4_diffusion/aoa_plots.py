import numpy as np
import os
import matplotlib.pyplot as plt
import glob
import pdb
import xarray as xr
import warnings
from climate_artist import horizontal_slice as plthor
from climate_artist import vertical_slice as pltvert
import climate_toolbox as ctb
import cftime


# =============================================================================


def time2day(time):
    start = '{}-{}-{}'.format(time[0].year, time[0].month, time[0].day)
    return cftime.date2num(time, 'days since {}'.format(start))


def aoa2clock(aoa, time):
    aoa_diff = np.zeros(aoa.shape)
    for k in range(len(time)):
        aoa_diff[k] = time[k] - aoa[k]
    return aoa_diff


# -----------------------------------------------------------------


def ssw_panels(runs, phissw=75, levssw=10, twoPanel=False, drawSSW=True, 
               histnum=0, ylim=None, xlim=None):
    '''
    Renders a plot similar to Gupta+2020 Fig 5

    Parameters
    ----------
    runs : str list
        list of CIME run directories contining the data
    phissw : float
        latitude at which to take data
    levssw : float
        pressure at which to take data
    twoPanel : bool
        whether or not to plot only the U,T and Age panels (omitting the clock tracer concentration)
    drawSSW : bool
        whether or not to plot SSW events as vertical dashed black lies
    histnum : int
        hsitoVry file group to read. Defaults to 0, in which case h0 files will
        be targeted.
    ylim : list of lists
        list of 4 length-2 lists, giving ylimits for axes 0, 1, 2, and 3
        Defaults to None, in which case it is ignored. Individual entries in this list
        may also be None.
    xlim : list
        list of length 2, giving shared xlimits for all axes
        Defaults to None, in which case it is ignored
    '''

    for i in range(len(runs)):
        
        dycore = runs[i].split('/')[-2].split('_')[0]
        grid = runs[i].split('/')[-2].split('_')[1]
        mod = runs[i].split('/')[-2].split('_')[-1]
        
        # --- compute or read data
        print('\n\n---------- working on run {} ----------'.format(runs[i].split('/')[-2]))
        dat = concat_data(runs[i], sel={'lat':phissw, 'lev':levssw}, mean=['lon'], 
                          sfx='_ssw', histnum=histnum)
        time = time2day(dat['time'].values)
        AOA2_diff = aoa2clock(dat['AOA2'], time)
        datHasNoT = False
 
        if(not twoPanel):
            fig = plt.figure(figsize=(5,7), constrained_layout=True)
            spec = fig.add_gridspec(4,1)
            ax1 = fig.add_subplot(spec[0,0])
            ax2 = fig.add_subplot(spec[1:3,0])
            ax3 = fig.add_subplot(spec[3,0])
        else:
            fig = plt.figure()
            ax1 = fig.add_subplot(211)
            ax3 = fig.add_subplot(212)
       
        ax0 = ax1.twinx()
        ax0.plot(time, dat['U'])
        try: ax1.plot(time, dat['T'], 'orange')
        except KeyError: 
            pass
            datHasNoT = True
        
        ax3.plot(time, AOA2_diff / 365, label='AOA2', color='m')
        if(histnum==0):
            ax3.plot(time, dat['AOA1'] / 365, label='AOA1', color='r')
            ax3.legend(fontsize=11, loc='lower right', frameon=False)
        
        if(not twoPanel):
            ax2.plot(time, dat['AOA2'], label='AOA2', color='r')

        if(ylim is not None):
            if(ylim[0] is not None): ax0.set_ylim(ylim[0])
            if(ylim[1] is not None): ax1.set_ylim(ylim[1])
            if(ylim[2] is not None and not twoPanel): ax2.set_ylim(ylim[2])
            if(ylim[3] is not None): ax3.set_ylim(ylim[3])
        if(xlim is not None):
            ax0.set_xlim(xlim)
            ax1.set_xlim(xlim)
            if(not twoPanel): ax2.set_xlim(xlim)
            ax3.set_xlim(xlim)
         
        try:
            try: ax1.set_ylim([np.mean(dat['T'])-40, np.mean(dat['T'])+40])
            except KeyError: pass
        except ValueError:
            pass
        try:
            ax0.set_ylim([np.median(dat['U'])-70, np.median(dat['U'])+70])
        except ValueError:
            pass

        # plot bashed lines for SSW's
        if(mod != 'mod0' and drawSSW):
            uylim = [np.mean(dat['U'])-50, np.mean(dat['U'])+50]
            for j in range(len(dat['U'])):
                if(dat['U'][j] < 0):
                    ax0.plot([time[j], time[j]], [ax0.get_ylim()[0], ax0.get_ylim()[1]], '--k', lw=0.7)
                    ax0.set_ylim([ax0.get_ylim()[0], ax0.get_ylim()[1]])
                    ax3.plot([time[j], time[j]], [ax3.get_ylim()[0], ax3.get_ylim()[1]], '--k', lw=0.7)
                    ax3.set_ylim([ax3.get_ylim()[0], ax3.get_ylim()[1]])
                    if(not twoPanel):
                        ax2.plot([time[j], time[j]], [ax2.get_ylim()[0], ax2.get_ylim()[1]], 
                                                                               '--k', lw=0.7)
                        ax2.set_ylim([ax2.get_ylim()[0], ax2.get_ylim()[1]])
        
        ax1.set_ylabel('T  [K]', fontsize=12)
        ax0.set_ylabel('u  [m/s]', fontsize=12)
        ax3.set_ylabel('Age of air  [years]', fontsize=12) 
        ax3.set_xlabel('time [years]', fontsize=12)
        ax1.axes.get_xaxis().set_ticks([])
        if(not twoPanel): 
            ax2.set_ylabel('Clock tracer concentration', fontsize=12)
            ax2.axes.get_xaxis().set_ticks([])
        if(datHasNoT):
            ax1.remove()

        ax0.set_title('{}_{} at {}deg, 10hPa'.format(dycore, grid, phissw), fontsize=14) 
        plt.savefig('figs/{}_{}_{}deg_ssw_{}_xlim{}_2p{}.png'.format(dycore, grid, 
                                      phissw, mod, xlim is None, twoPanel), dpi=300)


# =============================================================================


def aoa_profiles(runs, sfx, twoPanel=False, histnum=0, cvar='U'):
    '''
    Renders a plot similar to Gupta+2020 Fig.8, with a chosen variable on the color scale

    Parameters
    ----------
    runs : str list
        list of CIME run directories contining the data
    twoPanel : bool
        whether or not to render subpolots for both AOA1 and AOA2. If False (default), 
        plot only AOA1
    histnum : int
        hsitoVry file group to read. Defaults to 0, in which case h0 files will
        be targeted.
    cvar : str
        Variable to show on the colorscale
    '''
    
    for i in range(len(runs)):
        
        dycore = runs[i].split('/')[-2].split('_')[0]
        grid = runs[i].split('/')[-2].split('_')[1]
        mod = runs[i].split('/')[-2].split('_')[-1]
        
        # --- compute or read data
        print('\n\n---------- working on run {} ----------'.format(runs[i].split('/')[-2]))
        dat = concat_data(runs[i], sel={'time':slice('0015-01-01', '0025-01-01')}, mean=['lon'], 
                          sfx='_means', histnum=histnum)
        #time = time2day(dat['time'].values)
        #AOA2_diff = aoa2clock(dat['AOA2'], time)
        
        
        #fig = plt.figure(figsize=(10,5))
        fig = plt.figure()
        if(twoPanel):
            ax1 = fig.add_subplot(121)
            ax2 = fig.add_subplot(122)
        else:
            ax1 = fig.add_subplot(111)
 
        lat = dat['lat'].values
        lon = dat['lev'].values
        
        levels = np.linspace(1, 25, 25)
        ylim = [0.25,1000]

        var_dict = [{'var':dat['AOA1']/365, 'plotType':'contour', 'plotArgs':{'levels':levels}},
                    {'var':dat['U'], 'plotType':'contourf', 'colorArgs':{'label':'u  [m/s]'}}]
        pltvert(lat, lon, var_dict, ax=ax1, ylim=ylim)
        ax1.set_title('AOA1, 10 year mean')
        ax1.set_ylabel('lev  [hPa]', fontsize=12)
        ax1.set_xlabel('lat  [deg]', fontsize=12)
        
        if(twoPanel):
            var_dict[0]['var'] = AOA2_diff/365
            var_dict[1]['colorArgs']['ax']=ax2
            pltvert(lat, lon, var_dict, ax=ax2, ylim=ylim)
            ax2.set_title('AOA2, 10 year mean')
            ax2.set_xlabel('lat  [deg]', fontsize=12)
            
        plt.savefig('figs/sss_{}_{}_{}_{}.png'.format(dycore, grid, mod, sfx), dpi=300)





if __name__ == '__main__':            
    
    topdir = '/glade/scratch/jhollowed/CAM/cases_589/project4'
    RUNS = np.array(glob.glob('{}/*/run'.format(topdir)))


    # -------- concat all data --------
    
    if(0):
        for run in RUNS:
            run_name = run.split('/')[-2]
            dycore = run.split('/')[-2].split('_')[0]
            grid = run.split('/')[-2].split('_')[1]
            mod = run.split('/')[-2].split('_')[2]
            
            ctb.concat_resubs(run, sel={'lat':75, 'lev':10}, mean=['lon'], sfx='_ssw', 
                              histnum=0, overwrite=True, outFile='{}.CAT.SSW.nc'.format(run_name))
            ctb.concat_resubs(run, sel={'time':slice('0015-01-01', '0035-01-01')}, mean=['lon', 'time'], 
                              histnum=0, overwrite=True, outFile='{}.CAT.MEANS.nc'.format(run_name))
            pass
    
    mod3_ensRuns = RUNS[np.logical_and(['FV3_C48L72_mod3' in r for r in RUNS], 
                                       ['vt8' not in r for r in RUNS])]
    mod3_vt8_ensRuns = RUNS[['vt8' in r for r in RUNS]]
    mod4_ensRuns = RUNS[['FV3_C48L72_mod4' in r for r in RUNS]]
    SE_runs = RUNS[['SE' in r for r in RUNS]]

    mod3_vt8_means = [glob.glob('{}/*.CAT.MEANS.nc'.format(d)) for d in mod3_vt8_ensRuns]
    mod3_means = [glob.glob('{}/*.CAT.MEANS.nc'.format(d)) for d in mod3_ensRuns]
    mod4_means = [glob.glob('{}/*.CAT.MEANS.nc'.format(d)) for d in mod4_ensRuns]

    # compute ensemble mean, std

    if(0):
        ctb.ensemble_mean(mod3_means, std=True, avg_vars=['AOA1','U','V'], overwrite=True,
                          outFile = '{}/FV3_mod3_ensMean.nc'.format(topdir))
        ctb.ensemble_mean(mod3_vt8_means, std=True, avg_vars=['AOA1','U','V'], overwrite=True,
                          outFile = '{}/FV3_mod3_vt8_ensMean.nc'.format(topdir))
        ctb.ensemble_mean(mod4_means, std=True, avg_vars=['AOA1','U','V'], overwrite=True,
                          outFile = '{}/FV3_mod4_ensMean.nc'.format(topdir))

    # AOA slices for 3 SE, 3 FV3 ens means

    # SSW plots for 3 SE, 3 FV3 ens1

    # SSW tally plots for 3 FV3 ens means

    # Gupta fig 9 for 3 SE, 3 FV3 ens means




     

    #ssw_panels(baseRuns, 75, 10, False)
    #ssw_panels(SE_run, 75, 10, True, xlim = [5500,8500], ylim=[None, None,[1500,6000],[5,20]], histnum=1)
    #ssw_panels(ingredientRuns, 75, 10, True)

    #aoa_profiles(ingredientRuns, 'ingr')
    #aoa_profiles(baseRuns, 'base')
    
