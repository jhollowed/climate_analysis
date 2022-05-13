import numpy as np
import os
import sys
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


def ssw_panels(inFile, phissw=75, levssw=10, twoPanel=False, drawSSW=True, 
               ylim=None, xlim=None):
    '''
    Renders a plot similar to Gupta+2020 Fig 5

    Parameters
    ----------
    inFile : str
        location of file containing the data
    phissw : float
        latitude at which to take data
    levssw : float
        pressure at which to take data
    twoPanel : bool
        whether or not to plot only the U,T and Age panels (omitting the clock tracer concentration)
    drawSSW : bool
        whether or not to plot SSW events as vertical dashed black lies
    ylim : list of lists
        list of 4 length-2 lists, giving ylimits for axes 0, 1, 2, and 3
        Defaults to None, in which case it is ignored. Individual entries in this list
        may also be None.
    xlim : list
        list of length 2, giving shared xlimits for all axes
        Defaults to None, in which case it is ignored
    '''
    

    run_name = inFile.split('.nc')[0].split('/')[-1]
    dycore = run_name.split('_')[0]
    grid = run_name.split('_')[1]
    
    nl_names = []
    nl_vals = []
    if('__' in run_name):
        run_name = run_name.split('__RESUBMIT')[0]
        if('se_nu' in run_name):
            #nl_names.append('ν')
            #nl_vals.append('{:.2e}'.format(float(run_name.split('mod3__se_nu_')[-1].split('_')[0])))
            #nl_names.append('ν_div')
            #nl_vals.append('{:.2e}'.format(float(run_name.split('se_nu_div_')[-1].split('.CAT')[0])))
            nl_names.append('config')
            nl_vals.append(run_name.split(grid)[-1])
        if(dycore == 'FV3'):
            nl_names.append('hord_vt')
            nl_names.append('RF heating')
            if('vt' in run_name): nl_vals.append('8')
            else: nl_vals.append('10')
            if('mod3' in run_name): nl_vals.append('off')
            if('mod4' in run_name): nl_vals.append('on')
            ens = ' {}'.format(run_name.split('__')[-1].split('_')[0])
        else:
            ens = ''
    nl_str = ''
    for i in range(len(nl_names)):
        if(i == len(nl_names)-1): pp=''
        else: pp=', '
        nl_str = nl_str + '{}={}{}'.format(nl_names[i], nl_vals[i],pp)

    # --- read data
    print('\n\n---------- working on run {} ----------'.format(run_name))
    dat = xr.open_dataset(inFile)
    time = time2day(dat['time'].values)
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
    
    ax3.plot(time, dat['AOA1'] / 365, label='AOA1', color='r') 
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

    # plot dashed lines for SSW's
    uylim = [np.mean(dat['U'])-50, np.mean(dat['U'])+50]
    Tylim = [np.mean(dat['T'])-50, np.mean(dat['T'])+50]
    aoalim = ax3.get_ylim()
    if(not twoPanel):
        clocklim=ax2.get_ylim()

    if(drawSSW):
        for j in range(len(dat['U'])):
            if(dat['U'][j] < 0):
                ax0.plot([time[j], time[j]], [uylim[0], uylim[1]], '--k', lw=0.7)
                ax3.plot([time[j], time[j]], [aoalim[0], aoalim[1]], '--k', lw=0.7)
                if(not twoPanel):
                    ax2.plot([time[j], time[j]], [clocklim[0], clocklim[1]], '--k', lw=0.7)
    ax0.set_ylim([uylim[0], uylim[1]])
    ax1.set_ylim([Tylim[0], Tylim[1]])
    ax3.set_ylim([aoalim[0], aoalim[1]])
    if(not twoPanel): ax2.set_ylim([clocklim[0], clocklim[1]])
    
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

    ax0.set_title('{}_{}{} at {}deg, 10hPa\n{}'.format(dycore, grid, ens, phissw, nl_str), fontsize=14)
    plt.savefig('./figs/{}_{}deg_ssw.png'.format(run_name, phissw), dpi=300)


# =============================================================================


def aoa_profiles(inFile, twoPanel=False, cvar='U', cvar_label = 'u  [m/s]', concatFile=None, sfx=None):
    '''
    Renders a plot similar to Gupta+2020 Fig.8, with a chosen variable on the color scale

    Parameters
    ----------
    inFile : str
        location of file containing the data
    twoPanel : bool
        whether or not to render subpolots for both AOA1 and AOA2. If False (default), 
        plot only AOA1
    cvar : str or xarray DataArray
        Variable to show on the colorscale. If from data. If DataArray,
        plot this.
    cvar_label : str
        Label for colorbar
    sfx : str
        String to add as suffix to figure filename
    '''
    
    # --- extract metadata
    run_name = inFile.split('.nc')[0].split('/')[-1]
    dycore = run_name.split('_')[0]
    grid = run_name.split('_')[1]
    
    nl_names = []
    nl_vals = []
    if('__' in run_name):
        run_name = run_name.split('__RESUBMIT')[0]
        if('se_nu' in run_name):
            #nl_names.append('ν')
            #nl_vals.append('{:.2e}'.format(float(run_name.split('mod3__se_nu_')[-1].split('_')[0])))
            #nl_names.append('ν_div')
            #nl_vals.append('{:.2e}'.format(float(run_name.split('se_nu_div_')[-1].split('.CAT')[0])))
            nl_names.append('config')
            nl_vals.append(run_name.split(grid)[-1])
        if(dycore == 'FV3'):
            nl_names.append('hord_vt')
            nl_names.append('RF heating')
            if('vt' in run_name): nl_vals.append('8')
            else: nl_vals.append('10')
            if('mod3' in run_name): nl_vals.append('off')
            if('mod4' in run_name): nl_vals.append('on')
            ens = ' ensemble mean'
        else:
            ens = ''
    else:
        ens = ''

    nl_str = ''
    for i in range(len(nl_names)):
        if(i == len(nl_names)-1): pp=''
        else: pp=', '
        nl_str = nl_str + '{}={}{}'.format(nl_names[i], nl_vals[i],pp)
    
    # --- compute or read data
    print('\n\n---------- working on run {} ----------'.format(run_name))
    dat = xr.open_dataset(inFile)
    
    fig = plt.figure()
    if(twoPanel):
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)
    else:
        ax1 = fig.add_subplot(111)

    lat = dat['lat'].values
    lev = dat['lev'].values
    
    levels = np.linspace(1, 25, 25)
    ylim = [0.25,1000]

    if(isinstance(cvar, str)):
        color_var = dat[cvar]
    else:
        color_var = cvar

    var_dict = [{'var':dat['AOA1']/365, 'plotType':'contour', 'plotArgs':{'levels':levels}}]
    var_dict_color= [{'var':color_var, 'plotType':'contourf', 'colorArgs':{'label':cvar_label}}]
    pltvert(lat, lev, var_dict, ax=ax1, ylim=ylim)
    pltvert(color_var['lat'], color_var['lev'], var_dict_color, ax=ax1, ylim=ylim)
    ax1.set_title('AOA1, 20 year mean')
    ax1.set_ylabel('lev  [hPa]', fontsize=12)
    ax1.set_xlabel('lat  [deg]', fontsize=12)
    
    if(twoPanel):
        AOA2_diff = aoa2clock(dat['AOA2'], time)
        var_dict[0]['var'] = AOA2_diff/365
        var_dict[1]['colorArgs']['ax']=ax2
        pltvert(lat, lev, var_dict, ax=ax2, ylim=ylim)
        ax2.set_title('AOA2, 10 year mean')
        ax2.set_xlabel('lat  [deg]', fontsize=12)
    
    ax1.set_title('{}_{}{}\n{}'.format(dycore, grid, ens, nl_str), fontsize=14)
    plt.savefig('./figs/aoaSlice_{}{}.png'.format(run_name, sfx), dpi=300)


# =============================================================================


def aoa_levels(inFiles, errFiles, colors=None, sfx=None):
    '''
    Renders a plot similar to Gupta+2020 Fig.9

    Parameters
    ----------
    inFiles : str
        location of files containing the data
    errFiles : str
        location of files containing the spread of the data. Individual entries can be None
    colors : list
        list of line colors
    sfx : str
        String to add as suffix to figure filename
    '''
   
    #fig = plt.figure(figsize=(5,10))
    #ax1 = fig.add_subplot(311)
    #ax2 = fig.add_subplot(312)
    #ax3= fig.add_subplot(313)
    fig = plt.figure(figsize=(10,5))
    ax1 = fig.add_subplot(131)
    ax2 = fig.add_subplot(132)
    ax3= fig.add_subplot(133)
    axes = [ax1, ax2, ax3]

    if(colors is None):
        colors = ['g', 'b', 'c', 'k', 'r', 'orange', 'm']
    linestyles = {'FV3':'-', 'SE':'--'}
    labels=['SE increased total wind diffusion', '' ,'SE defaults', 'SE increased divergence damping', 'FV3 RF heating = off, hord_vt=8', 'FV3 RF heating = off, hord_vt=10', 'FV3 RF heating = on, hord_vt=8']

    for i in range(len(inFiles)):
        if(i == 1): continue
        inFile = inFiles[i]
   
        # --- extract metadata
        run_name = inFile.split('.nc')[0].split('/')[-1]
        dycore = run_name.split('_')[0]
        grid = run_name.split('_')[1]
        
        nl_names = []
        nl_vals = []
        if('__' in run_name or 'vt' in run_name):
            run_name = run_name.split('__RESUBMIT')[0]
            if('se_nu' in run_name):
                #nl_names.append('ν')
                #nl_vals.append('{:.2e}'.format(float(run_name.split('mod3__se_nu_')[-1].split('_')[0])))
                #nl_names.append('ν_div')
                #nl_vals.append('{:.2e}'.format(float(run_name.split('se_nu_div_')[-1].split('.CAT')[0])))
                nl_names.append('config')
                nl_vals.append(run_name.split(grid)[-1])
            if(dycore == 'FV3'):
                nl_names.append('hord_vt')
                nl_names.append('RF heating')
                if('vt' in run_name): nl_vals.append('8')
                else: nl_vals.append('10')
                if('mod3' in run_name): nl_vals.append('off')
                if('mod4' in run_name): nl_vals.append('on')
                ens = ' ensemble mean'
            else:
                ens = ''
        nl_str = ''
        for k in range(len(nl_names)):
            if(k == len(nl_names)-1): pp=''
            else: pp=', '
            nl_str = nl_str + '{}={}{}'.format(nl_names[k], nl_vals[k],pp)

        
        # --- compute or read data
        print('\n\n---------- working on run {} ----------'.format(run_name))
        dat = xr.open_dataset(inFile)
        if(errFiles[i] is not None): 
            err = xr.open_dataset(errFiles[i])
        else:
            err = None
        
        
        lat = dat['lat'].values
        AOA_lev1 = dat['AOA1'].sel({'lev':10}, method='nearest')
        AOA_lev2 = dat['AOA1'].sel({'lev':30}, method='nearest')
        AOA_lev3 = dat['AOA1'].sel({'lev':150}, method='nearest')
        AOA_allLev = [AOA_lev1, AOA_lev2, AOA_lev3]
        if(err is not None):
            AOA_lev1_err = err['AOA1'].sel({'lev':10}, method='nearest')
            AOA_lev2_err = err['AOA1'].sel({'lev':30}, method='nearest')
            AOA_lev3_err = err['AOA1'].sel({'lev':150}, method='nearest')
            AOA_allErr = [AOA_lev1_err, AOA_lev2_err, AOA_lev3_err]
        
        print('PLOTTING {} with COLOR {} at i={}'.format(inFile.split('/')[-1], colors[i], i))
        for j in range(len(axes)):
            if(err is not None):
                age_err = [AOA_allLev[j] + AOA_allErr[j], AOA_allLev[j] - AOA_allErr[j]]
                axes[j].fill_between(lat, age_err[0]/365, age_err[1]/365, color=colors[i], alpha=0.33)
            line = axes[j].plot(lat, AOA_allLev[j]/365, ls=linestyles[dycore], color=colors[i])[0]
            if(j == 0):             
                #line.set_label('{} {}'.format(dycore, nl_str))
                line.set_label(labels[i])

    ax1.set_title('AOA at different pressure levels, 20 year mean')
    
    ax1.set_ylabel('age  [years]', fontsize=12)
    ax1.annotate('10 hPa', xy=(0.95, 0.05), xycoords='axes fraction', fontsize=10,
                          horizontalalignment='right', verticalalignment='bottom')
    ax1.tick_params(labelbottom=False)
    ax1.set_xlim([-90, 90])
    lgd=ax1.legend(fontsize=9, bbox_to_anchor=(0, 1.1, 1, 0.2), loc='lower left', mode='expand')

    ax2.set_ylabel('age  [years]', fontsize=12)
    ax2.annotate('30 hPa', xy=(0.95, 0.05), xycoords='axes fraction', fontsize=10,
                           horizontalalignment='right', verticalalignment='bottom')
    ax2.tick_params(labelbottom=False)
    ax2.set_xlim([-90, 90])
    
    ax3.set_ylabel('age  [years]', fontsize=12)
    ax3.annotate('150 hPa', xy=(0.95, 0.05), xycoords='axes fraction', fontsize=10,
                            horizontalalignment='right', verticalalignment='bottom')
    ax3.set_xlabel('lat  [deg]', fontsize=12)
    ax3.set_xlim([-90, 90])

    plt.tight_layout()

    plt.show()
    plt.savefig('./figs/aoaLevels{}.png'.format(sfx), dpi=300, 
                bbox_extra_artists=(lgd,), bbox_inches='tight')

    
# =============================================================================


def ssw_counts(inFiles, colors=None, sfx=None):
    '''
    Renders a plot similar to Gupta+2020 Fig.6

    Parameters
    ----------
    inFiles : list
        location of files containing the data. elements can be list of strings, in which case
        that element will compute a single number of ssw events, with associated spread.
    colors : list
        list of bar colors
    sfx : str
        String to add as suffix to figure filename
    '''
   
    fig = plt.figure()
    ax1 = fig.add_subplot(111)

    if(colors is None):
        colors = ['r', 'o', 'm', 'g', 'b', 'c'] 
    width = 0.25
        
    for i in range(len(inFiles)):
        inFile = inFiles[i]
        nl_names = []
        nl_vals = []

        # is SE run
        if(isinstance(inFile, str)):
            run_name = inFile.split('.nc')[0].split('/')[-1]
            run_name = run_name.split('__RESUBMIT')[0]
            nl_names.append('ν')
            nl_vals.append('{:.2e}'.format(float(run_name.split('mod3__se_nu_')[-1].split('_')[0])))
            nl_names.append('ν_div')
            nl_vals.append('{:.2e}'.format(float(run_name.split('se_nu_div_')[-1].split('.CAT')[0])))
            dycore = 'SE'
            
            sswc = ctb.ssw_events(inFile)
            err = 0
        
        # is FV3 ensemble
        else:
            run_name = inFile.split('.nc')[0].split('/')[-1]
            nl_names.append('hord_vt')
            nl_names.append('RF heating')
            if('vt' in run_name): nl_vals.append('8')
            else: nl_vals.append('10')
            if('mod3' in run_name): nl_vals.append('off')
            if('mod4' in run_name): nl_vals.append('on')
            dycore='FV3'
            
            sswc_all = np.array([ctb.ssw_events(inf) for inf in inFile])
            sswc = np.mean(sswc_all)
            err = np.std(sswc_all)
             
        nl_str = ''
        for i in range(len(nl_names)):
            if(i == len(nl_names)-1): pp=''
            else: pp=', '
            nl_str = nl_str + '{}={}{}'.format(nl_names[i], nl_vals[i],pp)
        
        label = '{}, {}'.format(dycore, nl_str)
        
        
        time = time2day(dat['time'].values)
        total_time = int((time[-1] - time[0])/365) # total time length of averaging window in years
        
        ax1.bar(colors[i], sswc, yerr=err, label=label) 

    ax1.set_title('Average number of SSW events per {} years'.format(total_time))
    
    ax1.set_ylabel('Number of events', fontsize=12)
    ax1.get_xaxis().set_ticks([])
    ax1.get_xaxis().set_visible(False)

    plt.savefig('./figs/sswCounts{}.png'.format(sfx), dpi=300)


# ============================================================================================


if __name__ == '__main__':            
    
    topdir = '/glade/scratch/jhollowed/CAM/cases_589/project4'
    RUNS = np.array(glob.glob('{}/*/run'.format(topdir)))
    #sss = sys.argv[1]
    #RUNS = np.array(glob.glob('{}/{}/run'.format(topdir, sss)))
    ssw_files = []
    mean_files = []
    vort_files = []
    vort_wind_files = []

    # -------- concat all data --------
    
    if(1):
        for i in range(len(RUNS)):
            run = RUNS[i]
            run_name = run.split('/')[-2]
            run_name = run_name.split('__RESUBMIT')[0]
            dycore = run.split('/')[-2].split('_')[0]
            grid = run.split('/')[-2].split('_')[1]
            mod = run.split('/')[-2].split('_')[2]
            
            ssw_files.append('{}/{}.CAT.SSW.nc'.format(run, run_name))
            mean_files.append('{}/{}.CAT.MEANS.nc'.format(run, run_name))
            vort_wind_files.append('{}/{}.CAT.VORT_WINDS.nc'.format(run, run_name))
            vort_files.append('{}/{}.CAT.VORT.nc'.format(run, run_name))
            
            ctb.concat_run_outputs(run, sel={'lat':75, 'lev':10}, mean=['lon'], 
                              histnum=0, overwrite=False, outFile=ssw_files[i])
            ctb.concat_run_outputs(run, sel={'time':slice('0015-01-01', '0035-01-01')}, mean=['lon', 'time'], 
                              histnum=0, overwrite=False, outFile=mean_files[i])
            ctb.concat_run_outputs(run, sel={'time':slice('0015-01-01', '0035-01-01')}, mean=['time'], 
                              histnum=0, overwrite=False, outFile=vort_wind_files[i])
    
    # --- compute vorticity for all means for aoa slice plotting
    if(1):
        print('\n========== COMPUTING VORTICITY ==========')
        for i in range(len(vort_wind_files)):
            ctb.compute_vorticity(vort_wind_files[i], ncdf_out=vort_files[i], overwrite=False)
    mean_files = np.array(mean_files)
    ssw_files = np.array(ssw_files)
    
   
    # ----- group runs
    mod3_ensRuns = RUNS[np.logical_and(['FV3_C48L72_mod3' in r for r in RUNS], 
                                       ['vt8' not in r for r in RUNS])]
    mod3_vt8_ensRuns = RUNS[['vt8' in r for r in RUNS]]
    mod4_ensRuns = RUNS[['FV3_C48L72_mod4' in r for r in RUNS]]
    SE_runs = RUNS[['SE' in r for r in RUNS]]

    # --- get ensemble member means
    mod3_vt8_means = [glob.glob('{}/*.CAT.MEANS.nc'.format(d))[0] for d in mod3_vt8_ensRuns]
    mod3_means = [glob.glob('{}/*.CAT.MEANS.nc'.format(d))[0] for d in mod3_ensRuns]
    mod4_means = [glob.glob('{}/*.CAT.MEANS.nc'.format(d))[0] for d in mod4_ensRuns]
    
    mod3_vt8_means_vort = [glob.glob('{}/*.CAT.VORT.nc'.format(d))[0] for d in mod3_vt8_ensRuns]
    mod3_means_vort = [glob.glob('{}/*.CAT.VORT.nc'.format(d))[0] for d in mod3_ensRuns]
    mod4_means_vort = [glob.glob('{}/*.CAT.VORT.nc'.format(d))[0] for d in mod4_ensRuns]

    # --- compute ensemble mean, std of ensemble member means
    mod3_ensMean_fname = '{}/FV3_C48L72_mod3_ensMean.nc'.format(topdir)
    mod3_vt8_ensMean_fname = '{}/FV3_C48L72_mod3_vt8_ensMean.nc'.format(topdir)
    mod4_ensMean_fname = '{}/FV3_C48L72_mod4_ensMean.nc'.format(topdir)
    fv3_ensMeans = [mod3_ensMean_fname, mod3_vt8_ensMean_fname, mod4_ensMean_fname]
    fv3_ensMeans_err = ['{}.std'.format(f) for f in fv3_ensMeans]
    
    mod3_ensMean_fname_vort = '{}/FV3_C48L72_mod3_ensMean_vort.nc'.format(topdir)
    mod3_vt8_ensMean_fname_vort = '{}/FV3_C48L72_mod3_vt8_ensMean_vort.nc'.format(topdir)
    mod4_ensMean_fname_vort = '{}/FV3_C48L72_mod4_ensMean_vort.nc'.format(topdir)
    fv3_ensMeans_vort = [mod3_ensMean_fname, mod3_vt8_ensMean_fname, mod4_ensMean_fname]
    fv3_ensMeans_err_vort = ['{}.std'.format(f) for f in fv3_ensMeans]

    if(1):
        print('\n========== COMPUTING ENSEMBLE MEANS ==========')
        ctb.ensemble_mean(mod3_means, std=True, overwrite=False,
                          outFile = mod3_ensMean_fname)
        ctb.ensemble_mean(mod3_vt8_means, std=True, overwrite=False,
                          outFile = mod3_vt8_ensMean_fname)
        ctb.ensemble_mean(mod4_means, std=True, overwrite=False,
                          outFile = mod4_ensMean_fname)
        
        ctb.ensemble_mean(mod3_means_vort, std=True, overwrite=False,
                          outFile = mod3_ensMean_fname_vort)
        ctb.ensemble_mean(mod3_vt8_means_vort, std=True, overwrite=False,
                          outFile = mod3_vt8_ensMean_fname_vort)
        ctb.ensemble_mean(mod4_means_vort, std=True, overwrite=False,
                          outFile = mod4_ensMean_fname_vort)

    # --- gather files desired for aoa slice plotting
    se_files = mean_files[['SE' in f for f in mean_files]]
    aoa_slice_runs = np.hstack([se_files, fv3_ensMeans])
    aoa_slice_errs = np.hstack([[None]*len(se_files), fv3_ensMeans_err])

    # --- AOA slices for 3 SE, 3 FV3 ens means 
    if(0):
        for i in range(len(aoa_slice_runs)):
            aoa_profiles(aoa_slice_runs[i], cvar='U', cvar_label = 'u  [m/s]', 
                         twoPanel=False, sfx='_u')
            aoa_profiles(aoa_slice_runs[i], cvar=xr.open_dataset(vort_files[i])['VORT'], 
                         cvar_label = 'vorticity  [1/s]', twoPanel=False, sfx='_vort')
    
    # --- Gupta fig 9 for 3 SE, 3 FV3 ens means
    if(1):
        aoa_levels(aoa_slice_runs, aoa_slice_errs)


    # --- SSW plots for 3 SE, 3 FV3 ens1
    ssw_plot_runs = np.hstack([se_files, 
                               mean_files[['ens1' in f for f in mean_files]]])
    if(0):
        for i in range(len(aoa_slice_runs)):
            ssw_panels(aoa_slice_runs[i], phissw=75, levssw=10, twoPanel=False, drawSSW=True)
        
        
    # --- SSW tally plots for 3 FV3 ens means
    ssw_count_runs = [mod3_ensRuns, mod3_vt8_ensRuns, mod4_ensRuns, 
                      se_files[0], se_files[1], se_files[2], se_files[3]]
    if(0):
        ssw_counts(aoa_slice_runs, aoa_slice_errs)
