import os
import sys
import pdb
import glob
import Nio, Ngl
import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from cftime import DatetimeNoLeap
import matplotlib.ticker as ticker

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/..')
import climate_toolbox as ctb



# ===============================================================================================
# ===============================================================================================


def isolate_barnes_period_data(my_vars=None):
    '''
    Gather data for all CMIP6 variables in the period of 1980-2004, write out to 
    file per-variable

    Error messages containing 'Killed' means the progam has run out of memory; try
    again on an interactive slurm job

    Parameters
    ----------
    my_vars : str list
        list of variable names to isolate, if only a subsection is desired. 
        Default is None, in which case, all variables will be processed
     '''

    # don't consider any anomaly or climatology files that may be present here already
    cmip_files = CMIP_FILES[~np.array(['nomaly' in f for f in CMIP_FILES])]
    cmip_files = cmip_files[~np.array(['climatology' in f for f in cmip_files])]
    
    # gather metadata from file names
    realizations = np.array([int(f.split('_r')[-1].split('i')[0]) for f in cmip_files])
    variables = np.array([f.split('/')[-1].split('_')[0] for f in cmip_files])

    # only keep user-desired vars
    if(my_vars is not None):
        my_vars_mask = [var in my_vars for var in variables]
        cmip_files = cmip_files[my_vars_mask]
        realizations = realizations[my_vars_mask]
        variables = variables[my_vars_mask]
    
    var_names = np.unique(variables)
    
    # construct mask on files which contain any data in the pinatubo period
    # (data corresponding to the precise period will be extracted from these files later)
    times = np.array([f.split('/')[-1].split('_')[-1].split('.nc')[0].split('-') for f in cmip_files])
    for time in times:
        time[0] = int(time[0][:-2])
        time[1] = int(time[1][:-2])
    time_mask = [1 in t or t.tolist() == [0, 2] for t in np.searchsorted(PINATUBO_PERIOD, times)]

    print('Found {} files for vars {} on realizations {} over time periods {}'.format(
                len(cmip_files), var_names, np.unique(realizations), times[time_mask]))

    # loop over realiozations, variables
    for i in range(len(np.unique(realizations))):
        rlz_mask = realizations == i+1
        print('\n ---------- working on realization {}'.format(i+1))
        
        for j in range(len(var_names)):
            print('gathering data for variable {}'.format(var_names[j]))
           
            # identify files belogning to this variable, this realization, relevant time period
            var_mask = variables == var_names[j]
            file_mask = rlz_mask & var_mask & time_mask
            
            # this sort is important;
            # time slice bounds need to be unique values present in the dataset unless 
            # the sliced dimension is sorted monotonically. Sorting the data after reading,
            # which will require lots of memory, or we can just mak sure that we read them 
            # in in temporally ascending order. Here we assume that the filenames are 
            # identical except for the date range at the end, in which case this is all we 
            # need to do
            files = sorted(cmip_files[file_mask])
            print('{} files after masking'.format(len(files)))
 
            # concatenate contents of relevant time intervals to dataset...
            dsets = np.empty(len(files), dtype=xr.core.dataset.Dataset)
            for k in range(len(files)):
                dsets[k] = xr.open_dataset(files[k])
            var_alltime = xr.concat(dsets, dim='time')
           
            # isolate Pinatubo period
            var_out = var_alltime.sel(time=slice(str(PINATUBO_PERIOD[0]), str(PINATUBO_PERIOD[1])))
            
            # write out NetCDF file
            start = var_out['time'].values[0]
            end = var_out['time'].values[-1]
            startstr = '{}{}'.format(start.year, f'{start.month:02}')
            endstr = '{}{}'.format(end.year, f'{end.month:02}')
            dest = '{}/{}_mon_E3SM-1-0_{}_r{}_{}-{}.nc'.format(PINATUBO, var_names[j], MODEL_TYPE,
                                                                i+1, startstr, endstr) 
            print('writing out to {}'.format(dest))
            var_out.to_netcdf(dest, format='NETCDF4')

        
# ===============================================================================================


def compute_climatology(my_vars=None):
    '''
    Computes the monthly climatology of the CMIP6 data for the pinatubo period, and writes out to
    a NetCDF file
    
    Parameters
    ----------
    my_vars : str list
        list of variable names to isolate, if only a subsection is desired. 
        Default is None, in which case, all variables will be processed
    '''
    
    # don't consider any anomaly or climatology files that may be present here already
    pinatubo_files = PINATUBO_FILES[~np.array(['nomaly' in f for f in PINATUBO_FILES])]
    pinatubo_files = pinatubo_files[~np.array(['climatology' in f for f in pinatubo_files])]
    
    # gather metadata from file names
    realizations = np.array([int(f.split('/')[-1].split('_r')[-1].split('_')[0]) for f in \
                                                                                 pinatubo_files])
    variables = np.array([f.split('/')[-1].split('_')[0] for f in pinatubo_files])
    
    # only keep user-desired vars
    if(my_vars is not None):
        my_vars_mask = [var in my_vars for var in variables]
        pinatubo_files = pinatubo_files[my_vars_mask]
        realizations = realizations[my_vars_mask]
        variables = variables[my_vars_mask]
    
    var_names = np.unique(variables)
    
    for i in range(len(np.unique(realizations))):
        i += 1
        rlz_mask = realizations == i
        print('\n---------- working on realization {}'.format(i))
        
        for j in range(len(var_names)):
            print('computing monthly climatology for variable {}'.format(var_names[j]))
            
            # get this variable for this realization, open file
            var_mask = variables == var_names[j]
            file_mask = rlz_mask & var_mask
            files = pinatubo_files[file_mask]
            if(len(files) != 1):
                raise RuntimeError('More than one file found for {} in realization r{}:\n{}'.format(
                                    var_names[j], i, files))
            dset = xr.open_dataset(files[0])
           
            # compute climatology, write out
            out_file = '{}{}_climatology_{}'.format(files[0].split(var_names[j])[0],
                                               var_names[j], files[0].split('mon_')[-1])
            ctb.compute_climatology(dset, out_file)


# ===============================================================================================
            

def compute_anomalies(my_vars=None):
    '''
    Compute the monthly anomalies of variable x with respect to the mean climatology.
    Write out result per each variable, per each ensemble member, and the mean of the
    three members per each variable
    
    Parameters
    ----------
    my_vars : str list
        list of variable names to isolate, if only a subsection is desired. 
        Default is None, in which case, all variables will be processed
    '''
    
    # don't consider any anomaly files that may be present here already
    pinatubo_files = PINATUBO_FILES[~np.array(['nomaly' in f for f in PINATUBO_FILES])]
    
    # gather metadata from file names
    realizations = np.array([int(f.split('/')[-1].split('_r')[-1].split('_')[0]) for f in 
                             pinatubo_files])
    variables = np.array([f.split('/')[-1].split('_')[0] for f in pinatubo_files])
    
    # only keep user-desired vars
    if(my_vars is not None):
        my_vars_mask = [var in my_vars for var in variables]
        pinatubo_files = pinatubo_files[my_vars_mask]
        realizations = realizations[my_vars_mask]
        variables = variables[my_vars_mask]
    
    var_names = np.unique(variables)
    total_realz = len(np.unique(realizations))

    for j in range(len(var_names)):
        print('\n---------- computing monthly anomalies for variable {}'.format(var_names[j]))
        for i in range(total_realz):
            i += 1
            rlz_mask = realizations == i
            print('working on realization {}'.format(i))

            # get this variable for this realization, and climatology
            var_mask = variables == var_names[j]
            file_mask = rlz_mask & var_mask
            files = pinatubo_files[file_mask]
            if(len(files) != 2):
                raise RuntimeError('Expected 2 files, climatology and monthly avg time series,'\
                                   'for {} in realization r{}, but found {}:\n{}'.format(
                                    var_names[j], i, len(files), files))
            
            # read in
            mask = np.array(['climatology' in f for f in files])
            fclimatology = files[mask][0]
            fmonthly = files[~mask][0]
            dest = '{}{}_anomaly_{}'.format(fmonthly.split(var_names[j])[0],
                                            var_names[j], fmonthly.split('mon_')[-1])

            # compute anomalies, update mean
            # the groupby() call for the monthly object is very important! subtracting climatology
            # will have unexpected results otherwise, and cause the data volume to increase 12-fold
            # (skip this block if the current anomaly already exists as output;  good for runs
            # that have crashed)
            if(not os.path.isfile(dest)): 
                climatology = xr.open_dataset(fclimatology)
                monthly = xr.open_dataset(fmonthly)
               
                print('computing anomaly')
                anomaly = monthly.groupby('time.month') - climatology
                
                # write out anomaly
                monthly.close()
                climatology.close()
                anomaly.to_netcdf(dest, format='NETCDF4')
            else:
                print('anomaly for var {},{} already exists; reading in'.format(
                                                             var_names[j], i))
                anomaly = xr.open_dataset(dest)

            # update ensemble mean
            print('updating ensemble mean')
            if(i == 1):
                ensemble_mean = anomaly
            else:
                ensemble_mean += anomaly
                if(i == total_realz):
                    print('writing out ensemble mean')
                    ensemble_mean /= total_realz
                    dest = '{}{}_meanAnomaly_{}'.format(fmonthly.split(var_names[j])[0], 
                                                        var_names[j], fmonthly.split('_')[-1])
                    ensemble_mean.to_netcdf(dest, format='NETCDF4')
                    ensemble_mean.close()


# ===============================================================================================

def barnes_fig1(horizontal = True, showfig=False):
    '''
    Renders a fig comparable to Barnes+ Figure 1 for the following outputs sets
    - individual ensembles
    - ensemble mean

    Parameters
    ----------
    horizontal : bool
        Whether or not to render the plot rotated at 90 degrees with respect to 
        the form in Barnes+ (a horinzontal arrangement,read from right-left, rather 
        than top-bottom)
    showfig : bool
        Whether or not to display the figure, rather than saving to file. Defaults 
        to False
    '''
  
    # get files corresponding to ensemble member anomalies, and the ensemble mean anomaly...
    # the sorted() call ensures that the mean is plotted last, which is important to be able
    # to plot ensemble sign agreement as stippling... requires that ensemble anomaly file 
    # naming starts like 
    # 'variable_anomaly...', 
    # while the mean starts like 
    # 'variable_meanAnomaly...'
    u_files = sorted(glob.glob('{}/ua_*[aA]nomaly*'.format(PINATUBO)))
    t_files = sorted(glob.glob('{}/ta_*[aA]nomaly*'.format(PINATUBO)))
    assert len(u_files) == len(t_files), 'Found {} uf files, but only {} ta files'.format(
                                                                len(uf_files), len(ta_files))

    # define the period used for the Barnes+ Figure 1 plots
    tslice = slice('1991-07', '1992-02')

    # loop over members, mean
    for j in range(len(u_files)):
        uf = u_files[j]
        tf = t_files[j]

        # check that the u and t files are from the same simulation
        u_is_ensMean = 'meanAnomaly' in uf
        t_is_ensMean = 'meanAnomaly' in tf
        u_ensemble_num = uf.split('/')[-1].split('{}_r'.format(MODEL_TYPE))[-1].split('_')[0]
        t_ensemble_num = tf.split('/')[-1].split('{}_r'.format(MODEL_TYPE))[-1].split('_')[0]
        
        u_identifier = ['ensemble{}'.format(u_ensemble_num), 'ensembleMean'][u_is_ensMean]
        t_identifier = ['ensemble{}'.format(t_ensemble_num), 'ensembleMean'][t_is_ensMean]
        assert u_identifier == t_identifier, 'u identifier \'{}\' and t identifier \'{}\ '\
                                             'do not match; check files'.format(
                                                                         u_dentifier, t_identifier)
        identifier = u_identifier
        print('\n ===== Rendering Figure 1 for {}_{}...'.format(MODEL_TYPE, identifier))
        print('found {} u,t files'.format(len(u_files)))
        
        # get time slice, take zonal average
        print('opening datasets')
        u = (xr.open_dataset(uf)).sel(time=tslice)
        t = (xr.open_dataset(tf)).sel(time=tslice)
        print('taking zonal average')
        uu = u.mean('lon')
        tt = t.mean('lon')

        # set max wind anomaly for colro scale
        #max_level = int(np.max(uu['ua'])*0.9)
        max_level=12
        levels = np.linspace(-max_level, max_level, 18, dtype=int)

        # take union of sign of response with other ensemble members
        if(not u_is_ensMean):
            if(j == 0):
                ens_sign_union = (np.sign(uu['ua'])+1)/2
            else:
                old = ens_sign_union
                ens_sign_union = np.logical_and(ens_sign_union, (np.sign(uu['ua'])+1)/2)
                                                           

        # plot vertical cross section for each month 
        if(horizontal):
            order = [0, 1, 2, 3, 4, 5, 6, 7]
            give_ylabel = [0, 4]
            give_xlabel = [4, 5, 6, 7]
            nrow, ncol = 2, 4
            cbar_loc = 'right'
            figsize = (18, 9)
            labelfs = 13
            titlefs = 14
            tickfs = 11
        else:
            order = [0, 2, 4, 6, 1, 3, 5, 7]
            give_ylabel = [0, 1, 2, 3]
            give_xlabel = [3, 7]
            nrow, ncol = 4, 2
            cbar_loc = 'bottom'
            figsize=(9, 15)
            labelfs = 13
            titlefs = 14
            tickfs = 11

        fig, axes = plt.subplots(nrows=nrow, ncols=ncol, figsize=figsize, constrained_layout=True)
        for i in range(uu['ua'].shape[0]):

            print('plotting panel {}'.format(order[i]))
            uuu = uu['ua'][i]
            ttt = tt['ta'][i]
            lat = uu['lat']
            lev = uu['plev']
            time = uu['time'].values[i]
            ax = np.ravel(axes)[order[i]]

            # plot zonal velocity anomaly
            uc = ax.contourf(lat, lev/100, uuu, levels=levels, cmap='RdBu_r', extend='both')

            # ----- TMP: checking node magntiude positions -----
            #ax.plot([-75, -75], [min(lev/100), max(lev/100)], '--k')
            #ax.plot([-45, -45], [min(lev/100), max(lev/100)], '--k')
            #----------------

            # plot temperature anomaly on contours every 0.5 K, including zero
            t_levels = np.arange(np.floor(np.min(ttt)*2)/2, np.floor(np.max(ttt) * 2)/2, 0.5)
            zero = t_levels.tolist().index(0)
            ax.contour(lat, lev/100, ttt, levels = np.hstack([t_levels[0:zero], 
                                                              t_levels[zero+1:]]), 
                                           colors='k', alpha=0.55, linewidths=0.5)
            ax.contour(lat, lev/100, ttt, levels = [t_levels[zero]], colors='k', 
                                                        alpha=1, linewidths=0.75)
            #ax.clabel(cs, cs.levels, inline=True, fontsize=tickfs, fmt=lambda x:'{:.1f}'.format(x)) 
            
            
            # plot ensemble member agreement on sign of the anomaly as stippling, 
            # downsampling the zonal coordinate by a fraction lat_downsample
            if(u_is_ensMean):
               
                lat_downsample = 0.25
                lat_step = int(1/lat_downsample)
                lat = lat[1::lat_step]
                stippling = ens_sign_union[i][:, 1::lat_step]
                
                glat, glev = np.meshgrid(lat, lev/100)
                alphas = np.ravel(stippling).astype(float)
                colors = np.zeros((len(np.ravel(glat)), 4))
                colors[:,3] = alphas
                ax.scatter(np.ravel(glat), np.ravel(glev), color=colors, s=0.25)
          

            ax.set_ylim([10, 1000])
            ax.set_yscale('log')
            ax.invert_yaxis()
            ax.set_ylabel('Pressure (hPa)', fontsize=labelfs)
            ax.set_xlabel('Latitude (°N)', fontsize=labelfs)
            ax.set_title('{}-{}'.format(time.month, time.year), fontsize=titlefs)
            ax.tick_params(axis='both', which='major', labelsize=tickfs)
            
            if(i not in give_ylabel):
                ax.set_ylabel('')
            if(i not in give_xlabel):
                ax.set_xlabel('')
            
        suptitle = ['Ensemble {}'.format(identifier.split('ensemble')[-1]), 
                    'Ensemble Mean'][u_is_ensMean]
        fig.suptitle(suptitle, fontsize=11)
        
        # make the single for all subplots; since range is the same for all
        # just use the last created contour plot as the mappable
        cbar = fig.colorbar(uc, ax=axes.ravel().tolist(), location=cbar_loc, aspect=50)
        cbar.set_ticks(levels)
        cbar.ax.tick_params(labelsize=tickfs)
        if(horizontal):
            cbar.ax.set_ylabel('u  [m/s]', fontsize=labelfs)
        else:
            cbar.ax.set_xlabel('u  [m/s]', fontsize=labelfs)


        print('saving...')
        if(not showfig):
            plt.savefig('{}/CMIP6_BarnesFig1_{}_{}.png'.format(
                         FIGDIR, MODEL_TYPE, identifier), dpi=300)
        else:
            plt.show()


# ===============================================================================================


def verify_qbo():
    '''
    Verifys the the strong zonal wind anomaly signal present in the Barnes Fig.1 plots (function above)
    are indeed QBO signatures, by plotting 
    - several lat-lev cross sections across a time interval of a decade
    - time series of the equatoiral stratospheric zonal-mean zonal wind over the full domain
    for the following runs:
    - one AMIP ensemble member
    - one historical ensemble member
    '''
  
    # get files corresponding to ensemble member anomalies
    # the sorted() call is meant to target the first ensemble member; doesn't really matter
    uf = sorted(glob.glob('{}/ua_*anomaly*'.format(PINATUBO)))
    u_ensemble_num = [f.split('/')[-1].split('{}_r'.format(MODEL_TYPE))[-1].split('_')[0] for f in uf]
    identifier = ['ensemble{}'.format(n) for n in u_ensemble_num]

    # define the period used for vertical slices, and lat range used for time series
    tslice = slice('1985-01', '1995-01')
    latslice = slice(-5, 5)

    # ========== plot vertical slices for ens 1==========
    for j in range(len(uf)):
        if j > 0: continue
        print('\n ===== Rendering vertical slices for {}_{}...'.format(MODEL_TYPE, identifier[j]))
        
        # get time slice, take zonal average
        print('opening dataset')
        u = (xr.open_dataset(uf[j])).sel(time=tslice)
        print('taking zonal average')
        uu = u.mean('lon')

        # set max wind anomaly for color scale
        max_level = 30
        levels = np.linspace(-max_level, max_level, 18, dtype=int)

        order = [0, 1, 2, 3, 4, 5, 6, 7]
        give_ylabel = [0]
        nrow, ncol = 1, 8
        nplots = nrow*ncol
        cbar_loc = 'right'
        figsize = (18, 4)
        labelfs = 12
        titlefs = 14
        tickfs = 10

        fig, axes = plt.subplots(nrows=nrow, ncols=ncol, figsize=figsize, constrained_layout=True)
        for i in range(nplots):

            print('plotting panel {}'.format(order[i]))
            ti = int(len(u['time'])/(14))
            time = uu['time'].values[::ti][0:9][i]
            uuu = uu['ua'].values[::ti][0:9][i]
            lat = uu['lat']
            lev = uu['plev']
            ax = np.ravel(axes)[order[i]]

            # plot zonal velocity anomaly
            uc = ax.contourf(lat, lev/100, uuu, levels=levels, cmap='RdBu_r', extend='both') 
            
            ax.set_ylim([10, 100])
            ax.set_yscale('log')
            ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
            ax.yaxis.set_minor_formatter(ticker.ScalarFormatter())
            plt.ticklabel_format(axis='y', style='plain')
            ax.invert_yaxis()
            ax.set_ylabel('Pressure (hPa)', fontsize=labelfs)
            ax.set_xlabel('Latitude (°N)', fontsize=labelfs)
            ax.set_title('{}-{}'.format(time.month, time.year), fontsize=titlefs)
            ax.tick_params(axis='both', which='major', labelsize=tickfs) 
            ax.tick_params(axis='both', which='minor', labelsize=tickfs) 
            if(i not in give_ylabel):
                ax.set_ylabel('')
                ax.set_yticklabels([])
                ax.tick_params(axis='y', which='major', labelsize=0) 
                ax.tick_params(axis='y', which='minor', labelsize=0) 
        suptitle = 'Ensemble {}'.format(identifier[j].split('ensemble')[-1])
        fig.suptitle(suptitle, fontsize=11)
        
        # make single colorbar for all subplots; since range is the same for all, use the last 
        # created contour plot as the mappable
        cbar = fig.colorbar(uc, ax=axes.ravel().tolist(), location=cbar_loc, aspect=50)
        cbar.set_ticks(levels)
        cbar.ax.tick_params(labelsize=tickfs)
        cbar.ax.set_ylabel('u  [m/s]', fontsize=labelfs)
        print('saving...')
        plt.savefig('{}/CMIP6_QBOPanels_{}_{}.png'.format(FIGDIR, MODEL_TYPE, identifier[j]), dpi=300)
    
    
    # ========== plot time series ==========
    for j in range(len(uf)):
        print('\n ===== Rendering time series for {}_{}...'.format(MODEL_TYPE, identifier[j]))
        
        # get time slice, take zonal average
        # kinda slow so lets save the result for replotting
        savedest = '{}/uu_qboMean_{}_{}.nc'.format(FIGDATA, MODEL_TYPE, identifier[j])
        if(not os.path.isfile(savedest)): 
            print('opening dataset')
            u = (xr.open_dataset(uf[j])).sel(time=tslice)
            u = u.sel(lat=latslice)
            print('taking meridional average on +-5 deg.')
            uu = u.mean('lat')
            print('taking zonal average')
            uu = uu.mean('lon')
            uu.to_netcdf(savedest, format='NETCDF4')
        else:
            print('reading QBO means from file')
            uu = xr.open_dataset(savedest)

        # set max wind anomaly for color scale
        max_level = 40
        levels = np.linspace(-max_level, max_level, 18, dtype=int)
        
        figsize = (8, 3)
        labelfs = 10
        titlefs = 11
        tickfs = 9
        cmap = 'Spectral_r'
        #cmap = 'RdBu_r'

        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111)
        uuu = uu['ua']
        lev = uu['plev']
        time = uu['time']
        timei = np.arange(0, len(time))
        year = np.array([t.year for t in time.values])
        month = np.array([t.month for t in time.values])
        erupt_time = timei[np.logical_and(year==1991, month==7).tolist().index(1)]

        # plot zonal velocity anomaly
        uc = ax.contourf(timei, lev/100, uuu.T, levels=levels, cmap=cmap, extend='both')
        ax.plot([erupt_time]*2, [10,100], '--k', label='Pinatubo eruption', lw=1)
        
        ax.legend(fontsize=labelfs)
        ax.set_ylim([10, 100])
        ax.set_yscale('log')
        ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
        ax.yaxis.set_minor_formatter(ticker.ScalarFormatter())
        plt.ticklabel_format(axis='y', style='plain')
        ax.invert_yaxis()
        ax.set_ylabel('Pressure (hPa)', fontsize=labelfs)
        ax.set_xlabel('Time', fontsize=labelfs)
        
        ax.set_xticks(np.arange(0, len(time))[::12])
        tick_labels = [str(t.year) for t in time.values[::12]]
        ax.set_xticklabels(tick_labels)
        plt.xticks(rotation=90)
        ax.tick_params(axis='both', which='major', labelsize=tickfs) 
        ax.tick_params(axis='both', which='minor', labelsize=tickfs) 

        ax.set_title('Mean zonal wind anomaly', fontsize=titlefs)
        cbar = fig.colorbar(uc, ax=ax, aspect=20)
        cbar.set_ticks(uc.levels[::2])
        cbar.ax.tick_params(labelsize=tickfs)
        cbar.ax.set_ylabel('u  [m/s]', fontsize=labelfs)
        plt.tight_layout()
        print('saving...')
        plt.savefig('{}/CMIP6_QBO_{}_{}.png'.format(FIGDIR, MODEL_TYPE, identifier[j]), dpi=300)


# ===============================================================================================


def comp_poleward_node_mag(u, hem='S'):
    '''
    Computes the "poleward node magnitude" of the zonal wind anomaly, as defined in Barnes+

    Parameters
    ----------
    u : DataArray
        zonal wind anomaly
    hem : str
        Which hemisphere to search for the node. Either 'S' or 'N'. Default to 'S'

    Returns
    -------
    up : DataArray
        The poleward node magnitude, in m/s, as a function of pressure and time
    '''
   
    # define bounds for the poleward node
    s = [-1,1][['S', 'N'].index(hem)]
    bounds = sorted([s*45, s*75])
    latslice = slice(bounds[0], bounds[1])
    u = u.sel(lat=latslice)

    # take zonal mean, compute the node magnitude
    u = u.mean('lon')
    nodemag = u.max(dim='lat')
    return nodemag
    
   
# ===============================================================================================


def qbo_index_corr(qbo_lev=30, do_time_avg=True, qbo_avg='year', hem='S', showfig=False):
    '''
    Render a plot comparing the pinatubo response to the QBO index

    Parameters
    ----------
    qbo_lev : float
        Pressure level at which to take the zonal-mean zonal wind, in hPa. The standard
        choices for this index at 30 hPa (the default), and 50 hPa. Defaults to 30hPa. If
        this pressure level not present in the data, it will be linearly interpolated
    do_time_avg : bool
        Whether or not to perform time averaging. Defaults to True, in which case one point
        will be rendered on the plot per ensemble member, each index being averaged over the
        relevant time interval. If False, the data will be plotting as a time series, with
        each ensemble member contributing one "trajectory" in the phase plane
    showfig : bool
        Whether or not to display the figure, rather than saving to file. Defaults 
        to False
    hem : str
        Which hemisphere to use in the calculation of the poleward node magnitude. Must be
        either 'S' or 'N'. Defaults to 'S'.
    qbo_avg : str
        Defines what time range to use for the QBO averaging. If 'year', the average is taken
        over the year following the eruption. If 'eruption', the average is taken over the 2
        months following the eruption
    '''
    
    print('\n ---------- Doing QBO dependency analysis ----------')
    print('The QBO index will be taken at {} hPa\nPlot will be rendered with do_time_avg={}\
           \nQBO mean will be taken over period:{}\nHemisphere for node calculation:{}H'.format(
           qbo_lev, do_time_avg, qbo_avg, hem))
    
    # get files corresponding to ensemble member anomalies
    uf = glob.glob('{}/ua_*anomaly*'.format(PINATUBO))
    u_ensemble_num = [f.split('/')[-1].split('_')[-2].strip('r') for f in uf]
    identifier = ['ensemble{}'.format(n) for n in u_ensemble_num]
    
    # time ranges for QBO and poleward node mag averaging
    # per Barnes, this should cover the year following the eruption for the node mag
    # if do_time_avg=True, these slices need not coincide, and can use whatever 
    # definitions should give the best measure...
    # if do_time_avg=False, we are plotting trajectories, so set the time interval 
    # for each metric to one matching the Barnes+ Fig 1 panels
    if(do_time_avg):
        if(qbo_avg=='eruption'):
            qbo_tslice = slice('1991-06', '1992-08')
        if(qbo_avg=='year'):
            qbo_tslice = slice('1991-06', '1992-05')
        node_tslice = slice('1991-06', '1992-05')
    else:
        qbo_tslice = slice('1991-07', '1992-02')
        node_tslice = slice('1991-07', '1992-02')


    # create plot...
    fig = plt.figure()
    ax = fig.add_subplot(111)
    colors = ['r', 'b']
    labels = np.array(['historical', 'AMIP'])
    # make dummy plots for the legend...
    ax.plot([0],[4],'d',color=colors[0],label=labels[0])
    ax.plot([0],[4],'d',color=colors[1],label=labels[1])
    ax.plot([0],[4],'dw', ms=10)
    

    # loop through ensemble members...
    for i in range(len(uf)):
        print('working on {}'.format(identifier[i]))
        
        # read zonal wind
        u = xr.open_dataset(uf[i])
        
        # get qbo index, average over time window
        u_qbo = u.sel(time = qbo_tslice)
        qbo = ctb.QBO_index(u_qbo, qbo_lev)
        if(do_time_avg):
            qbo_tmean = qbo.mean('time')

        # get poleward node magnitude, vertically average, time average
        u_node = u.sel(time = node_tslice)
        node_mag = comp_poleward_node_mag(u_node, hem=hem)
        node_mag = node_mag.mean('plev')
        if(do_time_avg):
            node_mag_tmean = (node_mag.mean('time'))['ua']
        node_mag = node_mag['ua']

        # ---------- plot ----------
        color = colors[int('AMIP' in uf[i])] 
        if(do_time_avg):
            ax.plot(qbo_tmean, node_mag_tmean, 'dk', color=color)
        else:
            alphas = np.linspace(0.1, 0.8, len(qbo)-1)
            for j in range(len(qbo)-1):
                ax.plot([qbo[j], qbo[j+1]], [node_mag[j], node_mag[j+1]], 
                         '-', color=color, alpha=alphas[j])
            ax.plot(qbo[-1], node_mag[-1], '.', ms=10, color=color)

        #format
        ax.legend(fontsize=12, loc='upper left')
        if(qbo_avg=='eruption'):
            ax.set_xlabel('3-month QBO {}hPa Index  [m/s]'.format(qbo_lev), fontsize=13)
        if(qbo_avg=='year'):
            ax.set_xlabel('12-month QBO {}hPa Index  [m/s]'.format(qbo_lev), fontsize=13)
        ax.set_ylabel('{}H Poleward Node Magnitude  [m/s]'.format(hem), fontsize=13)

    # add zero line
    ylim=ax.get_ylim()
    ax.plot([0,0], ylim, '--k', lw=0.8, alpha=0.8)
    ax.set_ylim(ylim)

    # save figure
    if(not showfig):
        ss = np.array(['trajectory', 'mean'])[int(do_time_avg)]
        tt = np.array(['3', '12'])[int(qbo_avg=='year')] 
        out = '{}/CMIP6_{}moQBOcorr_{}hPa_{}_{}H_{}.png'.format(
               FIGDIR, tt, qbo_lev, MODEL_TYPE, hem, ss)
        print('saving figure to {}'.format(out.split('/')[-1]))
        plt.savefig(out, dpi=300)
    else:
        plt.show()
        
        

# ===============================================================================================


def nino_index_corr(do_time_avg=True, hem='S', showfig=False):
    '''
    Render a plot comparing the pinatubo response to the Nino3.4 index

    Parameters
    ----------
    do_time_avg : bool
        Whether or not to perform time averaging. Defaults to True, in which case one point
        will be rendered on the plot per ensemble member, each index being averaged over the
        relevant time interval. If False, the data will be plotting as a time series, with
        each ensemble member contributing one "trajectory" in the phase plane
    showfig : bool
        Whether or not to display the figure, rather than saving to file. Defaults 
        to False
    hem : str
        Which hemisphere to use in the calculation of the poleward node magnitude. Must be
        either 'S' or 'N'. Defaults to 'S'.
    '''
    
    print('\n ---------- Doing ENSO dependency analysis ----------')
    print('Plot will be rendered with do_time_avg={}, nodes computed in the {}H'.format(
           do_time_avg, hem))
    
    # get files corresponding to ensemble member anomalies
    uf = sorted(glob.glob('{}/ua_*anomaly*'.format(PINATUBO)))
    u_ensemble_num = [f.split('/')[-1].split('_')[-2].strip('r') for f in uf]
    identifier = ['ensemble{}'.format(n) for n in u_ensemble_num]
    
    sstf = sorted(glob.glob('{}/tos_*anomaly*'.format(PINATUBO)))
    sst_ensemble_num = [f.split('/')[-1].split('_')[-2].strip('r') for f in sstf]

    if(u_ensemble_num != sst_ensemble_num):
        raise RuntimeError('u and sst files don\'t seem to correspond... debug this... \
                            may have forgot to set MODEL_TYPE=\'historical\'')
    
    # time ranges for Nino and poleward node mag averaging
    # per Barnes, this should cover the year following the eruption for the node mag
    # if do_time_avg=True, these slices need not coincide, and can use whatever 
    # definitions should give the best measure...
    # if do_time_avg=False, we are plotting trajectories, so set the time interval 
    # for each metric to one matching the Barnes+ Fig 1 panels
    if(do_time_avg):
        nino_tslice = slice('1991-06', '1991-08')
        node_tslice = slice('1991-06', '1992-05')
    else:
        nino_tslice = slice('1991-07', '1992-02')
        node_tslice = slice('1991-07', '1992-02')


    # create plot...
    fig = plt.figure()
    ax = fig.add_subplot(111) 
    # make dummy plot for the legend...
    ax.plot([0],[4],'d',color='r',label='historical')
    ax.plot([0],[4],'dw', ms=10)

    # loop through ensemble members...
    for i in range(len(uf)):
        print('working on {}'.format(identifier[i]))
        
        # read zonal wind
        u = xr.open_dataset(uf[i])
        sst = xr.open_dataset(sstf[i])
        
        # get qbo index, average over time window
        sst = sst.sel(time = nino_tslice)
        nino = ctb.Nino34_index(sst)
        if(do_time_avg):
            nino_tmean = nino.mean('time')

        # get poleward node magnitude, vertically average, time average
        u_node = u.sel(time = node_tslice)
        node_mag = comp_poleward_node_mag(u_node, hem=hem)
        node_mag = node_mag.mean('plev')
        if(do_time_avg):
            node_mag_tmean = (node_mag.mean('time'))['ua']
        node_mag = node_mag['ua']

        # ---------- plot ----------
        if(do_time_avg):
            ax.plot(nino_tmean, node_mag_tmean, 'dk', color='r')
        else:
            alphas = np.linspace(0.1, 0.8, len(nino)-1)
            for j in range(len(nino)-1):
                ax.plot([nino[j], nino[j+1]], [node_mag[j], node_mag[j+1]], 
                         '-', color='r', alpha=alphas[j])
            ax.plot(nino[-1], node_mag[-1], '.', ms=10, color='r')

        #format
        ax.legend(fontsize=12, loc='upper left')
        ax.set_xlabel('3-month Nino 3.4 Index  [K]', fontsize=13)
        ax.set_ylabel('{}H Poleward Node Magnitude  [m/s]'.format(hem), fontsize=13)

    # add zero line
    ylim=ax.get_ylim()
    ax.plot([0,0], ylim, '--k', lw=0.8, alpha=0.8)
    ax.set_ylim(ylim)

    # save figure
    if(not showfig):
        ss = np.array(['trajectory', 'mean'])[int(do_time_avg)]
        out = '{}/CMIP6_NinoCorr_{}_{}H_{}.png'.format(FIGDIR, MODEL_TYPE, hem, ss)
        print('saving figure to {}'.format(out.split('/')[-1]))
        plt.savefig(out, dpi=300)
    else:
        plt.show()
                    
             
# ===============================================================================================
# ===============================================================================================


if __name__ == '__main__':
    
    MODEL_TYPE = sys.argv[1] # either 'AMIP', 'historical', or 'all'
    USAGE = sys.argv[2] # either 'process' or 'figs'
    if(MODEL_TYPE != 'all'):
        CMIP6 = '/nfs/turbo/cjablono2/hollowed/CMIP6/E3SM_{}/pinatubo_period'.format(MODEL_TYPE)
        PINATUBO = '/nfs/turbo/cjablono2/hollowed/CMIP6/E3SM_{}/pinatubo_period_exact'.format(
                                                                                       MODEL_TYPE)
    else:
        CMIP6 = '/nfs/turbo/cjablono2/hollowed/CMIP6/E3SM_*/pinatubo_period'
        PINATUBO = '/nfs/turbo/cjablono2/hollowed/CMIP6/E3SM_*/pinatubo_period_exact'

    FIGDIR = '/home/hollowed/repos/climate_analysis/CLDERA/figs'
    FIGDATA = '/home/hollowed/repos/climate_analysis/CLDERA/figdata'
    CMIP_FILES = np.array(glob.glob('{}/*.nc'.format(CMIP6)))
    PINATUBO_FILES = np.array(glob.glob('{}/*.nc'.format(PINATUBO)))
    PINATUBO_PERIOD = [1980, 2004]

    if USAGE == 'process':
        #if(len(glob.glob('{}/*'.format(PINATUBO))) == 0):
            isolate_barnes_period_data(my_vars=['tos'])
        #if(len(glob.glob('{}/*climatology*'.format(PINATUBO))) == 0):
            compute_climatology(my_vars=['tos'])
        #if(len(glob.glob('{}/*anomaly*'.format(PINATUBO))) == 0):
            compute_anomalies(my_vars=['tos'])
    elif USAGE == 'figs':
        showfig=False

        # should be run with MODEL_TYPE='all'
        #qbo_index_corr(do_time_avg=True, showfig=showfig, qbo_avg='eruption', hem='S')
        #qbo_index_corr(do_time_avg=True, showfig=showfig, qbo_avg='year', hem='S')
        #qbo_index_corr(do_time_avg=False, showfig=showfig, hem='S')
        #qbo_index_corr(do_time_avg=True, showfig=showfig, qbo_avg='eruption', hem='N')
        #qbo_index_corr(do_time_avg=True, showfig=showfig, qbo_avg='year', hem='N')
        #qbo_index_corr(do_time_avg=False, showfig=showfig, hem='N')
        
        # should be run with MODEL_TYPE='historical'
        nino_index_corr(do_time_avg=True, showfig=showfig, hem='S')
        nino_index_corr(do_time_avg=False, showfig=showfig, hem='S')
        nino_index_corr(do_time_avg=True, showfig=showfig, hem='N')
        nino_index_corr(do_time_avg=False, showfig=showfig, hem='N')
        
        # should be run with choice of MODEL_TYPE
        #barnes_fig1(showfig=showfig)
        #verify_qbo()
       
        
