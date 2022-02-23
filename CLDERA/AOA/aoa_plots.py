import numpy as np
import matplotlib.pyplot as plt
import glob
import pdb
import xarray as xr
from climate_artist import horizontal_slice as plthor
from climate_artist import vertical_slice as pltvert

runs = glob.glob('/glade/scratch/jhollowed/CAM/cases/aoa_runs/project3/*aoa/')

# =============================================================================


def ssw_panels():

    #for kk in [50, 55, 60, 65, 70, 75, 80, 85]:
    for kk in [75]:
        
        for i in range(len(runs)):
            
            fig = plt.figure(figsize=(5,7), constrained_layout=True)
            spec = fig.add_gridspec(4,1)
            ax1 = fig.add_subplot(spec[0,0])
            ax2 = fig.add_subplot(spec[1:3,0])
            ax3 = fig.add_subplot(spec[3,0])

            print('\n\n---------- working on run {} ----------'.format(runs[i].split('/')[-2]))
           
            dycore = runs[i].split('/')[-2].split('_')[0]
            grid = runs[i].split('/')[-2].split('_')[1]
           
            if(dycore == 'FV3'):
                hist = sorted(glob.glob('{}/run/*h0*regrid*'.format(runs[i])))
            elif(dycore == 'SE'):
                hist = sorted(glob.glob('{}/run/*h0*.nc'.format(runs[i])))
            AOA1_arrs = []
            AOA2_arrs = []
            clock_arrs = []
            u_arrs = []
            time = []
            
            AOA1 = np.array([])
            AOA2 = np.array([])
            clock = np.array([])
            u = np.array([])

            PHI = kk
            LEV = 10
            
            for j in range(len(hist)):
                print('\n---------- working on file {} ----------'.format(hist[j].split('/')[-1]))
                dat = xr.open_dataset(hist[j])
                aoa1dat = dat['AOA1'].sel(lat=PHI, method='nearest').sel(
                                                   lev=LEV, method='nearest').mean('lon')
                aoa2dat = dat['AOA2'].sel(lat=PHI, method='nearest').sel(
                                                   lev=LEV, method='nearest').mean('lon')
                cldat = dat['TT_LW'].sel(lat=PHI, method='nearest').sel(
                                                  lev=LEV, method='nearest').mean('lon')
                udat = dat['U'].sel(lat=PHI, method='nearest').sel(
                                             lev=LEV, method='nearest').mean('lon')
                AOA1 = np.hstack([AOA1, aoa1dat.values])
                AOA2 = np.hstack([AOA2, aoa2dat.values])
                clock = np.hstack([clock, cldat.values])
                u = np.hstack([u, udat.values])
                
                time_arr = dat['time']
                time.extend([t.day + (t.month-1)*30 + (t.year-1)*360 for t in time_arr.values])
             
            time = np.array(time)
           
            AOA2_diff = np.zeros(AOA2.shape)
            clock_diff = np.zeros(clock.shape)
            for k in range(len(time)):
                AOA2_diff[k] = time[k] - AOA2[k]
                clock_diff[k] = time[k] - clock[k]

            ax1.plot(time, u)
            
            ax2.plot(time, clock, label='AOA2', color='m')
            ax2.plot(time, AOA2, label='AOA1', color='r')
            ax2.plot([6000,6000], [8,8], label='TT_LW', color='c')
            ax2.legend(fontsize=11, loc='lower right', frameon=False)
            
            ax3.plot(time, clock_diff / 365, label='TT_LW', color='c')
            ax3.plot(time, AOA2_diff / 365, label='AOA2', color='m')
            ax3.plot(time, AOA1 / 365, label='AOA1', color='r')
            
            #ax1.set_xlim([5500, 8500])
            #ax2.set_xlim([5500, 8500])
            #ax3.set_xlim([5500, 8500])
            #ax2.set_ylim([1500, 6000])
            #ax3.set_ylim([5, 11])

            ax1.set_ylabel('u  [m/s]', fontsize=12)
            ax2.set_ylabel('Clock tracer concentration', fontsize=12)
            ax3.set_ylabel('Age of air  [years]', fontsize=12) 
            ax3.set_xlabel('time [years]', fontsize=12)
            ax2.axes.get_xaxis().set_ticks([])
            ax1.axes.get_xaxis().set_ticks([])

            ax1.set_title('{}_{} at {}deg, 10hPa'.format(dycore, grid, PHI), fontsize=14)
              
            plt.savefig('{}_{}_{}deg_ssw_full.png'.format(dycore, grid, PHI), dpi=300)

            
    


# =============================================================================


def aoa_profiles():
    
    for i in range(len(runs)):

        print('\n\n---------- working on run {} ----------'.format(runs[i].split('/')[-2]))
       
        dycore = runs[i].split('/')[-2].split('_')[0]
        grid = runs[i].split('/')[-2].split('_')[1]
       
        if(dycore == 'FV3'):
            hist = sorted(glob.glob('{}/run/*h0*regrid*'.format(runs[i])))
        elif(dycore == 'SE'):
            hist = sorted(glob.glob('{}/run/*h0*.nc'.format(runs[i])))
        AOA1_arrs = []
        AOA2_arrs = []
        clock_arrs = []
        u_arrs = []
        om_arrs = []
        time = []

        for j in range(len(hist)):
            print('\n---------- working on file {} ----------'.format(hist[j].split('/')[-1]))
            dat = xr.open_dataset(hist[j])
            om_arrs.append(dat['OMEGA'].mean('lon'))
            u_arrs.append(dat['U'].mean('lon'))
            AOA1_arrs.append(dat['AOA1'].mean('lon'))
            AOA2_arrs.append(dat['AOA2'].mean('lon'))
            clock_arrs.append(dat['TT_LW'].mean('lon'))
            time_arr = dat['time']
            time.extend([t.day + (t.month-1)*30 + (t.year-1)*360 for t in time_arr.values])
        om = xr.concat(om_arrs, dim='time')
        u = xr.concat(u_arrs, dim='time')
        AOA1 = xr.concat(AOA1_arrs, dim='time')
        AOA2 = xr.concat(AOA2_arrs, dim='time')
        clock = xr.concat(clock_arrs, dim='time')
        #time = xr.concat(time_arrs, dim='time')
         
        for k in range(len(time)):
            AOA2[k,:,:] = time[k] - AOA2[k,:,:]
            clock[k,:,:] = time[k] - clock[k,:,:]
        omm = om.mean('time')
        um = u.mean('time')
        AOA1m = AOA1.mean('time')
        AOA2m = AOA2.mean('time')
        clockm = clock.mean('time')
        
        fig = plt.figure(figsize=(10,5))
        ax1 = fig.add_subplot(131)
        ax2 = fig.add_subplot(132)
        ax3 = fig.add_subplot(133)
        
        lat = AOA1m['lat'].values
        lon = AOA1m['lev'].values
        
        levels = np.linspace(1, 15, 15)
        ylim = [0.1,1000]

        var_dict = [{'var':AOA1m/365, 'plotType':'contour', 'plotArgs':{'levels':levels}},
                    {'var':um, 'plotType':'contourf', 'colorFormatter':None}]
        pltvert(lat, lon, var_dict, ax=ax1, ylim=ylim)
        
        var_dict[0]['var'] = AOA2m/365
        var_dict[1]['colorArgs']['ax']=ax2
        pltvert(lat, lon, var_dict, ax=ax2, ylim=ylim)
        
        var_dict[1] = {'var':um, 'plotType':'contourf', 'colorArgs':{'label':'u  [m/s]'}}
        var_dict[0]['var'] = clockm/365
        pltvert(lat, lon, var_dict, ax=ax3, ylim=ylim)

        ax1.set_title('AOA1, {} year mean'.format(int(time[-1]/365)))
        ax2.set_title('AOA2, {} year mean'.format(int(time[-1]/365))) 
        ax3.set_title('TT_LW, {} year mean'.format(int(time[-1]/365)))
        
        ax1.set_ylabel('lev  [hPa]', fontsize=12)
        ax1.set_xlabel('lat  [deg]', fontsize=12)
        ax2.set_xlabel('lat  [deg]', fontsize=12)
        ax3.set_xlabel('lat  [deg]', fontsize=12)
            
        plt.savefig('figs/sss_{}_{}.png'.format(dycore, grid), dpi=300)

            
