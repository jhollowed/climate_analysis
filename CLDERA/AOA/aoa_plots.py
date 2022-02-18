import numpy as np
import matplotlib.pyplot as plt
import glob
import pdb
import xarray as xr

runs = glob.glob('/glade/scratch/jhollowed/CAM/cases/aoa_runs/project3/*aoa/')

for i in range(len(runs)):

    print('\n\n---------- working on run {} ----------'.format(runs[i].split('/')[-1]))
   
    dycore = runs[i].split('/')[-2].split('_')[0]
    grid = runs[i].split('/')[-2].split('_')[1]
   
    if(dycore == 'FV3'):
        hist = sorted(glob.glob('{}/run/*h0*regrid*'.format(runs[i])))
    elif(dycore == 'SE'):
        hist = sorted(glob.glob('{}/run/*h0*.nc'.format(runs[i])))
    AOA1_arrs = []
    AOA2_arrs = []
    clock_arrs = []
    time = []

    for j in range(len(hist)):
        print('\n---------- working on file {} ----------'.format(hist[j].split('/')[-1]))
        dat = xr.open_dataset(hist[i])
        AOA1_arrs.append(dat['AOA1'].mean('lon'))
        AOA2_arrs.append(dat['AOA2'].mean('lon'))
        clock_arrs.append(dat['TT_LW'].mean('lon'))
        time_arr = dat['time']
        time.extend([t.day + (t.month-1)*30 + ((t.year+(4*j))-1)*360 for t in time_arr.values])
    AOA1 = xr.concat(AOA1_arrs, dim='time')
    AOA2 = xr.concat(AOA2_arrs, dim='time')
    clock = xr.concat(clock_arrs, dim='time')
    #time = xr.concat(time_arrs, dim='time')
    
    for k in range(len(time)):
        AOA2[k,:,:] = time[k] - AOA2[k,:,:]
        clock[k,:,:] = time[k] - clock[k,:,:]

    fig = plt.figure(figsize=(8,5))
    ax1 = fig.add_subplot(131)
    ax2 = fig.add_subplot(132)
    ax3 = fig.add_subplot(133)
    
    AOA1m = AOA1.mean('time')
    AOA2m = AOA2.mean('time')
    clockm = clock.mean('time')
    
    levels = np.linspace(1, 10, 10)
    cmap = plt.cm.jet
    LAT, LEV = np.meshgrid(AOA1m['lat'].values, AOA1m['lev'].values)
    levels=10

    CS = ax1.contour(LAT, LEV, AOA1m/360, levels=levels)
    ax1.clabel(CS, CS.levels, inline=True, fontsize=10)
    ax1.set_title('AOA1, {} year mean'.format(int(time[-1]/360)))
    ax1.set_ylim([1, 1000])
    
    CS=ax2.contour(LAT, LEV, AOA2m/360, levels=levels)
    ax2.clabel(CS, CS.levels, inline=True, fontsize=10)
    ax2.set_title('AOA1, {} year mean'.format(int(time[-1]/360)))
    ax2.set_ylim([1, 1000])
    
    CS=ax3.contour(LAT, LEV, clockm/360, levels=levels)
    ax3.clabel(CS, CS.levels, inline=True, fontsize=10)
    ax3.set_title('AOA1, {} year mean'.format(int(time[-1]/360)))
    ax3.set_ylim([1, 1000])
    
    ax1.set_ylabel('Clock tracer value [years]', fontsize=12)
    ax1.set_xlabel('lat  [deg]', fontsize=12)
    ax2.set_xlabel('lat  [deg]', fontsize=12)
    ax3.set_xlabel('lat  [deg]', fontsize=12)
    ax1.invert_yaxis()
    ax2.invert_yaxis()
    ax3.invert_yaxis()
    ax1.set_yscale('log')
    ax2.set_yscale('log')
    ax3.set_yscale('log')
        
    #plt.show()
    plt.savefig('sss_{}_{}.png'.format(dycore, grid), dpi=300)

        
