import pdb
import sys
import glob
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

#VAR = 'T1000'
VAR = 'T050'
    
mc = '/global/cfs/cdirs/m4014/data/HSW/outputs/release_011423/netcdf/ens_stats_latlon/mean_climate.regrid.91x180_bilinear.nc'
mcd = xr.open_dataset(mc)
weights = np.cos(np.deg2rad(mcd['lat']))
weights.name = 'weights'
mcd = mcd[VAR].weighted(weights)
mcd = mcd.mean(['lon', 'lat'])


# ----------------------- LOW VAR JW1.00X FIXED-------------------
do_lowvar_jw_fixed = 1
N = 5

if(do_lowvar_jw_fixed):


    loc = '/pscratch/sd/j/jhollo/E3SM/E3SMv2_cases/sai_cases/verification_for_profiling/jw_fixed_runs'
    ensxm = []

    time0=xr.open_dataset(sorted(glob.glob('{}/*h2*linear.nc'.format(loc)))[0])['ndcur'].values
    time1=xr.open_dataset(sorted(glob.glob('{}/*h2*linear.nc'.format(loc)))[1])['ndcur'].values
    time = np.hstack([time0, time1])

    # just consider first 5 ensemble members
    for i in range(N):
        
        print('working on ens{:02d}...'.format(i+1))
        tmp = './tmp/ensxm_{:02d}_zonalmean_jw_fixed.nc'.format(i+1)
        try:
            print('reading...')
            ensx = xr.open_dataset(tmp)
        except FileNotFoundError:
            print('computing...')
            files=[xr.open_dataset(d) for d in 
                   sorted(glob.glob('{}/*ens{:02d}*h2*linear.nc'.format(loc,i+1)))]
            ensx = xr.concat(files, dim='time', data_vars=['T1000', 'T050']) 
            ensx = ensx.mean('lon')
            print('writing')
            ensx.to_netcdf(tmp)
       
        weights = np.cos(np.deg2rad(ensx['lat']))
        weights.name = 'weights'
        ensx = ensx[VAR].weighted(weights)
        ensx = ensx.mean('lat')
        ensxm.append(ensx)

        if(i == 0):
            ens_mean = xr.zeros_like(ensxm[0])
        ens_mean = ens_mean + ensxm[i]
    ens_mean = ens_mean / N

    # compute std
    sum_diffsq = xr.zeros_like(ens_mean)
    for i in range(N):
        diffsq = (ensxm[i] - ens_mean)**2
        sum_diffsq = sum_diffsq + diffsq
    std = (sum_diffsq / N) ** (1/2)

    fig_jw_fixed = plt.figure(figsize=(6,3))
    ax = fig_jw_fixed.add_subplot(111)
    ax.plot(time, ens_mean, '-', color='yellowgreen')
    ax.fill_between(time, ens_mean-std, ens_mean + std, color='yellowgreen', alpha=0.33)
    ax.grid()
    ax.set_xlabel('time [days]')
    ax.set_ylabel('{} [K]'.format(VAR))
    #ax.set_title('030123')
    ax.set_title('updated 1.00X ensemble mean {} response'.format(VAR))
    plt.tight_layout()
    plt.savefig('{}_lowvar_std_jw_fixed.png'.format(VAR), dpi=300)




# ----------------------- TEST JH PERLMUTTER -------------------
do_lowvar_pm = 0
N = 1

if(do_lowvar_pm):


    loc = '/pscratch/sd/j/jhollo/E3SM/E3SMv2_cases/sai_cases/tests/HSW_SAI_ne16pg2_L72_300day_PERLMUTTER_INJECTION_TEST/run'
    ensxm = []

    time = xr.open_dataset(sorted(glob.glob('{}/*h2*linear.nc'.format(loc)))[0])['ndcur'].values

    for i in range(N):
        
        print('working on ens{:02d}...'.format(i+1))
        tmp = './tmp/ensxm_{:02d}_zonalmean_pm.nc'.format(i+1)
        try:
            print('reading...')
            ensx = xr.open_dataset(tmp)
        except FileNotFoundError:
            print('computing...')
            files=[xr.open_dataset(d) for d in 
                   sorted(glob.glob('{}/*h2*linear.nc'.format(loc)))]
            
            ensx = files[0] 
            ensx = ensx.mean('lon')
            print('writing')
            ensx.to_netcdf(tmp)
       
        weights = np.cos(np.deg2rad(ensx['lat']))
        weights.name = 'weights'
        ensx = ensx[VAR].weighted(weights)
        ensx = ensx.mean('lat')
        ensxm.append(ensx)

        if(i == 0):
            ens_mean = xr.zeros_like(ensxm[0])
        ens_mean = ens_mean + ensxm[i]
    ens_mean = ens_mean / N

    # compute std
    sum_diffsq = xr.zeros_like(ens_mean)
    for i in range(N):
        diffsq = (ensxm[i] - ens_mean)**2
        sum_diffsq = sum_diffsq + diffsq
    std = (sum_diffsq / N) ** (1/2)

    fig_pm = plt.figure(figsize=(6,3))
    ax = fig_pm.add_subplot(111)
    ax.plot(time, ens_mean, '-', color='yellowgreen')
    ax.fill_between(time, ens_mean-std, ens_mean + std, color='yellowgreen', alpha=0.33)
    ax.grid()
    ax.set_xlabel('time [days]')
    ax.set_ylabel('{} [K]'.format(VAR))
    #ax.set_title('030123')
    ax.set_title('PM test ensemble mean surface response')
    plt.tight_layout()
    plt.savefig('T1000_lowvar_std_pm.png', dpi=300)


# ----------------------- LOW VAR AS1.20X-------------------
do_lowvar_as = 0
N = 5

if(do_lowvar_as):


    loc = '/pscratch/sd/j/jhollo/E3SM/E3SMv2_cases/sai_cases/asteyer_hsw_variability_runs'
    ensxm = []

    time0=xr.open_dataset(sorted(glob.glob('{}/*h2*linear.nc'.format(loc)))[0])['ndcur'].values
    time1=xr.open_dataset(sorted(glob.glob('{}/*h2*linear.nc'.format(loc)))[1])['ndcur'].values
    time = np.hstack([time0, time1])

    # just consider first 5 ensemble members
    for i in range(N):
        
        #if i+1 != int(sys.argv[1]): continue

        print('working on ens{:02d}...'.format(i+1))
        tmp = './tmp/ensxm_{:02d}_zonalmean_as.nc'.format(i+1)
        try:
            print('reading...')
            ensx = xr.open_dataset(tmp)
        except FileNotFoundError:
            print('computing...')
            files=[xr.open_dataset(d) for d in 
                   sorted(glob.glob('{}/*ens{:02d}*h2*linear.nc'.format(loc,i+1)))]
            
            for j in range(len(files)):
                files[j] = files[j].assign(T1000 = files[j]['T'].sel({'lev':1000}, method='nearest'))
                files[j] = files[j].assign(T050 = files[j]['T'].sel({'lev':50}, method='nearest'))
        
            ensx = xr.concat(files, dim='time', data_vars=['T1000', 'T050']) 
            ensx = ensx.mean('lon')
            print('writing')
            ensx.to_netcdf(tmp)
            #quit()
       
        weights = np.cos(np.deg2rad(ensx['lat']))
        weights.name = 'weights'
        ensx = ensx[VAR].weighted(weights)
        ensx = ensx.mean('lat')
        ensxm.append(ensx)

        if(i == 0):
            ens_mean = xr.zeros_like(ensxm[0])
        ens_mean = ens_mean + ensxm[i]
    ens_mean = ens_mean / N

    # compute std
    sum_diffsq = xr.zeros_like(ens_mean)
    for i in range(N):
        diffsq = (ensxm[i] - ens_mean)**2
        sum_diffsq = sum_diffsq + diffsq
    std = (sum_diffsq / N) ** (1/2)

    fig_as = plt.figure(figsize=(6,3))
    ax = fig_as.add_subplot(111)
    ax.plot(time, ens_mean, '-', color='yellowgreen')
    ax.fill_between(time, ens_mean-std, ens_mean + std, color='yellowgreen', alpha=0.33)
    ax.grid()
    ax.set_xlabel('time [days]')
    ax.set_ylabel('{} [K]'.format(VAR))
    #ax.set_title('030123')
    ax.set_title('low variability 1.20X ensemble mean surface response')
    plt.tight_layout()
    plt.savefig('T1000_lowvar_std_as.png', dpi=300)


# ----------------------- LOW VAR JW1.00X-------------------
do_lowvar_jw = 0
N = 5

if(do_lowvar_jw):


    loc = '/global/cfs/cdirs/m4014/data/HSW/outputs/release_030123/netcdf/low_var_mass_ens_latlon'
    ensxm = []

    time0=xr.open_dataset(sorted(glob.glob('{}/ens01_mass1.00X/*h2*'.format(loc)))[0])['ndcur'].values
    time1=xr.open_dataset(sorted(glob.glob('{}/ens01_mass1.00X/*h2*'.format(loc)))[1])['ndcur'].values
    time = np.hstack([time0, time1])

    # just consider first 5 ensemble members
    for i in range(N):

        print('working on ens{:02d}...'.format(i+1))
        tmp = './tmp/ensxm_{:02d}_zonalmean_jw.nc'.format(i+1)
        try:
            print('reading...')
            ensx = xr.open_dataset(tmp)
        except FileNotFoundError:
            print('computing...')
            files=[xr.open_dataset(d) for d in 
                   sorted(glob.glob('{}/ens{:02d}_mass1.00X/*h2*lin.nc'.format(loc,i+1)))]
            ensx = xr.concat(files, dim='time', data_vars=['T1000', 'T050']) 
            ensx = ensx.mean('lon')
            print('writing')
            ensx.to_netcdf(tmp)
       
        weights = np.cos(np.deg2rad(ensx['lat']))
        weights.name = 'weights'
        ensx = ensx[VAR].weighted(weights)
        ensx = ensx.mean('lat')
        ensxm.append(ensx)

        if(i == 0):
            ens_mean = xr.zeros_like(ensxm[0])
        ens_mean = ens_mean + ensxm[i]
    ens_mean = ens_mean / N

    # compute std
    sum_diffsq = xr.zeros_like(ens_mean)
    for i in range(N):
        diffsq = (ensxm[i] - ens_mean)**2
        sum_diffsq = sum_diffsq + diffsq
    std = (sum_diffsq / N) ** (1/2)

    fig_jw = plt.figure(figsize=(6,3))
    ax = fig_jw.add_subplot(111)
    ax.plot(time, ens_mean, '-', color='yellowgreen')
    ax.fill_between(time, ens_mean-std, ens_mean + std, color='yellowgreen', alpha=0.33)
    ax.grid()
    ax.set_xlabel('time [days]')
    ax.set_ylabel('{} [K]'.format(VAR))
    #ax.set_title('030123')
    ax.set_title('low variability 1.00X ensemble mean surface response')
    plt.tight_layout()
    plt.savefig('T1000_lowvar_std_jw.png', dpi=300)



# ----------------------- LOW VAR JH -------------------
do_lowvar = 1
N = 5

if(do_lowvar):


    loc = '/global/cfs/cdirs/m4014/data/HSW/outputs/release_030123/netcdf/low_var_ens_latlon'
    ensxm = []

    time0 = xr.open_dataset(sorted(glob.glob('{}/ens01/*h2*'.format(loc)))[0])['ndcur'].values
    time1 = xr.open_dataset(sorted(glob.glob('{}/ens01/*h2*'.format(loc)))[1])['ndcur'].values
    time = np.hstack([time0, time1])

    for i in range(N):
        

        print('working on ens{:02d}...'.format(i+1))
        tmp = './tmp/ensxm_{:02d}_zonalmean.nc'.format(i+1)
        try:
            print('reading...')
            ensx = xr.open_dataset(tmp)
        except FileNotFoundError:
            print('computing...')
            files=[xr.open_dataset(d) for d in sorted(glob.glob('{}/ens{:02d}/*h2*'.format(loc,i+1)))]
            ensx = xr.concat(files, dim='time', data_vars=['T1000', 'T050']) 
            ensx = ensx.mean('lon')
            print('writing')
            ensx.to_netcdf(tmp)
       
        weights = np.cos(np.deg2rad(ensx['lat']))
        weights.name = 'weights'
        ensx = ensx[VAR].weighted(weights)
        ensx = ensx.mean('lat')
        ensxm.append(ensx)

        if(i == 0):
            ens_mean = xr.zeros_like(ensxm[0])
        ens_mean = ens_mean + ensxm[i]
    ens_mean = ens_mean / N
    ens_mean = ens_mean

    # compute std
    sum_diffsq = xr.zeros_like(ens_mean)
    for i in range(N):
        diffsq = (ensxm[i] - ens_mean)**2
        sum_diffsq = sum_diffsq + diffsq
    std = (sum_diffsq / N) ** (1/2)

    fig = plt.figure(figsize=(6,3))
    ax = fig.add_subplot(111)
    ax.plot(time, ens_mean, '-', color='yellowgreen')
    ax.fill_between(time, ens_mean-std, ens_mean + std, color='yellowgreen', alpha=0.33)
    ax.grid()
    ax.set_xlabel('time [days]')
    ax.set_ylabel('{} [K]'.format(VAR))
    #ax.set_title('030123')
    ax.set_title('original low-variability ensemble mean {} response'.format(VAR))
    plt.tight_layout()
    plt.savefig('{}_lowvar_std.png'.format(VAR), dpi=300)


# ----------------------- HIGH VAR -------------------
do_hivar = 0

if(do_hivar):

    mc = '/global/cfs/cdirs/m4014/data/HSW/outputs/release_011423/netcdf/ens_stats_latlon/mean_climate.regrid.91x180_bilinear.nc'
    loc = '/global/cfs/cdirs/m4014/data/HSW/outputs/release_011423/netcdf/ens_members_latlon'
    ensxm = []

    time = xr.open_dataset(sorted(glob.glob('{}/ens01/*h2*near.nc'.format(loc)))[0])['ndcur'].values

    for i in range(5):
        
        print('working on ens{:02d}...'.format(i+1))
        print('computing...')
        ensx=xr.open_dataset(glob.glob('{}/ens{:02d}/*h2*near.nc'.format(loc,i+1))[0])[VAR]
        weights = np.cos(np.deg2rad(ensx['lat']))
        weights.name = 'weights'
        ensx = ensx.weighted(weights)
        ensx = ensx.mean(['lon','lat'])
        ensxm.append(ensx)

        if(i == 0):
            ens_mean = xr.zeros_like(ensxm[0])
        ens_mean = ens_mean + ensxm[i]
    ens_mean = ens_mean / 5
    ens_mean = ens_mean

    # compute std
    sum_diffsq = xr.zeros_like(ens_mean)
    for i in range(5):
        diffsq = (ensxm[i] - ens_mean)**2
        sum_diffsq = sum_diffsq + diffsq
    std = (sum_diffsq / 5) ** (1/2)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(time, ens_mean, '-', color='`yellowgreen')
    ax.fill_between(time, ens_mean-std, ens_mean+std, color='yellowgreen', alpha=0.33)
    ax.grid()
    ax.set_xlabel('time [days]')
    ax.set_ylabel('{} [K]'.format(VAR))
    ax.set_title('011423')

#plt.show()
