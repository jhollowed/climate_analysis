'''
Joe Hollowed Dec 2024
This script prepares the limvar data for computation of the RCTT.
This requires each limvar ensemble member to be concatenated with 
redicual circulation velocities from the simulation(s) that provided
the initialization files for the run. This script concats:
    - 8 years of limvar
    - 3 years of fullvar ens1, wich provided the limvar IC
    - 7 years of the CLDERA historical run, which provided the fullvar IC
This allows us to construct backward-trajectories through the residual
circulation for 10-years prior to the limvar initalization (and the eruptions), 
which will thus allow the RCTT to be known in the meridional plane for the
entire limvar time period.

usage:
RCTT_on_combines_limvar_fullvar.py {SO2 mass in Tg} {ensemble number}

a mass of 0Tg indicates the counterfactual simulations
'''
import sys
import pdb
import glob
import numpy as np
import xarray as xr
from datetime import datetime
from datetime import timedelta
import matplotlib.pyplot as plt
from mpi4py import MPI

# local imports
sys.path.insert(1, '/global/homes/j/jhollo/repos/RCTT')
from rctt import RCTT
    
# setup MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nranks = comm.Get_size()
printt = lambda s,end=None: print(s, end=end) if rank == 0 else None
    

def run(mass, ens):
    

    print('--- rank {} starting'.format(rank))
    sys.stdout.flush()
    comm.Barrier()

    # convenience functions
    to_datetime = lambda times: [datetime(t.year, t.month, t.day) for t in times.values]

    # get limvar data
    printt('reading limvar data for {}Tg ens{}...'.format(mass, ens))
    if(mass == 0): config = 'ens{}.cf.'.format(ens)
    else:          config = '{}Tg.ens{}.'.format(mass, ens)
    lv_loc  = '/pscratch/sd/j/jhollo/E3SM/historical_data/limvar/limvar_monthly'
    lv_temf = glob.glob('{}/*{}*TEM*_L45_monthlymean.nc'.format(lv_loc, config))[0]
    lv_zmf  = glob.glob('{}/*{}*zonalmeans_monthlymean.nc'.format(lv_loc, config))[0]
    lv_trop = xr.open_dataset(lv_zmf)['TROP_P']
    lv_tem  = xr.open_dataset(lv_temf)
    lv_time = lv_trop.time
    lv_vtem, lv_wtem = lv_tem['vtem'], lv_tem['wtem']

    # get fullvar data
    printt('reading fullvar data...')
    fv_loc  = '/pscratch/sd/j/jhollo/E3SM/historical_data/limvar/fullvar_monthly'
    fv_temf = glob.glob('{}/*TEM*_L45_monthlymean.nc'.format(fv_loc))[0]
    fv_zmf  = glob.glob('{}/*zonalmeans_monthlymean.nc'.format(fv_loc))[0]
    fv_trop = xr.open_dataset(fv_zmf)['TROP_P']
    fv_tem  = xr.open_dataset(fv_temf)
    fv_time = fv_trop.time
    fv_vtem, fv_wtem = fv_tem['vtem'], fv_tem['wtem']

    # get historical data
    printt('reading historical data...')
    hist_loc  = '/pscratch/sd/j/jhollo/E3SM/historical_data/limvar/historical_monthly'
    hist_temf = glob.glob('{}/*TEM*_L45_monthlymean.nc'.format(hist_loc))[0]
    hist_zmf  = glob.glob('{}/*zonalmeans_monthlymean.nc'.format(hist_loc))[0]
    hist_trop = xr.open_dataset(hist_zmf)['TROP_P']
    hist_tem  = xr.open_dataset(hist_temf)
    hist_time = hist_trop.time
    hist_vtem, hist_wtem = hist_tem['vtem'], hist_tem['wtem']

    # now shift fullvar and historical calendar times...
    # the limvar initial condition, which in limvar is 1991-6-1, originated from fullvar 1988-6-1
    # thus, we will shift the fullvar and historical calendars forward by 3 years
    printt('shifting calendars...')
    fv_time   = np.array([type(t)(t.year+3, t.month, t.day) for t in fv_time.values])
    fv_time   = xr.DataArray(fv_time, coords={'time':fv_time})
    hist_time = np.array([type(t)(t.year+3, t.month, t.day) for t in hist_time.values])
    hist_time = xr.DataArray(hist_time, coords={'time':hist_time})
    # assign the shfited times to the data arrays
    fv_vtem   = fv_vtem.assign_coords(time=fv_time)
    fv_wtem   = fv_wtem.assign_coords(time=fv_time)
    fv_trop   = fv_trop.assign_coords(time=fv_time)
    hist_vtem = hist_vtem.assign_coords(time=hist_time)
    hist_wtem = hist_wtem.assign_coords(time=hist_time)
    hist_trop = hist_trop.assign_coords(time=hist_time)

    # remove duplicate times from dataset edges
    printt('removing time duplicates from dataset edges...')
    fv_vtem   = fv_vtem.sel(time=slice(fv_time[0], lv_time[0]-timedelta(days=1)))
    fv_wtem   = fv_wtem.sel(time=slice(fv_time[0], lv_time[0]-timedelta(days=1)))
    fv_trop   = fv_trop.sel(time=slice(fv_time[0], lv_time[0]-timedelta(days=1)))
    hist_vtem = hist_vtem.sel(time=slice(hist_time[0], fv_time[0]-timedelta(days=1)))
    hist_wtem = hist_wtem.sel(time=slice(hist_time[0], fv_time[0]-timedelta(days=1)))
    hist_trop = hist_trop.sel(time=slice(hist_time[0], fv_time[0]-timedelta(days=1)))
    fv_time   = fv_vtem.time
    hist_time = hist_vtem.time

    # sanity check
    debug=False
    if(debug and rank==0):
        printt('plotting sanity check...')
        fig, ax = plt.subplots(3)
        vargs = {'lat':45, 'plev':100, 'method':'nearest'}
        targs = {'lat':45, 'method':'nearest'}
        ax[0].plot(to_datetime(lv_time), lv_vtem.sel(**vargs), '-k')
        ax[0].plot(to_datetime(fv_time), fv_vtem.sel(**vargs), '-r')
        ax[0].plot(to_datetime(hist_time), hist_vtem.sel(**vargs), '-b')
        ax[0].set_ylabel('vtem')
        ax[1].plot(to_datetime(lv_time), lv_wtem.sel(**vargs), '-k')
        ax[1].plot(to_datetime(fv_time), fv_wtem.sel(**vargs), '-r')
        ax[1].plot(to_datetime(hist_time), hist_wtem.sel(**vargs), '-b')
        ax[0].set_ylabel('wtem')
        ax[2].plot(to_datetime(lv_time), lv_trop.sel(**targs), '-k')
        ax[2].plot(to_datetime(fv_time), fv_trop.sel(**targs), '-r')
        ax[2].plot(to_datetime(hist_time), hist_trop.sel(**targs), '-b') 
        ax[0].set_ylabel('trop_p')
        plt.show()
    comm.Barrier()

    # concatenate limvar, fullvar, historical in time
    printt('concatenating datasets in time...')
    vtem = xr.concat([hist_vtem, fv_vtem, lv_vtem], dim='time')
    wtem = xr.concat([hist_wtem, fv_wtem, lv_wtem], dim='time')
    trop = xr.concat([hist_trop, fv_trop, lv_trop], dim='time')
    time = vtem.time
    
    # replace tropopause with flat 700 hPa boundary
    surface_trop = True
    if(surface_trop):
        trop = trop/trop * (700*100)

    # configure RCTT launch points
    launch_lats  = np.arange(-86, 87, 2)
    launch_lats  = xr.DataArray(launch_lats, coords={'lat':launch_lats})
    launch_plev  = np.logspace(0.5, 2.5, 20)
    launch_plev  = xr.DataArray(launch_plev, coords={'plev':launch_plev})

    # configure RCTT launch times; disribute amongst ranks
    sys.stdout.flush()
    comm.Barrier()
    launch_times = time.sel(time=time.time.dt.year.isin(np.linspace(1991,2000,10)))
    launch_times = launch_times.sel(time=launch_times.time.dt.month.isin([1,4,7,10]))
    ntimes = len(launch_times)
    launch_times = np.array_split(launch_times, nranks)
    launch_times = launch_times[rank]
    print('--- rank {} has {}/{} launch times; {}--{}'.format(rank, len(launch_times), ntimes, 
          launch_times.values[0].strftime("%Y-%m-%d"),
          launch_times.values[-1].strftime("%Y-%m-%d")))
    sys.stdout.flush()
    comm.Barrier()

    # call RCTT
    printt('calling RCTT...')
    outdir    = '/pscratch/sd/j/jhollo/E3SM/historical_data/limvar/RCTT/rank_outputs/'
    outprefix = 'limvar_{}Tg_ens{}_fvlvconcat_surfTrop{}_rank{}'.format(
                                      mass, ens, int(surface_trop), rank)
    age_limit, resday = 10, 5
    rctt = RCTT(vtem, wtem, trop, outdir=outdir, outprefix=outprefix, quiet=rank!=0)
    rctt.launch(launch_lats, launch_plev, launch_times, 
                overwrite=False, resday=resday, age_limit=age_limit)
    print('--- rank {} finished'.format(rank))
    sys.stdout.flush()
    comm.Barrier()

    # concatenate outputs of all ranks to single files
    if(rank == 0):
        rank_outputs = np.array(glob.glob('{}/{}*'.format(outdir, outprefix.split('_rank')[0])))
        rank_rctt_outputs = np.array([f for f in rank_outputs if 'RCTT' in f.split('/')[-1]])
        rank_traj_outputs = np.array([f for f in rank_outputs if 'Traj' in f.split('/')[-1]])
        rank_rctt_sort = np.argsort([int(f.split('/')[-1].split('_rank')[-1].split('_')[0]) \
                                     for f in rank_rctt_outputs])
        rank_traj_sort = np.argsort([int(f.split('/')[-1].split('_rank')[-1].split('_')[0]) \
                                     for f in rank_rctt_outputs])
        rank_rctt_outputs = rank_rctt_outputs[rank_rctt_sort]
        rank_traj_outputs = rank_traj_outputs[rank_traj_sort]
        print('concatenating RCTT, trajectory files on rank 0...')
        rctt_concat = xr.concat([xr.open_dataset(f) for f in rank_rctt_outputs], dim='time')
        traj_concat = xr.concat([xr.open_dataset(f) for f in rank_traj_outputs], dim='time')
        concat_times = [ti.strftime("%Y-%m-%d") for ti in rctt_concat.time.values]
        print('writing out concatenating files on rank 0...')
        outdir    = '/pscratch/sd/j/jhollo/E3SM/historical_data/limvar/RCTT/concat_outputs'
        prefix = outprefix.split('_rank')[0]
        rctt_concat_file = '{}/{}_RCTT_{}--{}_ageLimit{}_res{}.nc'.format(outdir, prefix, 
                                    concat_times[0], concat_times[-1], age_limit, resday)
        traj_concat_file = '{}/{}_RCTT_{}--{}_ageLimit{}_res{}.nc'.format(outdir, prefix, 
                                    concat_times[0], concat_times[-1], age_limit, resday)
        rctt_concat.to_netcdf(rctt_concat_file)
        traj_concat.to_netcdf(traj_concat_file)
    comm.Barrier()
    printt('all done')

    if(rank==0 and 0):
        tx, ty = traj['trajectories_lat'], traj['trajectories_plev']
        printt(np.nanmin((ttimes/365).values))
        printt(np.nanmax((ttimes/365).values))
        cc = plt.contourf(ttimes.lat, ttimes.plev, ttimes[-1].T/365, 
                          cmap='viridis', levels=np.arange(10), extend='both')
        targs = {'lat':76, 'plev':40, 'method':'nearest'}
        plt.plot(tx.sel(**targs), ty.sel(**targs), '-r')
        plt.gca().set_yscale('log')
        plt.gca().invert_yaxis()
        plt.colorbar(cc)
        plt.show()
        pdb.set_trace()

if(__name__ == '__main__'): 
    # command line args
    ens = sys.argv[1] 
    mass = sys.argv[2]
    run(mass, ens)
    #for mass in [10]:
    #    run(mass, ens)
