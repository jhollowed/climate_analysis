import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import pdb
import climate_toolbox as ctb

def heating_comp(dat1f, dat2f):
    
    dat1 = xr.open_dataset(dat1f)
    dat2 = xr.open_dataset(dat2f)
    #datdiff = dat1['T'] - dat2['T']
    #datdiff.to_netcdf('Tdiff.nc')

    colors = plt.cm.plasma([0.2, 0.4, 0.6, 0.8])
    psel = [30, 50, 70, 100]
    fig = plt.figure(figsize=(8, 4))
    ax = fig.add_subplot(111)

    for i in range(len(psel)):
        p = psel[i]

        # zonal mean at pressure level
        T1 = dat1['T'].sel({'lev':p}, method='nearest').mean('lon') 
        T2 = dat2['T'].sel({'lev':p}, method='nearest').mean('lon') 

        # mean near tropics
        T1 = T1.sel({'lat':slice(-30, 30)})
        T2 = T2.sel({'lat':slice(-30, 30)})

        weights1 = np.cos(np.deg2rad(T1.lat))
        weights1.name = 'weights'
        weights2 = np.cos(np.deg2rad(T2.lat))
        weights2.name = 'weights'
        T1_weighted = T1.weighted(weights1)
        T2_weighted = T2.weighted(weights2)
        T1 = T1_weighted.mean('lat')
        T2 = T2_weighted.mean('lat')
        #T1 = T1.mean('lat')
        #T2 = T2.mean('lat')
        
        Tdiff = T1-T2
        t = ctb.time2day(T1['time'].values)
        ax.plot(t, -Tdiff, label='{} hPa'.format(p), color=colors[i])
        ax.set_xlim([0, 30])
        ax.set_xlabel('time [days]', fontsize=14)
        ax.set_ylabel('Temperature difference [K]', fontsize=14)
    ax.legend()
    plt.savefig('./figs/stratHeating.png', dpi=150)
    #plt.show()

    




if __name__ == '__main__':

    datout = '/global/cscratch1/sd/jhollo/E3SM/E3SMv2_cases/sai_cases/E3SM_ne16_L72_FIDEAL_SAI'
    dat1 = '{}_90days_stratHeating/run/'\
           'E3SM_ne16_L72_FIDEAL_SAI_90days_stratHeating.eam.h0.0001-01-01-00000.regrid.2x2.nc'\
           .format(datout)
    dat2 = '{}/run/E3SM_ne16_L72_FIDEAL_SAI.eam.h0.0001-01-01-00000.regrid.2x2.nc'.format(datout)

    heating_comp(dat1, dat2)
