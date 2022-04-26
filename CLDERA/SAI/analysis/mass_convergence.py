# eruption_movies.py
# Joe Hollowed
# CLIMATE 589, 4/15/2022
#
# renders video frames of the eruption in horizontal cross, vertical cross, 
# and AzimuthalEquidistant projection

import numpy as np
import xarray as xr
import matplotlib as mpl
import artist_utils as claut
import climate_toolbox as ctb
import matplotlib.pyplot as plt
import pdb

params = {'text.usetex': False, 'mathtext.fontset': 'stixsans'}
plt.rcParams.update(params)

# ============================================================


def mass_convergence(runf, title, savedest):
   

    # read data, transform time coordinate to ndays
    dat = xr.open_dataset(runf)
    td = ctb.time2day(dat['time'])
    mass_atm = (dat['SAI_MASS']).sum(['lon', 'lat', 'lev'])
    mass_so2 = (dat['SAI_MASS'] * dat['SAI_SO2']).sum(['lon', 'lat', 'lev'])
    mass_ash = (dat['SAI_MASS'] * dat['SAI_ASH']).sum(['lon', 'lat', 'lev'])

    MSO2 = 2e10
    Mash = 2e10
     
    # ---------- plot ----------
    print('\n\n=============== {}'.format(title)) 

    fig = plt.figure(figsize=(5,6))
    ax = fig.add_subplot(211)
    ax.plot(td, mass_atm / 1e19, '-k')
    ax.set_ylabel(r'mass [kg] $\times 10^{-19}$', fontsize=12)
    ax.tick_params(axis="x", direction='in')
    ax.set_xticklabels([])
    ax.set_title('{}'.format(title), fontsize=14)
    
    ax = fig.add_subplot(212)
    ax.plot(td, mass_so2/MSO2, '-b', label='SO2')
    #ax.plot(td, mass_ash, '-r', label='ash mass')
    ax.plot(td, np.ones(len(td)), '--k', lw=0.75)
    #ax.plot(td, np.ones(len(td))*Mash, '--r', label='expected total ash mass') 
    ax.plot([0,0], [0,0], '-k', label='total mass of atmosphere')
    ax.legend(frameon=False)
    ax.set_xlabel('time [days]', fontsize=12)
    ax.set_ylabel('tracer mass / true tracer mass', fontsize=12)
    ax2 = ax.secondary_xaxis('top')
    ax2.tick_params(axis='x', direction='in')
    ax2.set_xticklabels([])

    plt.tight_layout()
    plt.subplots_adjust(hspace=0)
    
    

    plt.show()
    pdb.set_trace()
    plt.savefig('{}/{}_t{}.png'.format(savedest, title, tsnap), dpi=150)
        





if(__name__ == '__main__'):

    inp = '/glade/scratch/jhollowed/CAM/cases/sai_runs/SE_ne16L72_whs_sai_fix0_tau0_nsplit1_nodiff1'\
          '/run/SE_ne16L72_whs_sai_fix0_tau0_nsplit1_nodiff1.cam.h0.0001-01-01-00000.nc'
    mass_convergence(inp, 'SE ne30L72, CAM HSW \nno diffusion, constant injection for 48 hr', '.')
    
