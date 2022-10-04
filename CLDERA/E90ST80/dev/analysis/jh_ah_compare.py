import numpy as np
import matplotlib.pyplot as plt
import metpy.constants as const
from metpy.units import units as u
import xarray as xr
import pdb
from climate_artist import vertical_slice as pltvert
from climate_artist import horizontal_slice as plthor
import artist_utils as aut
import cartopy.crs as ccrs


def compare_q(var, f1, l1, l2, trange='30day', native=None, varIs2d=False, unit='ppb'):

    if(trange == '30day'):
        if(var == 'E90'):
            levels = np.linspace(0, 40, 5)
            levdiff = np.linspace(0.1, 3.1, 7)
            cdiff = 'b'
            cmap = 'YlGnBu'
            levslice = {'lev':slice(39, 1000)}
        elif(var == 'ST80_25'):
            levels = np.linspace(0, 200, 11)
            levdiff = np.linspace(-0.6, 1, 9)
            cdiff = 'r'
            cmap = 'YlOrRd'
            levslice = {'lev':slice(55, 160)}
        elif(var == 'SFE90'):
            levels = 10
            levdiff = 10
            cdiff = 'b'
            cmap = 'YlGnBu'
    elif(trange == '1year'):
        if(var == 'E90'):
            levels = np.linspace(0, 150, 13)
            levdiff = np.linspace(0, 5, 11)
            cdiff = 'b'
            cmap = 'YlGnBu'
            levslice = {'lev':slice(20, 800)}
        elif(var == 'ST80_25'):
            levels = np.linspace(0, 200, 11)
            levdiff = np.linspace(-0.6, 1, 9)
            cdiff = 'r'
            cmap = 'YlOrRd'
            levslice = {'lev':slice(55, 160)}
        elif(var == 'SFE90'):
            levels = 10
            levdiff = 10
            cdiff = 'b'
            cmap = 'YlGnBu'


    print('-------- reading data')
    d = xr.open_dataset(f1)
    time = d['ndcur']
    lat = d['lat']
    lon = d['lon']
    varj = '{}j'.format(var)

    print('-------- taking diff')
    if(not varIs2d): 
        maxdim_n = ['ncol','lev']
        maxdim   = ['lat','lon','lev']
    else:
        maxdim_n = ['ncol']
        maxdim = ['lat','lon']
    
    if(native is not None):
        print('using native grid for global diff')
        dn = xr.open_dataset(native)
        dq = (dn[varj] - dn[var])
        dqmax = np.abs(dq).max(maxdim_n) * 1e9
    else:
        print('using remapped data for global diff')
        dq = (d[varj] - d[var]) 
        dqmax = np.abs(dq).max(maxdim) * 1e9
    
    if(not varIs2d):
        d1 = d[varj][-1].mean('lon') * 1e9
        d2 = d[var][-1].mean('lon') * 1e9
        dqm = (d2 - d1)
        d1 = d1.sel(levslice)
        d2 = d2.sel(levslice)
        dqm = dqm.sel(levslice)
        lev = d1['lev']
    else:
        d1 = d[varj][-1]
        d2 = d[var][-1]
        dqm = (d2 - d1)

    if(not varIs2d):
        x, y       = lat, lev
        xlab, ylab = 'lat', 'p[hPa]'
        plotter    = pltvert
        argsp      = {}
        argsq      = {'plot_zscale':False}
        argsqc     = {'plot_zscale':False, 'inverty':False}
    else:
        x, y       = lon, lat
        xlab, ylab = 'lon', 'lat'
        plotter    = plthor
        argsp      = {'projection':ccrs.Miller()}
        argsq      = {}
        argsqc     = {}


    print('-------- plotting diff')
    fig = plt.figure(figsize=(18, 6))
    ax1 = fig.add_subplot(131)
    ax2 = fig.add_subplot(132, **argsp)
    ax3 = fig.add_subplot(133, **argsp)
    ax1.plot(time, dqmax, '-', color=cdiff)
    ax1.set_ylabel('global max. gridpoint {} difference\n (|{} impl - {} impl|) [{}]'.format(
                    var, l2, l1, unit))
    ax1.set_xlabel('time [days]')
    ax1.grid()

    print('-------- plotting q')
    pltargs = {'levels':levels, 'cmap':cmap, 'zorder':0}
    cArgs = {'orientation':'horizontal', 'location':'top', 
             'label':'{} zonal-mean {}\n ({} impl) [{}]'.format(trange, var, l1, unit),
             'aspect':30,'format':'%d'}
    pltargs_c = {'levels':levels, 'colors':'k', 'linewidths':0.6, 'zorder':1}
    cArgs_c = {'fmt':'%d'}
    var_dict = [{'var':d1, 'plotType':'contourf', 'plotArgs':pltargs, 'colorArgs':cArgs}]
    var_dict_c = [{'var':d1, 'plotType':'contour', 'plotArgs':pltargs_c, 'colorArgs':cArgs_c}]
    cf = plotter(x, y, var_dict, ax=ax2, annotation='', **argsq)
    plotter(x, y, var_dict_c, ax=ax2, annotation='', **argsqc)
    ax2.set_xlabel(xlab)
    ax2.set_ylabel(ylab)
    if(type(levels) != int):
        cf[0].set_ticks(levels)
    
    print('-------- plotting qdiff')
    pltargs = {'levels':levdiff, 'cmap':cmap, 'zorder':0}
    cArgs = {'orientation':'horizontal', 'location':'top', 
             'label':'{} zonal-mean {} difference\n ({} impl - {} impl) [{}]'.format(trange, var, l2, 
              l1, unit), 'aspect':30,'format':'%.1f'}
    pltargs_c = {'levels':levdiff, 'colors':'k', 'linewidths':0.6, 'zorder':1}
    cArgs_c = {'fmt':'%.1f'}
    var_dict = [{'var':dqm, 'plotType':'contourf', 'plotArgs':pltargs, 'colorArgs':cArgs}]
    var_dict_c = [{'var':dqm, 'plotType':'contour', 'plotArgs':pltargs_c, 'colorArgs':cArgs_c}]    
    cf = plotter(x, y, var_dict, ax=ax3, annotation='', **argsq)
    plotter(x, y, var_dict_c, ax=ax3, annotation='', **argsqc) 
    ax3.set_xlabel(xlab)
    ax3.set_ylabel('')
    if(type(levdiff) != int):
        cf[0].set_ticks(levdiff)

    if(var == 'ST80_25'):
        aut.insert_labelled_tick(ax2, 'y', 150, label=None)
        aut.insert_labelled_tick(ax3, 'y', 150, label=None)

    plt.tight_layout()
    plt.savefig('{}_1year.png'.format(var), dpi=250)


if __name__ == '__main__':

    remap = '/global/homes/j/jhollo/repos/climate_analysis/CLDERA/E90ST80/cases/lowres_output/E3SM_F1850_ne4pg2_oQU480_L72_e90st80_bothimpl/run/E3SM_F1850_ne4pg2_oQU480_L72_e90st80_bothimpl.eam.h0.0001-01-01-00000.regrid.25x48_aave.nc'
    native = '/global/homes/j/jhollo/repos/climate_analysis/CLDERA/E90ST80/cases/lowres_output/E3SM_F1850_ne4pg2_oQU480_L72_e90st80_bothimpl/run/E3SM_F1850_ne4pg2_oQU480_L72_e90st80_bothimpl.eam.h0.0001-01-01-00000.nc'
    #trange='1year'
    trange='30day'
    compare_q('E90', remap, 'EAM', 'MOZART', trange, native)
    compare_q('ST80_25', remap, 'EAM', 'MOZART', trange, native)
    #compare_q('SFE90', remap, 'EAM', 'MOZART', trange, native, True, 'kg/m^2/s')
    plt.show()

