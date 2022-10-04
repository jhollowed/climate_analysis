import pdb
import glob
import textwrap
import numpy as np
import xarray as xr
import artist_utils as art
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import metpy.constants as const
from metpy.units import units as u
from climate_artist import vertical_slice as pltvert
from climate_artist import horizontal_slice as plthor


# =========================================================================


class pv_verify:
    '''
    Performs verification tests on our PV implementation in E3SM

    Parameters
    ----------
    histfile : str
        history file to use in verification tests
    name : str, oiptional
        string to use in plot titles etc.
    save_dest : str, optional
        destination for output figures and data. Defaults to None, in 
        which case figures are rendered to the display only and not 
        written to a file
    show_fig : bool, optional
        whether or not to display figures. Defaults to False.
    
    Attributes
    ----------
    data : xarrat Dataset
        Data read from the hist file
    pv_phys : xarray DataArray
        PV on the physics grid
    pv_dyn : xarray DataArray
        PV on the dynamics grid
    
    Methods
    -------
    global_maxmin_dyn_phys():
        Plots the global max and min of the dynamics and physics PV 
        fields as a time series
    zonal_mean_dyn_phys():
        Plots a zonal mean of the fractional disagreement between the 
        (remapped) dynamics and physics PV
    '''
    
    def __init__(self, histfile, name='', save_dest=None, show_fig=False):  
        self.name     = name
        self.save_dest = save_dest
        self.show_fig = show_fig
        self.data_read = False
        self.data    = xr.open_dataset(histfile)
        self.time = self.data['ndcur']

    def read_data(self):
        print('reading PV, DYN_PV for {}...'.format(self.name))
        self.pv_phys = self.data['PV']
        self.pv_dyn  = self.data['DYN_PV']
        self.data_read = True


    def global_maxmin_dyn_phys(self, infig=None):
        '''
        Plots the fractional difference in the global max and min of the dynamics 
        and physics PV fields as a time series

        Parameters
        ----------
        infig : matplotlib Figure object, optional
            figure on which to plot. If provided, it is assumed that the axes to plot the
            max and min time series to are fig.axes[0]
            Defaults to None, in which case a new figure will be created.

        Returns
        -------
        fig : matplotlib Figure object
            if the input arg 'fig' was provided, this function returns that object rather
            than rendering images of or displaying any figures.
        '''
        
        # ---- create or load figure
        if infig is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)
        else:
            fig = infig
            ax = fig.axes[0]

        # ---- read data, compute min/max
        if(self.save_dest is not None):
            minmax_tmpfile = '{}/tmp_global_minmax_{}.npy'.format(self.save_dest, name)
            try:
                [max_pv_phys, min_pv_phys, max_pv_dyn, min_pv_dyn] = \
                                           np.load(minmax_tmpfile, allow_pickle=True)
                print('read {} global min/max PV from {}...'.format(self.name, minmax_tmpfile))
            except FileNotFoundError:
                print('minmax tmp file not found; computing global'\
                      'max.min for {}...'.format(self.name))
                if not self.data_read: self.read_data()
                max_pv_phys = self.pv_phys.max(['ncol', 'lev'])
                min_pv_phys = self.pv_phys.min(['ncol', 'lev'])
                max_pv_dyn = self.pv_phys.max(['ncol', 'lev'])
                min_pv_dyn = self.pv_phys.min(['ncol', 'lev'])
                print('writing global minmax to {}...'.format(minmax_tmpfile))
                outdat = np.array([max_pv_phys, min_pv_phys, max_pv_dyn, min_pv_dyn], dtype=object)
                np.save(minmax_tmpfile, outdat, allow_pickle=True)
        else:
            if not self.data_read: self.read_data()
            # scale from m2 K/kg/s to PVU
            print('save_dest not specified; computing global max.min for {}...'.format(self.name))
            max_pv_phys = self.pv_phys.max(['ncol', 'lev'])
            min_pv_phys = self.pv_phys.min(['ncol', 'lev'])
            max_pv_dyn = self.pv_dyn.max(['ncol', 'lev'])
            min_pv_dyn = self.pv_dyn.min(['ncol', 'lev'])
        
        pdiff_max_pv = (max_pv_phys - max_pv_dyn) / max_pv_dyn
        pdiff_min_pv = (min_pv_phys - min_pv_dyn) / min_pv_dyn
        pdb.set_trace()
        
        # ---- plot
        ax.plot(self.time, pdiff_max_pv, 'r', lw=2, label='global max')
        ax.plot(self.time, pdiff_min_pv, 'r--', lw=2, label='global min')

        ax.set_xlabel('time [days]')
        ax.set_ylabel('global PV extrema fractional difference [%]')
        ax.set_title(self.name)
        ax.grid()
        ax.legend()
        plt.tight_layout()

        if infig is None:
            if self.save_dest is not None:
                plt.savefig('{}/{}_global_minmax.png'.format(self.save_dest, name), dpi=300)
            if self.show_fig:
                plt.show()
        else:
            return fig
        

    def zonal_mean_dyn_phys(self, remapped_histfiles, t, t2=None):
        '''
        Plots a zonal mean of the fractional disagreement between the
        (remapped) dynamics and physics PV

        Parameters
        ----------
        remapped_histfile : len 2 str list, optional
            lat-lon remapped history files. It is expected that the PV on the physics and 
            dynamics grid were reampped separately, and live on two separate netcdf files. 
            The first entry in this list should be the remapped data from the physics grid
            (including 'PV'), and the second should be the remapped data from the dynamics 
            grid (including 'DYN_PV'). It is assumed that the remapped lat,lon coordinates 
            are the same for each dataset.
        t : float
            time at which to take the zonal mean, in days
        t2 : float, optional
            if this is provided, then plot the time-mean zonal-mean over the window [t, t2], in days
        '''
        
        # ---- configure time slice
        tsel = t
        tstr = 'day {}'.format(t)
        if t2 is not None: 
            tsel = slice(t, t2)
            tstr = 'time mean days {}-{}'.format(t, t2)
        tsel = {'time':tsel}
        
        # ---- label figure
        label = '{} physics-dynamics grid fractional zonal-mean PV difference, {} [%]'.format(
                                                                               self.name, tstr)
        # ---- read remapped data, write pdiff to tmp file, read if exists
        remap_phys_data = xr.open_dataset(remapped_histfiles[0])
        remap_dyn_data  = xr.open_dataset(remapped_histfiles[1])
        remap_lat       = remap_dyn_data['lat']
        remap_lev       = remap_dyn_data['lev']
        x, y = remap_lat, remap_lev

        remap_phys_data = remap_phys_data.assign_coords(time=self.time)
        remap_dyn_data = remap_dyn_data.assign_coords(time=self.time)
       
        if(self.save_dest is not None):
            pdiff_tmpfile = '{}/tmp_pv_pdiff_data_{}.nc'.format(self.save_dest, name)
            try:
                pdiff = xr.open_dataset(pdiff_tmpfile)
                print('read {} from {}...'.format(label, pdiff_tmpfile))
            except FileNotFoundError:
                print('pdiff tmp file not found; computing {}...'.format(label))
                remap_pv_phys   = remap_phys_data['PV']
                remap_pv_dyn    = remap_dyn_data['DYN_PV']
                zm_pv_phys = remap_pv_phys.sel(tsel).mean(['time', 'lon'])
                zm_pv_dyn = remap_pv_dyn.sel(tsel).mean(['time', 'lon'])
                pdiff = (zm_pv_phys - zm_pv_dyn) / zm_pv_dyn
                print('writing {} to {}...'.format(label, pdiff_tmpfile))
                pdiff.to_netcdf(pdiff_tmpfile)
        else:
            print('no save_dest specified; computing {}...'.format(label))
            remap_pv_phys   = remap_phys_data['PV']
            remap_pv_dyn    = remap_dyn_data['DYN_PV']
            zm_pv_phys = remap_pv_phys.sel(tsel).mean(['time', 'lon'])
            zm_pv_dyn = remap_pv_dyn.sel(tsel).mean(['time', 'lon'])
            pdiff = (zm_pv_phys - zm_pv_dyn) / zm_pv_dyn

        # ---- plot 
        levels = 10
        cmap = plt.cm.rainbow
        label = textwrap.fill(label, int(len(label)/2))
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        pltargs = {'levels':levels, 'cmap':cmap, 'zorder':0}
        cArgs = {'orientation':'horizontal', 'location':'top', 'label':label, 
                 'aspect':30,'format':'%d'}
        pltargs_c = {'levels':levels, 'colors':'k', 'linewidths':0.6, 'zorder':1}
        cArgs_c = {'fmt':'%d'}
        var_dict = [{'var':pdiff, 'plotType':'contourf', 'plotArgs':pltargs, 'colorArgs':cArgs}]
        var_dict_c = [{'var':pdiff, 'plotType':'contour', 'plotArgs':pltargs_c, 'colorArgs':cArgs_c}]
        cf = pltvert(x, y, var_dict, ax=ax2, annotation='')
        pltvert(x, y, var_dict_c, ax=ax2, annotation='', invert_y=False, plot_zscale=False)
        ax2.set_xlabel('lat [hPa]')
        ax2.set_ylabel('p [hPa]')
        if(type(levels) != int):
            cf[0].set_ticks(levels)
        plt.tight_layout()
        
        if self.save_dest is not None:
            plt.savefig('{}/{}_pv_pdiff_{}.png'.format(self.save_dest, name, tstr), dpi=300)
        if(self.show_fig):
            plt.show()
        
        return


# ==========================================================================


topdir = '/global/cscratch1/sd/jhollo/E3SM/E3SMv2_cases/pv_cases'
rundirs = glob.glob('{}/E3SM*/run/'.format(topdir))
for rundir in rundirs:
    histfile = glob.glob('{}/*eam.h0*00.nc'.format(rundir))[0]
    remapped_histfiles = sorted(glob.glob('{}/*eam.h0*regrid*.nc'.format(rundir)))[::-1]
    name = histfile.split('/')[-3]
    pvv = pv_verify(histfile, name, save_dest='.', show_fig=True)
    pvv.global_maxmin_dyn_phys()
    pvv.zonal_mean_dyn_phys(remapped_histfiles, 20, 30)
