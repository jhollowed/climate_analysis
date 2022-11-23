import pdb
import glob
import textwrap
import warnings
import numpy as np
import xarray as xr
import artist_utils as art
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import metpy.constants as const
from metpy.units import units as u
import matplotlib.colors as colors
from climate_artist import vertical_slice as pltvert
from climate_artist import horizontal_slice as plthorz


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
        self.name      = name
        self.save_dest = save_dest
        self.show_fig  = show_fig
        self.data_read = False

        self.data   = xr.open_dataset(histfile)
        self.time   = self.data['ndcur']
        self.ncol   = self.data['ncol']
        self.lev    = self.data['lev']
        self.lat    = self.data['lat']
        self.lon    = self.data['lon']
        try:
            self.ncol_d = self.data['ncol_d']
            self.lat_d  = self.data['lat_d']
            self.lon_d  = self.data['lon_d']
        except KeyError:
            warnings.warn('no lat_d, lon_d, or ncol_d found in dataset; '\
                          'assuming phys and dyn grids are the same.')
            self.ncol_d = self.ncol
            self.lat_d = self.lat
            self.lon_d = self.lon

        self.data = self.data.assign_coords(time=self.time)

    def read_data(self):
        print('reading PV, DYN_PV for {}...'.format(self.name))
        self.pv_phys = self.data['PV']     * 1e6  # scale to PVU
        self.pv_dyn  = self.data['DYN_PV'] * 1e6  # scale to PVU
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
       
        plt.ion()
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
                max_pv_dyn = self.pv_dyn.max(['ncol_d', 'lev'])
                min_pv_dyn = self.pv_dyn.min(['ncol_d', 'lev'])
                print('writing global minmax to {}...'.format(minmax_tmpfile))
                outdat = np.array([max_pv_phys, min_pv_phys, max_pv_dyn, min_pv_dyn], dtype=object)
                np.save(minmax_tmpfile, outdat, allow_pickle=True)
        else:
            if not self.data_read: self.read_data()
            # scale from m2 K/kg/s to PVU
            print('save_dest not specified; computing global max.min for {}...'.format(self.name))
            max_pv_phys = self.pv_phys.max(['ncol', 'lev'])
            min_pv_phys = self.pv_phys.min(['ncol', 'lev'])
            max_pv_dyn = self.pv_dyn.max(['ncol_d', 'lev'])
            min_pv_dyn = self.pv_dyn.min(['ncol_d', 'lev'])
        
        #pdiff_max_pv = 100 * (max_pv_phys - max_pv_dyn) / max_pv_dyn
        #pdiff_min_pv = 100 * (min_pv_phys - min_pv_dyn) / min_pv_dyn
        # scale difference to PVU
        pdiff_max_pv = (max_pv_phys - max_pv_dyn) * 1e-6
        pdiff_min_pv = (min_pv_phys - min_pv_dyn) * 1e-6

        # ---- plot
        ax.plot(self.time, pdiff_max_pv, 'r', lw=2, label='global max')
        ax.plot(self.time, pdiff_min_pv, 'r--', lw=2, label='global min')

        ax.set_xlabel('time [days]')
        ax.set_ylabel('global PV extrema difference [PVU]')
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
        
    
    def levels_maxmin_dyn_phys(self):
        '''
        Plots the fractional difference in the horizontal global max and min of the dynamics 
        and physics PV fields as a time series per vertical level
        '''
       
        fig = plt.figure(figsize=(10,5))
        ax = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)

        # ---- read data, compute min/max
        if(self.save_dest is not None):
            minmax_tmpfile = '{}/tmp_levels_minmax_{}.npy'.format(self.save_dest, name)
            try:
                [max_pv_phys, min_pv_phys, max_pv_dyn, min_pv_dyn] = np.load(minmax_tmpfile)
                print('read {} levels level min/max PV from {}...'.format(self.name, minmax_tmpfile))
            except FileNotFoundError:
                print('minmax tmp file not found; computing global'\
                      'level max/min for {}...'.format(self.name))
                if not self.data_read: self.read_data()
                max_pv_phys = self.pv_phys.max('ncol')
                min_pv_phys = self.pv_phys.min('ncol')
                max_pv_dyn = self.pv_dyn.max('ncol_d')
                min_pv_dyn = self.pv_dyn.min('ncol_d')
                print('writing global minmax to {}...'.format(minmax_tmpfile))
                outdat = np.array([max_pv_phys, min_pv_phys, max_pv_dyn, min_pv_dyn])
                np.save(minmax_tmpfile, outdat)
        else:
            if not self.data_read: self.read_data()
            # scale from m2 K/kg/s to PVU
            print('save_dest not specified; computing global level max/min '\
                  'for {}...'.format(self.name))
            max_pv_phys = self.pv_phys.max('ncol')
            min_pv_phys = self.pv_phys.min('ncol')
            max_pv_dyn = self.pv_dyn.max('ncol_d')
            min_pv_dyn = self.pv_dyn.min('ncol_d')
        
        #pdiff_max_pv = 100 * (max_pv_phys - max_pv_dyn) / max_pv_dyn
        #pdiff_min_pv = 100 * (min_pv_phys - min_pv_dyn) / min_pv_dyn
        # scale difference to PVU
        pdiff_max_pv = (max_pv_phys - max_pv_dyn) * 1e-6
        pdiff_min_pv = (min_pv_phys - min_pv_dyn) * 1e-6
        
        # ---- plot
        Y, X = np.meshgrid(self.lev, self.time)
       
        c = ax.pcolormesh(X, Y, pdiff_max_pv, cmap=plt.cm.viridis, shading='nearest')
        c2 = ax2.pcolormesh(X, Y, pdiff_min_pv, cmap=plt.cm.viridis, shading='nearest')
        
        #cb = plt.colorbar(c, ax=ax, label='horizontal PV maximum fractional difference [%]')
        #cb2 = plt.colorbar(c2, ax=ax2, label='horizontal PV minimum fractional difference [%]')
        cb = plt.colorbar(c, ax=ax, label='horizontal PV maximum difference [PVU]')
        cb2 = plt.colorbar(c2, ax=ax2, label='horizontal PV minimum difference [PVU]')

        ax.set_yscale('log')
        ax.invert_yaxis()
        ax2.set_yscale('log')
        ax2.invert_yaxis()

        ax.set_xlabel('time [days]')
        ax.set_ylabel('lev [hPa]')
        ax2.set_xlabel('time [days]')
        ax2.set_ylabel('lev [hPa]')
        ax.grid()
        ax2.grid()
        fig.suptitle(self.name)
        plt.tight_layout()

        if self.save_dest is not None:
            plt.savefig('{}/{}_levels_minmax.png'.format(self.save_dest, name), dpi=300)
        if self.show_fig:
            plt.show() 


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
        # ---- read remapped data, write pdiff to tmp file, read if exists
        remap_phys_data = xr.open_dataset(remapped_histfiles[0])
        remap_dyn_data  = xr.open_dataset(remapped_histfiles[1])
        remap_lat       = remap_dyn_data['lat']
        remap_lev       = remap_dyn_data['lev']
        x, y            = remap_lat, remap_lev
        remap_phys_data = remap_phys_data.assign_coords(time=self.time)
        remap_dyn_data  = remap_dyn_data.assign_coords(time=self.time)
       
        if(self.save_dest is not None or overwrite):
            
            zm_pv_phys_tmpfile = '{}/tmp_zm_pv_phys_data_{}.npy'.format(self.save_dest, name)
            try:
                zm_pv_phys = np.load(zm_pv_phys_tmpfile)
                print('read from {}...'.format(zm_pv_phys_tmpfile))
            except FileNotFoundError:
                print('zm_pv_phys tmp file not found; computing...')
                remap_pv_phys   = remap_phys_data['PV'] * 1e-6  # scale to pvu
                zm_pv_phys      = remap_pv_phys.sel(tsel).mean(['time', 'lon'])
                print('writing to {}...'.format(zm_pv_phys_tmpfile))
                np.save(zm_pv_phys_tmpfile, zm_pv_phys)
            
            zm_pv_dyn_tmpfile = '{}/tmp_zm_pv_dyn_data_{}.npy'.format(self.save_dest, name)
            try:
                zm_pv_dyn = np.load(zm_pv_dyn_tmpfile)
                print('read from {}...'.format(zm_pv_dyn_tmpfile))
            except FileNotFoundError:
                print('zm_pv_dyn tmp file not found; computing...')
                remap_pv_dyn   = remap_dyn_data['DYN_PV'] * 1e-6  # scale to pvu
                zm_pv_dyn      = remap_pv_dyn.sel(tsel).mean(['time', 'lon'])
                print('writing to {}...'.format(zm_pv_dyn_tmpfile))
                np.save(zm_pv_dyn_tmpfile, zm_pv_dyn)
            
            pdiff_tmpfile = '{}/tmp_zm_pv_pdiff_data_{}.npy'.format(self.save_dest, name)
            try:
                pdiff      = np.load(pdiff_tmpfile)
                print('read from {}...'.format(pdiff_tmpfile))
            except FileNotFoundError:
                print('pdiff tmp file not found; computing...')
                pdiff = (zm_pv_phys - zm_pv_dyn)
                print('writing to {}...'.format(pdiff_tmpfile))
                np.save(pdiff_tmpfile, pdiff)
        else:
            print('no save_dest specified; computing {}...'.format(label))
            remap_pv_phys = remap_phys_data['PV'] * 1e-6  # scale to pvu
            remap_pv_dyn  = remap_dyn_data['DYN_PV'] * 1e-6  # scale to pvu
            zm_pv_phys    = remap_pv_phys.sel(tsel).mean(['time', 'lon'])
            zm_pv_dyn     = remap_pv_dyn.sel(tsel).mean(['time', 'lon'])
            pdiff         = (zm_pv_phys - zm_pv_dyn)

        pdb.set_trace()
        # ---- plot 
        #levels = np.linspace(-10, 10, 14)
        levels = 12
        levelsd = np.array([0, 0.1, 0.5, 0.75, 1, 1.5, 2, 3, 4, 5, 10, 20, 40, 60, 80, 100, 200])
        levelsd = np.hstack([-levelsd, levelsd])
        levelsp = np.array([0, 0.1, 0.5, 0.75, 1, 1.5, 2, 3, 4, 5, 10, 20, 40, 60, 80, 100, 200])
        levelsp = np.hstack([-levelsp, levelsp])
        cmap = plt.cm.rainbow
        
        fig = plt.figure(figsize=(15, 5))
        ax = fig.add_subplot(133)
        ax1 = fig.add_subplot(131)
        ax2 = fig.add_subplot(132)
       
        # ---- pv dyn
        pltargs = {'levels':levelsp, 'cmap':cmap, 'zorder':0}
        cArgs = {'orientation':'horizontal', 'location':'top', 'aspect':30,'format':'%.2f',
                 'label':'dyn grid zonal-mean PV, {}'.format(tstr)}
        var_dict =[{'var':zm_pv_dyn*1e-6,'plotType':'contourf','plotArgs':pltargs,'colorArgs':cArgs}]
        cf = pltvert(x, y, var_dict, ax=ax1, annotation='', xlabel='lat [hPa]', ylabel='p [hPa]')
        
        # ---- pv phys
        pltargs = {'levels':levelsd, 'cmap':cmap, 'zorder':0}
        cArgs = {'orientation':'horizontal', 'location':'top', 'aspect':30,'format':'%.2f',
                 'label':'phys grid zonal-mean PV, {}'.format(tstr)}
        var_dict=[{'var':zm_pv_phys*1e-6,'plotType':'contourf','plotArgs':pltargs,'colorArgs':cArgs}]
        cf = pltvert(x, y, var_dict, ax=ax2, annotation='', xlabel='lat [hPa]', ylabel='p [hPa]')
        
        # ---- pv diff
        label = '{} physics-dynamics grid zonal-mean PV difference, {} [PVU]'.format(
                                                                      self.name, tstr)
        label = textwrap.fill(label, int(len(label)/2))
        pltargs = {'levels':levels, 'cmap':cmap, 'zorder':0}
        cArgs = {'orientation':'horizontal', 'location':'top', 'label':label, 
                 'aspect':30,'format':'%.2f'}
        var_dict = [{'var':pdiff, 'plotType':'contourf', 'plotArgs':pltargs, 'colorArgs':cArgs}]
        cf = pltvert(x, y, var_dict, ax=ax, annotation='', xlabel='lat [hPa]', ylabel='p [hPa]')
        
        fig.suptitle(self.name)
        plt.tight_layout()
        
        if self.save_dest is not None:
            plt.savefig('{}/{}_pv_pdiff_{}.png'.format(self.save_dest, name, tstr), dpi=300)
        if(self.show_fig):
            plt.show()        
        return
    

    def hslice_dyn_phys(self, t, t2=None, plevel=100):
        '''
        Plots the disagreement between the dynamics and physics PV at a certain level

        Parameters
        ----------
        t : float
            time at which to take the zonal mean, in days
        t2 : float, optional
            if this is provided, then plot the time-mean zonal-mean over the window [t, t2], in days
        plevel : float, optional
            pressure level to display in hPa (will display nearest value to plevel). Default is 100 hPa
        '''
        
        # ---- configure time slice
        tsel = t
        tstr = 'day {}'.format(t)
        tmethod = 'nearest'
        if t2 is not None: 
            tsel = slice(t, t2)
            tstr = 'time mean days {}-{}'.format(t, t2)
            tmethod = None
        tsel = {'time':tsel}
       
        # ---- read and/or compute data arrays
        if(self.save_dest is not None):
            
            pvp_tmpfile   = '{}/tmp_horz_pvp_lev{}_{}.nc'.format(self.save_dest, plevel, self.name)
            try:
                pvp = xr.open_dataset(pvp_tmpfile)['PV']
                print('read from {}...'.format(pvp_tmpfile))
            except FileNotFoundError:
                print('pvp lev={} tmp file not found; computing...'.format(plevel))
                if(not self.data_read):
                    self.read_data()
                pvp = self.pv_phys.sel({'lev':plevel}, method='nearest')
                pvp = pvp.sel(tsel, method=tmethod)
                if(t2 is not None): pvp = pvp.mean('time')
                print('writing to {}...'.format(pvp_tmpfile))
                pvp.to_netcdf(pvp_tmpfile) 

            pvd_tmpfile   = '{}/tmp_horz_pvd_lev{}_{}.nc'.format(self.save_dest, plevel, self.name)
            try:
                pvd = xr.open_dataset(pvd_tmpfile)['DYN_PV']
                print('read from {}...'.format(pvd_tmpfile))
            except FileNotFoundError:
                print('pvd lev={} tmp file not found; computing...'.format(plevel))
                if(not self.data_read):
                    self.read_data()
                pvd = self.pv_dyn.sel({'lev':plevel}, method='nearest')
                pvd = pvd.sel(tsel, method=tmethod)
                if(t2 is not None): pvd = pvd.mean('time')
                print('writing to {}...'.format(pvd_tmpfile))
                pvd.to_netcdf(pvd_tmpfile)
        else:
            print('no save_dest specified; computing {}...'.format(label))
            if(not self.data_read()):
                seldf.read_data()
            pvp = self.pv_phys.sel({'lev':plevel}, method='nearest')
            pvp = pvp.sel({'time':tsel}, method=tmethod)
            pvd = self.pv_dyn.sel({'lev':plevel}, method='nearest')
            pvd = pvd.sel({'time':tsel}, method=tmethod)
            if(t2 is not None): 
                pvd = pvd.mean('time')
                pvp = pvp.mean('time')
        
        # ---- plot
        xp = self.lon
        yp = self.lat
        xd = self.lon_d
        yd = self.lat_d
        
        levelsdiff = 12
        levelsd = np.array([1, 3, 10, 30, 100, 300, 1000])
        levelsd = np.hstack([-levelsd[::-1], [0], levelsd])
        levelsp = np.array([1, 3, 10, 30, 100, 300, 1000])
        levelsp = np.hstack([-levelsp[::-1], [0], levelsp])
        if(plevel==950):
            print('===== USING 950 LEVELS')
            levelsdiff = 12
            levelsd = np.array([1, 2, 3, 4, 5, 6, 7])
            levelsd = np.hstack([-levelsd[::-1], [0], levelsd])
            levelsp = np.array([0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4])
            levelsp = np.hstack([-levelsp[::-1], [0], levelsp])
        
        cmap = plt.cm.rainbow
        norm = colors.SymLogNorm(linthresh=2, linscale=1)
        
        fig = plt.figure(figsize=(10, 5))
        #ax1 = fig.add_subplot(121, projection=ccrs.PlateCarree())
        #ax2 = fig.add_subplot(122, projection=ccrs.PlateCarree())
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)

        # ---- pv dyn
        pltargs = {'levels':levelsd, 'cmap':cmap, 'zorder':0, 'norm':norm}
        cArgs = {'orientation':'horizontal', 'location':'top', 'aspect':30,'format':'%.1f',
                 'label':'dyn grid PV [PVU]'}
        var_dict =[{'var':pvd,'plotType':'tricontourf','plotArgs':pltargs,'colorArgs':cArgs}]
        cf = plthorz(xd, yd, var_dict, ax=ax1, include_contour_labels=False, 
                     xlabel='lon [deg]', ylabel='lat [deg]')
        ax1.set_ylim([-90, 90])
        ax2.set_xlim([0, 360])
        cf[0].ax.tick_params(rotation=90)
        
        # ---- pv phys
        #pltargs = {'levels':levelsp, 'cmap':cmap, 'norm':norm}
        pltargs = {'levels':levelsp, 'cmap':cmap}
        cArgs = {'orientation':'horizontal', 'location':'top', 'aspect':30,'format':'%.1f',
                 'label':'phys grid PV [PVU]'}
        var_dict =[{'var':pvp,'plotType':'tricontourf','plotArgs':pltargs,'colorArgs':cArgs}]
        cf = plthorz(xp, yp, var_dict, ax=ax2, include_contour_labels=False,
                     xlabel='lon [deg]', ylabel='lat [deg]')
        ax2.set_ylim([-90, 90])
        ax2.set_xlim([0, 360])
        cf[0].ax.tick_params(rotation=90)
        
        # ---- pv diff
        fig.suptitle('{}, lev={:.2f}, {}'.format(self.name, float(pvp['lev']), tstr))
        fig.tight_layout()
        
        if self.save_dest is not None:
            plt.savefig('{}/{}_pv_hor_lev{}_{}.png'.format(
                        self.save_dest, self.name, plevel, tstr.replace(' ', '_')), 
                        dpi=300)
        if(self.show_fig):
            plt.show()        
        return


# ==========================================================================


topdir = '/global/cscratch1/sd/jhollo/E3SM/E3SMv2_cases/pv_cases'
#rundirs = glob.glob('{}/E3SM*/run/'.format(topdir))
#rundirs = glob.glob('{}/E3SM*F1850*/run/'.format(topdir))
rundirs = glob.glob('{}/E3SM*ne16pg2*FIDEAL*/run/'.format(topdir))
print(rundirs)
for rundir in rundirs:
    histfile = glob.glob('{}/*eam.h0*00.nc'.format(rundir))[0]
    remapped_histfiles = sorted(glob.glob('{}/*eam.h0*regrid*bilinear.nc'.format(rundir)))[::-1]
    name = histfile.split('/')[-3].split('E3SM_')[-1].split('_PV')[0].replace('_L72_', '_')
    pvv = pv_verify(histfile, name, save_dest='.', show_fig=True)
    #pvv.global_maxmin_dyn_phys()
    #pvv.levels_maxmin_dyn_phys()
    #pvv.zonal_mean_dyn_phys(remapped_histfiles, 20, 30)
    #pvv.hslice_dyn_phys(20, 30)
    #pvv.hslice_dyn_phys(30, plevel=100)
    pvv.hslice_dyn_phys(30, plevel=950)
plt.show()
