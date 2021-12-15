import Nio
import Ngl
import pdb
import glob
import subprocess
from util import * 
import numpy as np
import matplotlib.pyplot as plt


class cyclone_analysis:
    '''
    Creates an object which performs and visualizes analysis on tropical cyclone 
    CESM run outputs

    Parameters
    ----------
    runs : list of strings
        list of CESM history files containing the outputs of tropical cyclone 
        idealized test case runs. Full path required.
    names : list of strings
        names to give each run (for plotting, file naming). Defaults to None, 
        in which case the file names will be used (which are probably long and ugly)
    '''
    def __init__(self, runs, names=None, cmap=plt.cm.gist_rainbow):
        self._file_names = runs
        self.files = [Nio.open_file(f, 'r') for f in runs]
        if names is None:
            self.run_names = [s.split('/')[-1] for s in runs]
        else: 
            self.run_names = names


    def time_evolution(self, lev=28, cmap=plt.cm.gist_rainbow, ls=None):
        '''
        Computes and visualizes the time evolution for each cyclone in the input runs,
        including:
        - minimum surface pressure
        - maximum wind speed
        - radius of maximum wind
        This corresponds to Fig.5 of Reed+Jablonowski 2010.
        The compouted data is stored as 1D time series in .npy files in ./output
        The figures are rendered as pdfs to ./output

        Parameters
        ----------
        lev : int
            the level at which to measure the evolution of the max. wind speed and 
            radius of max. wind. Maybe to implement: specify by height instead of lev
        cmap : matplotlib colormap instance
            cmap to use to color code each run in the rendered figure
        ls : list of strings
            plot style formatting string per-run. Defaults to None, in which case all 
            plots will use the linestyle '-'
        '''

        # create the three panels
        fig = plt.figure(figsize=(12, 9))
        ax1 = fig.add_subplot(311)
        ax2 = fig.add_subplot(312)
        ax3 = fig.add_subplot(313)
        colors = cmap(np.linspace(0.1, 0.9, len(self.files)))
        l = (self.files[0].variables['lev'][:])[lev] # assume these are all the same...
        ax1.set_ylabel('Minimum Surface Pressure (hPa)', fontsize=11)  
        ax2.set_ylabel('Maximum wind speed (m/s) at {:.2f} hPa'.format(l), fontsize=11)
        ax3.set_xlabel('days', fontsize=12)
        ax3.set_ylabel('RMW (km)', fontsize=12)
        if ls is None:
            ls = ['-']*len(self.files)
 
        # loop over all input runs
        for i in range(len(self.files)):   
            print('working on file {}...'.format(self._file_names[i].split('/')[-1]))

            # ================== MIN P_S ====================
            PS = self.files[i].variables['PS'][:]
            ps_file = 'output/ps_{}.npy'.format(self.run_names[i])
            try:
                print('reading min PS for {}'.format(self.run_names[i]))
                minPS = np.load(ps_file)
            except FileNotFoundError:
                print('computing min PS for {}'.format(self.run_names[i]))
                minPS = np.min(PS, axis=1)
                np.save(ps_file, minPS)
            days = self.files[i].variables['time'][:]
            pp = ax1.plot(days, minPS/100, ls[i], color=colors[i], label=self.run_names[i])

            # ================== MAX WIND ====================
            U = self.files[i].variables['U'][:]
            V = self.files[i].variables['V'][:]
            vv = np.linalg.norm([U[:,lev,:], V[:,lev,:]], axis=0)
            vv_file = 'output/vv_{}.npy'.format(self.run_names[i])
            try:
                print('reading max wind speed for {}'.format(self.run_names[i]))
                max_vv = np.load(vv_file)
            except FileNotFoundError:
                print('computing max wind speed for {}'.format(self.run_names[i]))
                max_vv = np.max(vv, axis=1)
                np.save(vv_file, max_vv)
            ax2.plot(days, max_vv, ls[i], color=colors[i])
        
            # ================== RADIUS OF MAX WIND ====================
            rmw_file = 'output/rmw_{}.npy'.format(self.run_names[i])
            try:
                print('reading RMW for {}'.format(self.run_names[i]))
                rmw = np.load(rmw_file)
            except FileNotFoundError:
                print('computing RMW for {}'.format(self.run_names[i]))
                # get location of pressure minimum
                ps_idx = np.argmin(PS, axis=1)
                cyclone_lat = (self.files[i].variables['lat'][:])[ps_idx].compressed()
                cyclone_lon = (self.files[i].variables['lon'][:])[ps_idx].compressed()
 
                # get location of max. wind
                mw_idx = np.argmax(vv, axis=1)
                mw_lat = (self.files[i].variables['lat'][:])[mw_idx].compressed()
                mw_lon = (self.files[i].variables['lon'][:])[mw_idx].compressed()
 
                # get angular separation of points
                rmw = np.array([great_circle(cyclone_lat[i], cyclone_lon[i], mw_lat[i], mw_lon[i]) 
                                for i in range(len(days))])
                rmw /= 1000
                np.save(rmw_file, rmw)   
     
            #ax3.plot(days, rmw, ls[i], color=colors[i], label=self.run_names[i])
            ax3.plot(days, rmw, ls[i], color=colors[i])
       

        # done; save
        ax1.legend(bbox_to_anchor=(1.04, 1))
        plt.tight_layout()
        plt.savefig('output/cyclone_evolution.eps', format='eps')
    
    
    def vertical_structure(self, var, lev=28, cmap=plt.cm.viridis):
        '''
        Visualizes the vertical strcuture of a cyclone for the variable var, 
        '''
        return
    


# ============================================================================================
# ============================================================================================


         
if __name__ == '__main__':

    ff = sorted(glob.glob('/scratch/cjablono_root/cjablono1/hollowed/cesm2.2/'\
                          'cyclone_tests/clones/*/run/*h0*.nc'))
    #ff.append('/scratch/cjablono_root/cjablono1/hollowed/cesm2.2/cyclone_tests/'
    #          'root/cesm2.2.ne60.L30.RJ12/run/cesm2.2.ne60.L30.RJ12.cam.h0.0001-01-01-00000.nc')
    #names = ['FV3: C192, n_sponge=0', 'FV3: C192, sg_adj=1800, n_sponge=30', 
    #        'FV3: C192, sg_adj=3600, n_sponge=30', 'SE: ne30']
    names = [s.split('RJ12')[-1].split('.')[0][2:] for s in ff]
    aa = cyclone_analysis(ff, names)
    ls = ['-']*len(names)
    for i in range(len(names)):
        if(names[i][-1] == '6'): ls[i] = '--'
    aa.time_evolution(ls = ls, cmap = plt.cm.jet)


