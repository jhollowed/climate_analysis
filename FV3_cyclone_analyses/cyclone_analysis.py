import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
import astropy.constants as const
import Nio
import Ngl
import subprocess
from math import radians
import glob
import pdb

class cyclone_analysis:
    def __init__(self, runs, names=None, cmap=plt.cm.gist_rainbow):
        self.files = [Nio.open_file(f, 'r') for f in runs]
        if names is None:
            self.run_names = [s.split('/')[-1] for s in runs]
        else: 
            self.run_names = names
        self.colors = cmap(np.linspace(0.1, 0.9, len(self.files)))
 
    def rl_fig_5(self, lev=28):

        fig = plt.figure(figsize=(5, 9))
        ax1 = fig.add_subplot(311)
        ax2 = fig.add_subplot(312)
        ax3 = fig.add_subplot(313)
 
        for i in range(len(self.files)):   
           
           PS = self.files[i].variables['PS'][:]
           ps_file = 'ps_{}.npy'.format(self.run_names[i])
           try:
               print('reading min PS for {}'.format(self.run_names[i]))
               minPS = np.load(ps_file)
           except FileNotFoundError:
               print('computing min PS for {}'.format(self.run_names[i]))
               minPS = np.min(PS, axis=1)
               np.save(ps_file, minPS)
           days = self.files[i].variables['time'][:]
           ax1.plot(days, minPS/100, '-', color=self.colors[i], label=self.run_names[i])
           ax1.set_ylabel('Minimum Surface Pressure (hPa)', fontsize=11)
  

           U = self.files[i].variables['U'][:]
           V = self.files[i].variables['V'][:]
           vv = np.linalg.norm([U[:,lev,:], V[:,lev,:]], axis=0)
           vv_file = 'vv_{}.npy'.format(self.run_names[i])
           try:
               print('reading max wind speed for {}'.format(self.run_names[i]))
               max_vv = np.load(vv_file)
           except FileNotFoundError:
               print('computing max wind speed for {}'.format(self.run_names[i]))
               max_vv = np.max(vv, axis=1)
               np.save(vv_file, max_vv)
           l = (self.files[i].variables['lev'][:])[lev]
           ax2.plot(days, max_vv, '-', color=self.colors[i], label=self.run_names[i])
           ax2.set_ylabel('Maximum wind speed (m/s) at {:.2f} hPa'.format(l), fontsize=11)
        
           rmw_file = 'rmw_{}.npy'.format(self.run_names[i])
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

           ax3.plot(days, rmw, '-', color=self.colors[i], label=self.run_names[i])
           ax3.set_xlabel('days', fontsize=12)
           ax3.set_ylabel('RMW (km)', fontsize=12)
           ax3.legend()
        
        plt.tight_layout()
        plt.savefig('cyclone.eps', format='eps')
    

def great_circle(lat1, lon1, lat2, lon2):
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    return const.R_earth.to(u.m).value * \
           np.arccos((np.sin(lat1) * np.sin(lat2) + np.cos(lat1) * np.cos(lat2) * np.cos(lon1-lon2)))
            
         
if __name__ == '__main__':

    print('CHECK: {}'.format(great_circle(90, 0, -90, 0)/(const.R_earth.to(u.m).value*np.pi)))
    ff = sorted(glob.glob('/scratch/cjablono_root/cjablono1/hollowed/cesm2.2/'\
                          'cyclone_tests/clones/*/run/*h0*.nc'))
    ff.append('/scratch/cjablono_root/cjablono1/hollowed/cesm2.2/cyclone_tests/'
              'root/cesm2.2.ne60.L30.RJ12/run/cesm2.2.ne60.L30.RJ12.cam.h0.0001-01-01-00000.nc')
    names = ['FV3: C192, n_sponge=0', 'FV3: C192, sg_adj=1800, n_sponge=30', 
            'FV3: C192, sg_adj=3600, n_sponge=30', 'SE: ne30']
    aa = cyclone_analysis(ff, names)
    aa.rl_fig_5()


