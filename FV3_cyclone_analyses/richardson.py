import numpy as np
import matplotlib.pyplot as plt
import astropy.units as uu
import astropy.constants as const
import Nio
import Ngl
import subprocess
from math import radians
import glob
import pdb
import sys
from guppy import hpy

sys.path.append('/home/hollowed/repos/ncl_exports')
import wrappers

class richardson:
    def __init__(self, runs, names=None, cmap=plt.cm.gist_rainbow):
        self.files = [Nio.open_file(f, 'wr') for f in runs]
        if names is None:
            self.run_names = [s.split('/')[-1] for s in runs]
        else: 
            self.run_names = names
        self.colors = cmap(np.linspace(0.1, 0.9, len(self.files)))
        
        # add string giving file loaction to attributes of NioFile object 
        # (temporary fix; this doesn't seem to exist natively yet)
        for i in range(len(self.files)):
            self.files[i].location = runs[i]
   

    def grad(self, u, v, fname, status=True):
       
        try:
            dudv = np.load(fname)
        except FileNotFoundError:
            dim = np.shape(u)
            dudv = np.zeros(dim)
            for t in range(dim[0]):
                if(status): print('grad: {}/{}'.format(t, dim[0]))
                for x in range(dim[2]):
                    for y in range(dim[3]):
                        dudv[t,:,x,y] = np.gradient(u[t,:,x,y], v[t,:,x,y], edge_order=2)
            np.save(fname, dudv)
        return dudv


    def Rc(self, p):
        
        Rc = np.zeros(np.shape(p))
        
        #in Pa
        Rc[p < 40000] = 0.25
        Rc[p > 60000] = 1.0
        
        ramp_mask = np.logical_and(p>40000, p<60000)
        pr = p[ramp_mask]
        m, b = (0.75/20000), -1.25
        Rc[ramp_mask] = m*pr + b

        return Rc


    def richardson_z(self):

        for i in range(len(self.files)):
            if i==0: continue
            fig = plt.figure()
            ax1 = fig.add_subplot(111)
            
            print('reading variables for {}'.format(self.run_names[i]))
  
            f = self.files[i]
            pdb.set_trace()
            print('1')
            T = f.variables['T'][:]
            print('2')
            Q = f.variables['Q'][:]
            print('3')
            PS = f.variables['PS'][:]
            print('4')
            hyam = f.variables['hyam'][:]
            print('5')
            hybm = f.variables['hybm'][:]
            print('6')
            P0 = f.variables['P0'].get_value()
            print('7')
            Rd = 287.058
            cp = 1005
            g = 9.8065
            
            print('computing P for {}'.format(self.run_names[i]))
            
            # performs P = hyam*P0 + hybm*PS
            P = np.swapaxes([hybm[k] * PS for k in range(len(hybm))], 0, 1)
            for k in range(len(hybm)):
                P[:,k,:,:] += hyam[k]*P0
            
            print('computing N^2 for {}'.format(self.run_names[i]))
            
            # get virtual temp, pot. temp
            Tv = T*(1+0.61*Q)
            thetav = Tv*(P0/P)**(Rd/cp)
            f.variables['THETAV'] = thetav

            # interpolate pot temp. to height positions, approx gradient
            levels = np.linspace(100, 60000, 100)
            try:
                dtdz = np.load('tmp/dtdz_{}.npy'.format(i))
            except FileNotFoundError:
                thetav_z = wrappers.vertical_interp(f, 'THETAV', 'z', levels)
                dtdz = np.gradient(thetav_z, levels, edge_order=2, axis=0)
                np.write(dtdz, 'tmp/dtdz_{}.npy'.format(i))
            N2 = g/thetav * dtdz

            print('computing Ri for {}'.format(self.run_names[i]))
            
            # interpolate u,v to height positions, approx gradient
            try:
                vv = np.load('tmp/vv_{}.npy'.format(i))
            except FileNotFoundError:
                u = wrappers.vertical_interp(f, 'U', 'z', levels)
                dudz = np.gradient(u, levels, edge_order=2, axis=0)
                v = wrappers.vertical_interp(f, 'V', 'z', levels)
                dvdz = np.gradient(v, levels, edge_order=2, axis=0)
                
                vv = np.linalg.norm([dudz, dvdz], axis=0)**2
                np.write(vv, 'tmp/vv_{}.npy'.format(i))
            
            # take mag squared, compute Ri
            Ri = N2/vv

            pdb.set_trace()
    
    
    def richardson_p(self):

        h = hpy()

        for i in range(len(self.files)):
            if i==0: continue
            fig = plt.figure()
            ax1 = fig.add_subplot(111)
            f = self.files[i]
            
            print('localizing cyclone')
            # take rough localization mid-way through the simulation, draw (nn x nnl) degree box
            # (nn in number of grid points)
            PS = f.variables['PS'][:]
            clat, clon = np.unravel_index(np.argmin(PS[15]), np.shape(PS[15]))
            nn = 30
            lonf = 1.5
            nnl = int(nn*lonf)
           
            print('reading variables for {}'.format(self.run_names[i]))
            PS = PS[:,clat-nn:clat+nn, clon-nnl:clon+nnl]
            T = f.variables['T'][:,:,clat-nn:clat+nn, clon-nnl:clon+nnl]
            Q = f.variables['Q'][:,:,clat-nn:clat+nn, clon-nnl:clon+nnl]
            U = f.variables['U'][:,:,clat-nn:clat+nn, clon-nnl:clon+nnl]
            V = f.variables['V'][:,:,clat-nn:clat+nn, clon-nnl:clon+nnl]
            Z3 = f.variables['Z3'][:,:,clat-nn:clat+nn, clon-nnl:clon+nnl]
            hyam = f.variables['hyam'][:]
            hybm = f.variables['hybm'][:]
            P0 = f.variables['P0'].get_value()
            Rd = 287.058
            cp = 1005
            g = 9.8065
            eps = 0.61
 
            print('computing P') 
            # performs P = hyam*P0 + hybm*PS
            P = np.swapaxes([hybm[k] * PS for k in range(len(hybm))], 0, 1)
            for k in range(len(hybm)):
                P[:,k,:,:] += hyam[k]*P0
            
            print('computing pressure vertical gradient') 
            dpdz = self.grad(P, Z3, 'tmp/dpdz_{}.npy'.format(i))
            RHO = dpdz/(-g)
           
            print('computing Ri') 
            # get virtual temp, pot. temp, thetav and velocity gradients
            w = Q/(1-Q)
            Tv = T*(1+eps*w)
            thetav = Tv*(P0/P)**(Rd/cp)
            dtdp = self.grad(thetav, P, 'tmp/dtdp_{}.npy'.format(i))
            dudp = self.grad(U, P, 'tmp/dudp_{}.npy'.format(i))
            dvdp = self.grad(V, P, 'tmp/dvdp_{}.npy'.format(i))
            vv = np.linalg.norm([dudp, dvdp], axis=0)**2 

            # get richardson number, critical richardson number
            Ri = -dtdp/(thetav*RHO * vv)
            Rc = self.Rc(P)
            mixing = Ri < Rc

            pdb.set_trace()


if __name__ == '__main__':

    ff = sorted(glob.glob('/scratch/cjablono_root/cjablono1/hollowed/cesm2.2/'\
                          'cyclone_tests/clones/*_30/run/*h0*regrid*.nc'))
    names = ['FV3: C192, sg_adj=1800, n_sponge=30', 'FV3: C192, sg_adj=3600, n_sponge=30']
    aa = cyclone_analysis(ff, names)
    aa.richardson_p()


