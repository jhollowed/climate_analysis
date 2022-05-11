import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pdb
import matplotlib.ticker as ticker
from metpy.units import units as u
import metpy.constants as const

efoldf = 1
efolda = 40
efolds = 4
sigmab = 0.7
t00 = 200
onemsig = 1 - sigmab
kf = 1/(86400*efoldf)
ka = 1/(86400*efolda)
ks = 1/(86400*efolds)
cappa = 2/7
cpair = (const.Cp_d.to(u.J/u.K/u.kg)).m

efoldaa = 40
kaa = 1/(86400*efoldaa)
pi = 4*np.arctan(1)
pih = 2*np.arctan(1)
phi0   = 60*np.pi/180
dphi0  = 15*np.pi/180
a0     = 2.65/dphi0

psurf_ref = 1000 * 100
ptop_ref = 10
aeq    = 10000 / psurf_ref

apole  = 200 / psurf_ref
apole_max  = 10 / psurf_ref

rair = (const.Rd.to(u.J/u.K/u.kg)).m
gravit = const.g.m
lapsew = -3.345e-03
constw = rair*lapsew/gravit
lapsec =  2.00e-03
constc = rair*lapsec/gravit

pref_mid_norm = (np.logspace(np.log10(ptop_ref), np.log10(100000), 500)) / psurf_ref
pmid = pref_mid_norm * psurf_ref
clat = np.linspace(-90, 90, 500)*np.pi/180

coslat = np.zeros(len(clat))
sinsq = np.zeros(len(clat))
cossq = np.zeros(len(clat))
cossqsq = np.zeros(len(clat))
for i in range(len(clat)):
      coslat[i] = np.cos(clat[i])
      sinsq[i] = np.sin(clat[i])*np.sin(clat[i])
      cossq[i] = coslat[i]*coslat[i]
      cossqsq[i] = cossq[i]*cossq[i]

teq = np.zeros((len(clat), len(pmid)))
teq_max = np.zeros((len(clat), len(pmid)))

for k in range(len(pref_mid_norm)):
    for i in range(len(clat)):
        if(pref_mid_norm[k] > sigmab):
            kt = ka + (ks - ka)*cossqsq[i]*(pref_mid_norm[k] - sigmab)/onemsig
        else:
            kt = ka
        
        acoslat = abs(np.arccos(coslat[i]))
        trefc   = 315 - (60 * sinsq[i])

        trefa=(trefc - 10*cossq[i]*np.log10((pmid[k]/psurf_ref)))*(pmid[k]/psurf_ref)**cappa
        trefa   = max(t00,trefa)
        trefa_max = trefa

        if (pref_mid_norm[k] < aeq):
            trefa = t00*((pmid[k]/10000))**constc
            trefa_max = trefa

            p0strat = aeq - (aeq - apole)*0.5*(1 + np.tanh(a0*(acoslat - phi0)))
            if (pref_mid_norm[k] < p0strat):
                trefa = trefa + t00*( ((pref_mid_norm[k]/p0strat))**constw - 1)
            
            p0strat_max = aeq - (aeq - apole_max)*0.5*(1 + np.tanh(a0*(acoslat - phi0)))
            if (pref_mid_norm[k] < p0strat_max):
                trefa_max = trefa_max + t00*( ((pref_mid_norm[k]/p0strat_max))**constw - 1)
        
        teq[i, k] = trefa
        teq_max[i, k] = trefa_max

PHI, P = np.meshgrid(clat, pref_mid_norm)
P = P*psurf_ref/100
PHI = PHI*180/np.pi

fig = plt.figure(figsize=(8,5))
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

CS = ax1.contour(PHI, P, teq.T, colors='k', levels=np.arange(150, 320, 10))
ax1.set_ylim([min(pmid)/100, max(pmid)/100])
ax1.invert_yaxis()
ax1.set_yscale('log')
ax1.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '{:g}'.format(y)))
ax1.set_ylabel('p  [hPa]', fontsize=11)
ax1.set_xlabel('φ  [deg]', fontsize=11)
ax1.clabel(CS, fmt='%.0f', manual=False, levels = [200, 220, 170])
ax1.plot([0,0], [10,10], label='rayleigh friction coefficient')
ax1.plot([-90, 90], [3, 3], '--r', label='Williamson+97 model top')
ax1.plot([-90, 90], [1, 1], '--b', label='stratopause')
ax1.legend(loc='upper center', bbox_to_anchor=(0.5,1.22), ncol=1, frameon=False, fontsize=10)

CS = ax2.contour(PHI, P, teq_max.T, colors='k', levels=np.arange(150, 320, 10))
ax2.set_ylim([min(pmid)/100, max(pmid)/100])
ax2.invert_yaxis()
ax2.set_yscale('log')
ax2.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '{:g}'.format(y)))
ax2.set_xlabel('φ  [deg]', fontsize=11)
ax2.clabel(CS, fmt='%.0f', manual=False, levels = [200, 220, 170])
ax2.plot([-90, 90], [1, 1], '--b')
ax2.yaxis.tick_right()


cutoff = 100/psurf_ref # 100 Pa in sigma
rev = 1/(86400*3) # 1/3day
pih = np.pi/2

num = pih*np.log(cutoff/pref_mid_norm)
den = np.log(cutoff/pref_mid_norm[0])
kr = rev * (np.sin(num/den))**2
for i in range(len(kr)):
    if(pref_mid_norm[i] > cutoff): kr[i] = 0

ax3 = ax2.twiny()
ax3.plot(kr, pmid/100, '', lw=2)
ax3.set_xlabel('RF coefficient [1/s]', fontsize=10)


plt.tight_layout()
plt.savefig('WHS.png', dpi=300)

