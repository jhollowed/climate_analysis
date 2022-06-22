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

RES=500
pref_mid_norm = (np.logspace(np.log10(ptop_ref), np.log10(100000), RES)) / psurf_ref
pmid = pref_mid_norm * psurf_ref
clat = np.linspace(-90, 90, RES)*np.pi/180

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

fig = plt.figure()
ax2 = fig.add_subplot(111)

CS = ax2.contour(PHI, P, teq_max.T, colors='k', levels=np.arange(150, 320, 10))

rev = 1/(86400*3) # 1/3day
pih = np.pi/2


# ----- for RF colormap plot ----
cutoff = 1 # 1 hPa
num = pih*np.log(cutoff/P)
den = np.log(cutoff/np.min(P))
kr = rev * (np.sin(num/den))**2
for i in range(kr.shape[0]):
    for j in range(kr.shape[1]):
        if(P[i][j] > cutoff): kr[i][j] = 0

vmin, vmax = 0, np.max(kr)
dv = vmax*2
N = 256
hot_r = plt.cm.get_cmap('hot_r', N)
newcolors = hot_r(np.linspace(0, 0.8, N))
newmap = mpl.colors.ListedColormap(newcolors)

cc = ax2.pcolor(PHI, P, kr, cmap=newmap, zorder=0)
cbar = fig.colorbar(cc, label='RF coefficient [1/s]', orientation='vertical', ax=ax2)
cbar.ax.tick_params(labelsize=9)
#cax.xaxis.set_label_position('top')
#cax.xaxis.set_ticks_position('top')
#cax.set_xlim([0, np.max(kr)])

#ax2.plot([-90, 90], [3, 3], '--r', label='Williamson+97 model top')
#ax2.plot([-90, 90], [1, 1], '--b', label='stratopause')
#ax2.legend(loc='upper center', bbox_to_anchor=(0.5,1.15), ncol=1, frameon=False, fontsize=10)



ax2.set_ylim([min(pmid)/100, max(pmid)/100])
ax2.invert_yaxis()
ax2.set_yscale('log')
ax2.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '{:g}'.format(y)))
ax2.set_xlabel('φ  [deg]', fontsize=11)
ax2.clabel(CS, fmt='%.0f', manual=False, levels = [200, 220, 170])


plt.tight_layout()
#plt.show()
plt.savefig('WHS_new.png', dpi=300)

