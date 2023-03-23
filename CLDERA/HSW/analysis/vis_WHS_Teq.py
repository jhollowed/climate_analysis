import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pdb
import matplotlib.ticker as ticker
from metpy.units import units as u
import metpy.constants as const
from climate_artist import vertical_slice as pltvert


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

pref_mid_norm = (np.logspace(np.log10(ptop_ref), np.log10(100000), 1000)) / psurf_ref
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
teq_trunc = np.zeros((len(clat), len(pmid)))
teq_trunc1 = np.zeros((len(clat), len(pmid)))
teq_trunc2 = np.zeros((len(clat), len(pmid)))

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
        trefa_trunc = trefa
        trefa_trunc1 = trefa

        if (pref_mid_norm[k] < aeq):
            trefa = t00*((pmid[k]/10000))**constc
            trefa_max = trefa
            pp = (max(pref_mid_norm[k], apole))
            trefa_trunc = t00*((pp*psurf_ref/10000))**constc
            trefa_trunc2 = t00*((pp*psurf_ref/10000))**constc

            p0strat = aeq - (aeq - apole)*0.5*(1 + np.tanh(a0*(acoslat - phi0)))
            if (pref_mid_norm[k] < p0strat):
                trefa = trefa + t00*( ((pref_mid_norm[k]/p0strat))**constw - 1)
            
            p0strat_max = aeq - (aeq - apole_max)*0.5*(1 + np.tanh(a0*(acoslat - phi0)))
            if (pref_mid_norm[k] < p0strat_max):
                trefa_max = trefa_max + t00*( ((pref_mid_norm[k]/p0strat_max))**constw - 1)
            
            p0strat_trunc = aeq - (aeq - apole)*0.5*(1 + np.tanh(a0*(acoslat - phi0)))
            if (pref_mid_norm[k] < p0strat_trunc):
                trefa_trunc = trefa_trunc + t00*( (( pp/p0strat))**constw - 1)
                trefa_trunc2 = trefa_trunc2 + t00*( (( pp/p0strat))**constw - 1)
         
        teq[i, k] = trefa
        teq_max[i, k] = trefa_max
        teq_trunc[i, k] = trefa_trunc
        teq_trunc1[i, k] = trefa_trunc1
        teq_trunc2[i, k] = trefa_trunc2

PHI, P = np.meshgrid(clat, pref_mid_norm)
P = P*psurf_ref/100
PHI = PHI*180/np.pi

fig = plt.figure(figsize=(12,5))
ax1 = fig.add_subplot(131)
ax2 = fig.add_subplot(132)
ax3 = fig.add_subplot(133)

CS = ax1.contour(PHI, P, teq.T, colors='k', levels=np.arange(150, 320, 10))
ax1.set_ylim([min(pmid)/100, max(pmid)/100])
ax1.invert_yaxis()
ax1.set_yscale('log')
ax1.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '{:g}'.format(y)))
ax1.set_ylabel('p  [hPa]', fontsize=11)
ax1.set_xlabel('φ  [deg]', fontsize=11)
ax1.clabel(CS, fmt='%.0f', manual=False, levels = [160, 170, 180, 200, 220, 240, 260, 280, 300])
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
ax2.clabel(CS, fmt='%.0f', manual=False, levels = [160, 170, 180, 200, 220, 240, 260, 280, 300])
ax2.plot([-90, 90], [1, 1], '--b')
ax2.yaxis.tick_right()

CS = ax3.contour(PHI, P, teq_trunc.T, colors='k', levels=np.arange(150, 320, 10))
ax3.set_ylim([min(pmid)/100, max(pmid)/100])
ax3.invert_yaxis()
ax3.set_yscale('log')
ax3.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '{:g}'.format(y)))
ax3.set_xlabel('φ  [deg]', fontsize=11)
ax3.clabel(CS, fmt='%.0f', manual=False, levels = [160, 170, 180, 200, 220, 240, 260, 280, 300])
ax3.plot([-90, 90], [1, 1], '--b')
ax3.yaxis.tick_right()


from metpy.units import units as u
cutoff = 100/psurf_ref # 100 Pa in sigma
rev = 1/((3 * u.day).to(u.s)) # 1/3day in s
pih = np.pi/2

num = pih*np.log(cutoff/pref_mid_norm)
den = np.log(cutoff/pref_mid_norm[0])
kr = rev * (np.sin(num/den))**2
for i in range(len(kr)):
    if(pref_mid_norm[i] > cutoff): kr[i] = 0 * 1/u.s

ax3 = ax2.twiny()
ax3.plot(kr.to(1/u.day), pmid/100, '-r', lw=2)
ax3.set_xlabel('RF coefficient [1/day]', fontsize=10, color='r')

# --------------- better figure -----------------------
levels2 = np.hstack([np.arange(100,200,10), np.arange(200, 320, 10)])
levels1 = np.arange(140, 320, 20)
fig = plt.figure(figsize=(7,5))
ax = fig.add_subplot(111)

plotRF = False
if plotRF:
    axx = ax.twiny()
    axx.plot(kr.to(1/u.day), pmid/100, '-r', lw=2)
    axx.set_xlabel('RF damping timescale [1/day]', fontsize=12, color='r')
    axx.tick_params(axis='x', colors='red')
    axx.set_xlim([-0.003, 0.35])

cArgs_c = {'fmt':'%d', 'manual':((-43.2, 228), (-26.57, 305), (0, 399), (26.87, 533), 
                                 (40.55, 814), (0, 870))}
pltargs_c = {'levels':levels1, 'colors':'k', 'linewidths':1.5, 'zorder':1}
var_dict_c = [{'var':teq_trunc1.T, 'plotType':'contour', 'plotArgs':pltargs_c, 'colorArgs':cArgs_c}]
cc = pltvert(clat * 180/np.pi, pref_mid_norm*psurf_ref/100, var_dict_c, ax=ax, 
        plot_zscale=True, annotation='', xlabel='latitude', ylabel='p  [hPa]')
for label in cc[0]:
    label.set_rotation(0)

cArgs_c = {'fmt':'%d', 'manual':((0, 105), (0, 29), (0, 16), (0, 9), (0, 6),(0, 2), 
                                 (-77.467,44), (-77.467, 17), (-80, 6), (-81, 0.3),
                                 (77.467,44), (77.467, 17), (80, 6), (81, 0.3))}
pltargs_c = {'levels':levels2, 'colors':'k', 'linewidths':1.5, 'zorder':1}
var_dict_c = [{'var':teq_trunc2.T, 'plotType':'contour', 'plotArgs':pltargs_c, 'colorArgs':cArgs_c}]
cc = pltvert(clat * 180/np.pi, pref_mid_norm*psurf_ref/100, var_dict_c, ax=ax, 
        plot_zscale=False, inverty=False, annotation='', xlabel='latitude', ylabel='p  [hPa]')


plt.tight_layout()
fig.savefig('WHS_newTeq_noRF.png', dpi=300)
plt.show()

