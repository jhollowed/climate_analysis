import pdb
import numpy as np
import metpy.constants as const
from metpy.units import units as u
import matplotlib.pyplot as plt

Tmax = 315 * u.K
dTy = 60 * u.K
sb = 5.670374419e-8 * u.W/(u.m**2 * u.K**4)
Re = const.Re
eps = 1 # emissivity

lat = np.linspace(-np.pi/2, np.pi/2, 1000)
latd = lat * 180/np.pi
Teq = (Tmax - dTy*np.sin(lat)**2)

Ilw = eps * sb * Teq**4

integ_lw = eps * sb * 4*np.pi * Re**2 * (7736342625*u.K**4)

S = integ_lw / (np.pi*Re)**2
Isw = S * np.cos(lat)

integ_sw = (np.pi*Re)**2 * S

# -------- heating balance fig

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(latd, Ilw, '-r', label='outgoing longwave radiation')
ax.plot(latd, Isw, '-b', label='absorbed solar radiation')
ax.plot(latd, Isw - Ilw, '-k', lw=1, label='Net radiation')
ax.plot([-90, 90], [0, 0], '--k', lw=0.75)
ax.set_xlim([-90, 90])
ax.grid(alpha=0.4)
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
ax.legend(loc='upper left', frameon=False)
ax.set_xlabel('Latitude')
ax.set_ylabel('Flux Density (W/m2)')
ax.set_ylim([-250, 725])

plt.savefig('toa.png', dpi=300)


# --------- heating functions with AOD, vmr

g = const.g
cp = const.Cp_d
cell_area= (200*u.km)**2
dp_surf = 20*u.hPa # 23 hPa layer at the surface ~ 200m ~ 3 levels in E3SM
dp_strat = 8*u.hPa
m_strat = (dp_strat * cell_area/g).to(u.kg)
m_surf = (dp_surf * cell_area/g).to(u.kg)
p_strat = (m_strat / (cell_area * 4*u.km)).to(u.kg/u.m**3)
blw = 1e-2 * u.m**2/u.kg
bsw = 100 * u.m**2/u.kg


lat = 0 * np.pi/180
Teq = (Tmax - dTy*np.sin(lat)**2)
Ilw = eps * sb * Teq**4
Isw = S * np.cos(lat)

# strat tuning estimate
dT = 0.3*u.K/u.day
blw = m_strat*dT*cp/(cell_area*Ilw) * (g/(1e-4*dp_strat))
blw = blw.to(u.m**2/u.kg)
pdb.set_trace()

# surf tuning estimate
dT = -0.02*u.K/u.day
bsw = (0.2*cell_area/(2e7*u.kg)).to(u.m**2/u.kg)
#zeta = 0.0005 # zetaiciency of heating rate to surface -> heatinr rate to atm
zeta = 1 / (1/cp * cell_area/m_surf * 1/dT * Isw * (np.exp(-0.2) - 1))
zeta = zeta.to_base_units()


aod = np.linspace(0, 0.2, 100)
#aod = np.linspace(0, 10, 1000)
dI_surf = Isw * (np.exp(-aod) - 1)
s_surf  = zeta * (cell_area * dI_surf / (m_surf)).to(u.J/(u.kg*u.s))
dT_surf = (s_surf/cp).to(u.K/u.day)
dI_surf_approx = Isw * -(aod)
s_surf_approx  = zeta * (cell_area * dI_surf_approx / (m_surf)).to(u.J/(u.kg*u.s))
dT_surf_approx = (s_surf_approx/cp).to(u.K/u.day)

#q = np.logspace(-10, -2, 1000)
q = np.logspace(-10, -2, 100)
dI_strat = Ilw * (1 - np.exp(-blw * q * dp_strat / g))
s_strat = (cell_area * dI_strat / (m_strat)).to(u.J/(u.kg*u.s))
dT_strat = (s_strat/cp).to(u.K/u.day)
dI_strat_approx = Ilw * (blw * q * dp_strat / g)
s_strat_approx = (cell_area * dI_strat_approx / (m_strat)).to(u.J/(u.kg*u.s))
dT_strat_approx = (s_strat_approx/cp).to(u.K/u.day)

plt.ion()
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(aod, dT_surf, '-r')
ax.plot(aod, dT_surf_approx, '--r')
ax.grid(alpha=0.4)
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(aod, dI_surf, '-r')
ax.plot(aod, dI_surf_approx, '--r')
ax.grid(alpha=0.4)
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
fig = plt.figure(figsize=(6, 2))
ax = fig.add_subplot(111)
ax.plot(aod, (np.abs((dI_surf_approx-dI_surf) / dI_surf) * 100), '--k', label='linear approx')
ax.set_xlabel('AOD')
ax.set_ylabel('relative error [%]')
ax.grid(alpha=0.4)
plt.tight_layout()
plt.savefig('aod_approx_err.png', dpi=300)

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
ax2.plot(q, dT_strat, '-b')
ax2.plot(q, dT_strat_approx, '--b')
ax2.set_xscale('log')
ax2.grid(alpha=0.4)
ax2.set_ylim([0, 6])
ax2.yaxis.set_ticks_position('both')
ax2.xaxis.set_ticks_position('both')
fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
ax2.plot(q, dI_strat, '-b')
ax2.plot(q, dI_strat_approx, '--b')
ax2.set_xscale('log')
ax2.grid(alpha=0.4)
ax2.set_ylim([0, 6])
ax2.yaxis.set_ticks_position('both')
ax2.xaxis.set_ticks_position('both')
fig = plt.figure(figsize=(6, 2))
ax = fig.add_subplot(111)
ax.plot(q, (np.abs((dI_strat_approx-dI_strat) / dI_strat) * 100), '--k', label='linear approx')
ax.set_xlabel('stratospheric sulfate mixing ratio')
ax.set_ylabel('relative error [%]')
ax.legend(fancybox=False)
ax.grid(alpha=0.4)
ax.set_xscale('log')
plt.tight_layout()
plt.savefig('q_approx_err.png', dpi=300)

fig = plt.figure(figsize=(6, 3))
ax = fig.add_subplot(111)
aod = np.linspace(0, 5, 100)
dI_surf = Isw * (np.exp(-aod) - 1)
ax.plot(aod, dI_surf, '-k')
ax.set_xlabel('AOD')
ax.set_ylabel('Attenuated SW flux\ndensity  [W/m2]')
ax.grid(alpha=0.4)
ax.set_xlim([0, 4.5])
plt.tight_layout()
plt.savefig('aod_sat.png', dpi=300)




pdb.set_trace()
