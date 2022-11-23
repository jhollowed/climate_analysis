import matplotlib.pyplot as plt
import numpy as np
from metpy.units import units as u
import metpy.constants as const
import pdb

# ------- local heating

dTstrat = 0.3 * u.K/u.day
qstar = 9.5e-7
q0 = 1e-12
qplotmax = -4
#pp = np.abs(-12 - np.log10(qstar))
q = np.logspace(-12, qplotmax, 100)
gamma = [8, 10, 12, 14, 16]

# ---- linear form
s = (q/qstar * const.Cp_d * dTstrat).to(u.J/u.kg/u.s)
sd = (s / const.Cp_d).to(u.K/u.day)
plt.plot(q, sd, 'r', label='original heating form')

# ---- log form
#s2 = (np.log10(q)/np.log10(qstar) * const.Cp_d * dT).to(u.J/u.kg/u.s)
#s2 = ((np.log10(q) - np.log10(qstar) + pp) * (1/pp) * const.Cp_d * dT).to(u.J/u.kg/u.s)
#s2 = ((np.log10(q/qstar) + pp) * (1/pp) * const.Cp_d * dT)
for i in range(len(gamma)):
    s = ((1 - np.log10(q/qstar)/np.log10(q0/qstar))**gamma[i] * const.Cp_d * dTstrat).to(u.J/u.kg/u.s)
    sdd = (s / const.Cp_d).to(u.K/u.day)
    plt.plot(q, sdd, label='gamma={}'.format(gamma[i]))

plt.xscale('log')
plt.plot([qstar, qstar], plt.ylim(), '--k', lw=0.75)
plt.xlabel('mixing ratio')
plt.ylabel('local heating rate [K/day]')
plt.grid()
plt.legend()


# ---------- aod cooling

plt.figure()
dTsurf = -0.012 * u.K/u.day
taustar = 1.3e7
tau0 = 1e3
#pp = np.abs(-12 - np.log10(qstar))
tau = np.logspace(2, 8, 100)

# ---- linear form
s = (tau/taustar * const.Cp_d * dTsurf).to(u.J/u.kg/u.s)
sd = (s / const.Cp_d).to(u.K/u.day)
plt.plot(tau, sd, 'r', label='original cooling form')

# ---- log form
for i in range(len(gamma)):
    s = ((1 - np.log10(tau/taustar)/np.log10(tau0/taustar))**gamma[i] * 
              const.Cp_d * dTsurf).to(u.J/u.kg/u.s)
    sdd = (s / const.Cp_d).to(u.K/u.day)
    plt.plot(tau, sdd, label='gamma={}'.format(gamma[i]))

plt.xscale('log')
plt.plot([taustar, taustar], plt.ylim(), '--k', lw=0.75)
plt.xlabel('AOD (column mass burden value) (dimensionless)')
plt.ylabel('local cooling rate [K/day]')
plt.grid()
plt.legend()

plt.show()
