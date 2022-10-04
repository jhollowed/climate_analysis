import matplotlib.pyplot as plt
import numpy as np
from metpy.units import units as u
import metpy.constants as const
import pdb

# ------- local heating

dTstrat = 0.3 * u.K/u.day
qstar = 9.5e-7
q0 = 1e-12
#pp = np.abs(-12 - np.log10(qstar))
q = np.logspace(-12, -5, 100)
gamma = [1, 2, 3, 4]

# ---- linear form
s = (q/qstar * const.Cp_d * dTstrat).to(u.J/u.kg/u.s)
# ---- log form
#s2 = (np.log10(q)/np.log10(qstar) * const.Cp_d * dT).to(u.J/u.kg/u.s)
#s2 = ((np.log10(q) - np.log10(qstar) + pp) * (1/pp) * const.Cp_d * dT).to(u.J/u.kg/u.s)
#s2 = ((np.log10(q/qstar) + pp) * (1/pp) * const.Cp_d * dT)
s2 = ((1 - np.log10(q/qstar)/np.log10(q0/qstar))**gamma[0] * const.Cp_d * dTstrat).to(u.J/u.kg/u.s)
# ---- log^2 form
s3 = ((1 - np.log10(q/qstar)/np.log10(q0/qstar))**gamma[1] * const.Cp_d * dTstrat).to(u.J/u.kg/u.s)
# ---- log^3 form
s4 = ((1 - np.log10(q/qstar)/np.log10(q0/qstar))**gamma[2] * const.Cp_d * dTstrat).to(u.J/u.kg/u.s)
# ---- log^4 form
s5 = ((1 - np.log10(q/qstar)/np.log10(q0/qstar))**gamma[3] * const.Cp_d * dTstrat).to(u.J/u.kg/u.s)

sd = (s / const.Cp_d).to(u.K/u.day)
s2d = (s2 / const.Cp_d).to(u.K/u.day)
s3d = (s3 / const.Cp_d).to(u.K/u.day)
s4d = (s4 / const.Cp_d).to(u.K/u.day)
s5d = (s5 / const.Cp_d).to(u.K/u.day)


plt.plot(q, sd, 'r', label='original heating form')
plt.plot(q, s2d, 'b', label='log heating form')
plt.plot(q, s3d, 'g', label='log^2 heating form')
plt.plot(q, s4d, 'm', label='log^3 heating form')
plt.plot(q, s5d, 'k', label='log^4 heating form')
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
gamma = [1, 2, 3, 4]

# ---- linear form
s = (tau/taustar * const.Cp_d * dTsurf).to(u.J/u.kg/u.s)
# ---- log form
s2 = ((1 - np.log10(tau/taustar)/np.log10(tau0/taustar))**gamma[0] * const.Cp_d * dTsurf).to(u.J/u.kg/u.s)
# ---- log^2 form
s3 = ((1 - np.log10(tau/taustar)/np.log10(tau0/taustar))**gamma[1] * const.Cp_d * dTsurf).to(u.J/u.kg/u.s)
# ---- log^3 form
s4 = ((1 - np.log10(tau/taustar)/np.log10(tau0/taustar))**gamma[2] * const.Cp_d * dTsurf).to(u.J/u.kg/u.s)
# ---- log^4 form
s5 = ((1 - np.log10(tau/taustar)/np.log10(tau0/taustar))**gamma[3] * const.Cp_d * dTsurf).to(u.J/u.kg/u.s)

sd = (s / const.Cp_d).to(u.K/u.day)
s2d = (s2 / const.Cp_d).to(u.K/u.day)
s3d = (s3 / const.Cp_d).to(u.K/u.day)
s4d = (s4 / const.Cp_d).to(u.K/u.day)
s5d = (s5 / const.Cp_d).to(u.K/u.day)


plt.plot(tau, sd, 'r', label='original cooling form')
plt.plot(tau, s2d, 'b', label='log cooling form')
plt.plot(tau, s3d, 'g', label='log^2 cooling form')
plt.plot(tau, s4d, 'm', label='log^3 cooling form')
plt.plot(tau, s5d, 'k', label='log^4 cooling form')
plt.xscale('log')
plt.plot([taustar, taustar], plt.ylim(), '--k', lw=0.75)
plt.xlabel('AOD (column mass burden value) (dimensionless)')
plt.ylabel('local cooling rate [K/day]')
plt.grid()
plt.legend()

plt.show()
