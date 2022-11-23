import matplotlib.pyplot as plt
import numpy as np
from metpy.units import units as u
import metpy.constants as const
import pdb

# ------- local heating

dTstrat = 0.3 * u.K/u.day
qstar = 9.5e-7
q0 = 1e-12
qplotmax = 0
dTplotmax = 1
#pp = np.abs(-12 - np.log10(qstar))
q = np.logspace(-12, qplotmax, 100)
gamma = [8, 10, 12]

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

#--- sigmoid form
k = [1, 2, 3, 4, 5]
for i in range(len(k)):
    sig = lambda qq: 2 / (1 + np.exp(-k[i]*np.log10(qq/qstar)))
    dsig = lambda qq: \
           k[i]*2*np.exp(-k[i]*np.log10(qq/qstar)) / (1 + np.exp(-k[i]*np.log10(qq/qstar)))**2
    s = sig(q)
    #s[q > qstar] = ((dsig(qstar)) * np.log10(q/qstar) + 1)[q>qstar]
    #s = s * np.abs((np.log10(qstar) - 7)/(np.log10(q) - 2))
    #s = s + 0.5 * np.log10(q)-np.log10(qstar)
    s = (s * dTstrat * const.Cp_d).to(u.J/u.kg/u.s)
    sddd = (s / const.Cp_d).to(u.K/u.day)
    plt.plot(q, sddd, '--', label='k={}'.format(k[i]))

#--- asinh form
#b = [6, 8, 100]
#q0 = 1e-10
#for i in range(len(b)):
#    s = (dTstrat * (1 - (np.arcsinh(b[i]*np.log10(qstar/q)) / np.arcsinh(b[i]*np.log10(qstar/q0)))) \
#         * np.abs(np.log10(qstar)/np.log10(q)) * const.Cp_d).to(u.J/u.kg/u.s)
#    sdddd = (s / const.Cp_d).to(u.K/u.day)
#    plt.plot(q, sdddd, ':', label='b={}'.format(b[i]))

plt.xscale('log')
plt.plot([qstar, qstar], plt.ylim(), '--k', lw=0.75)
plt.xlabel('mixing ratio')
plt.ylabel('local heating rate [K/day]')
plt.ylim([-0.2, dTplotmax])
plt.grid()
plt.legend()

# tmp
plt.show()

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
