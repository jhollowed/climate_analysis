import matplotlib.pyplot as plt
import numpy as np
from metpy.units import units as u
import metpy.constants as const
import pdb

dT = 0.3 * u.K/u.day
qstar = 9.5e-7
q0 = 1e-12
#pp = np.abs(-12 - np.log10(qstar))
q = np.logspace(-12, -5, 100)

s = (q/qstar * const.Cp_d * dT).to(u.J/u.kg/u.s)
#s2 = (np.log10(q)/np.log10(qstar) * const.Cp_d * dT).to(u.J/u.kg/u.s)
#s2 = ((np.log10(q) - np.log10(qstar) + pp) * (1/pp) * const.Cp_d * dT).to(u.J/u.kg/u.s)
#s2 = ((np.log10(q/qstar) + pp) * (1/pp) * const.Cp_d * dT)
s2 = ((1 - np.log10(q/qstar)/np.log10(q0/qstar)) * const.Cp_d * dT).to(u.J/u.kg/u.s)

sd = (s / const.Cp_d).to(u.K/u.day)
s2d = (s2 / const.Cp_d).to(u.K/u.day)


plt.plot(q, sd, 'r', label='original heating form')
plt.plot(q, s2d, 'b', label='log heating form')
plt.xscale('log')
plt.plot([qstar, qstar], plt.ylim(), '--k', lw=0.75)
plt.xlabel('mixing ratio')
plt.ylabel('heating rate for 1 kg air mass [K/day]')
plt.grid()
plt.legend()

plt.show()
