import pdb
import numpy as np
import matplotlib.pyplot as plt
import metpy.constants as const
import astropy.constants as aconst
from metpy.units import units as u

molecules_to_moles = 1/(aconst.N_A.value * (1/u.mol))
molar_mass = 28 * (u.g/u.mol)
cflx = 2.7736e11 * molecules_to_moles * (u.cm**-2 * u.s**-1)   # from Abalos
cflx = (cflx * molar_mass).to(u.kg * u.m**-2 * u.s**-1)

dt = 1800 * u.s
dp = 1.5 * u.hPa
tmp = dt * const.g/dp

t = np.linspace(0, 90*86400, 10000) * u.s
q = np.zeros(len(t))

for i in range(len(t) - 1):
    if(i%100 == 0): print(i)
    q[i+1] = q[i] + (q[i] * -1/(90*86400*u.s) * dt)
    q[i+1] += q[i] + (cflx * tmp).to_base_units()

plt.ion()
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(t/3600, q * 1e9, '-r')
ax.set_xlabel('time  [hr]')
ax.set_ylabel('E90 mixing ratio [ppb]')
plt.show()
pdb.set_trace()


