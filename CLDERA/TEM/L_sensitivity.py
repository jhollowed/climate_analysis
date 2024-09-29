import numpy as np
import matplotlib.pyplot as plt

#sph_err = np.load('/global/homes/j/jhollo/repos/PyTEMDiags/maxerr.npy')
sph_err = np.load('./data/maxerr.npy')
LL = LL = [25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 325, 350, 375, 400, 425, 450]

plt.plot(LL, sph_err, '-d')
plt.plot(LL, np.ones(len(LL))*0.8857934713616029, '--r')
plt.plot(LL, np.ones(len(LL))*1.2368602741338992, '--r')

plt.ylim([-2, 5])

plt.show()

