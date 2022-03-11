import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pdb
import matplotlib.ticker as ticker

psurf_ref = 100000 # Pa
pref_mid_norm = (np.logspace(-2, 3, 1000)*100)/psurf_ref
cutoff = 100/psurf_ref # 100 Pa in sigma
rev = 1/(86400*3) # 1/3day
pih = np.pi/2

num = pih*np.log(cutoff/pref_mid_norm)
den = np.log(cutoff/pref_mid_norm[0])
kr = rev * (np.sin(num/den))**2
for i in range(len(kr)):
    if(pref_mid_norm[i] > cutoff): kr[i] = 0

plt.plot(kr, pref_mid_norm, 'k', lw=2)
plt.ylim([min(pref_mid_norm), max(pref_mid_norm)])
plt.gca().invert_yaxis()
plt.yscale('log')
plt.gca().yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '{:g}'.format(y)))
plt.ylabel('p/p0', fontsize=14)
plt.xlabel('RF coefficient k_r  [1/s]', fontsize=14)
plt.savefig('ttt.png', dpi=300)
