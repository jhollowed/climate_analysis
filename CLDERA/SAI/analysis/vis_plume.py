
import numpy as np
import matplotlib.pyplot as plt
import metpy.constants as const
from metpy.units import units as u
import pdb
import math
import cartopy.crs as ccrs
import matplotlib.ticker as mticker
import metpy.calc as mc
from matplotlib.ticker import ScalarFormatter
import scipy
import matplotlib as mpl
from climate_toolbox import climate_artist as cla

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica"]})
#mpl.rcParams['text.latex.preamble'] = [
#       r'\usepackage{siunitx}',   # i need upright \micro symbols, but you need...
#       r'\sisetup{detect-all}',   # ...this to force siunitx to actually use your fonts
#       r'\usepackage{helvet}',    # set the normal font here
#       r'\usepackage{sansmath}',  # load up the sansmath so that math -> helvet
#       r'\sansmath'               # <- tricky! -- gotta actually tell tex to use!
#]
params = {'text.usetex': False, 'mathtext.fontset': 'stixsans'}
plt.rcParams.update(params)

# ========== define parameters ==========
d2r = np.pi/180                
r2d = 180/np.pi

lat0 = 15.15 * d2r
lon0 = 120.35 * d2r
z0 = 25000 * u.m

dr = 100000 * u.m
#dz = 7500 * u.m
dz = 5000 * u.m
rs = dr*3
zs = dz*3

tf = 172800 * u.s
tau = -np.log(0.05)/tf
a = const.Re.to(u.m)

M_SO2 = 2e10 * u.kg
M_ash = 2e10 * u.kg
k_SO2 = 1/2592000 * 1/u.s
k_ash = 1/86400 * 1/u.s



# ========== define variables ==========
lat = np.linspace(10, 20, 300) * d2r
lon = np.linspace(115, 125, 302) * d2r
LAT,LON = np.meshgrid(lat, lon)


z = np.linspace(0, 50000, 303) * u.m
t = np.linspace(0, 2*tf.m, 100) * u.s

LAT, LON, Z = np.meshgrid(lat, lon, z.m)
Z = Z * u.m

R = (a+z0) * np.arccos( np.sin(LAT)*np.sin(lat0) + np.cos(LAT)*np.cos(lat0) * np.cos(np.abs(LON-lon0)) ) #great circle

# ========== construct source term ==========
H = np.exp(-(1/2)*((R)/dr)**2)      # horizontal
V = np.exp(-(1/2)*((Z-z0)/dz)**2)   # vertical

Te = np.exp(-tau*t)                 # temporal (exponential)
Tc = t <= tf                        # temporal (constant)

# -----  truncate
H  *= (R<=rs).astype(int)
Z  *= (abs(Z-z0)<=zs).astype(int)
Te *= (t <= tf).astype(int)

# ----- normalization

# time constants
alphae = ((1-np.exp(-tau*tf))/tau)
alphac = tf

# mass scaling factor
dne = alphae * np.sqrt(2*np.pi**3) * dr**2 * dz
dnc = alphac * np.sqrt(2*np.pi**3) * dr**2 * dz

# H, V integrations
Hint = (1 - np.exp(-rs**2/(2*dr**2)))
Vint = (math.erf(z0/(np.sqrt(2)*dz)) - math.erf((z0-zs)/(np.sqrt(2)*dz)))
HVint = Hint * Vint

# norm
Ae_SO2 = (M_SO2/dne) * HVint 
Ac_SO2 = (M_SO2/dnc) * HVint 
Ae_ash = (M_ash/dne) * HVint 
Ac_ash = (M_ash/dnc) * HVint 

# ----- exp source at t=0
f = H * V
fe_SO2 = Ae_SO2 * f

# ========== peak concentration over time==========
fpe_SO2 = Ae_SO2 * Te         # fpe = forcing peak, exponential in t
fpc_SO2 = Ac_SO2 * Tc
fpe_ash = Ae_ash * Te        # fpe = forcing peak, exponential in t
fpc_ash = Ac_ash * Tc

# solve ode at peak 
# rho(r=0, z=z0, t) = -k * rho(0, 0, t) + f(0, 0, t)
# with init value rho(0, 0, 0) = 0

# old numerical solution
print('scipy')

ud = u.kg/u.m**3
tb = lambda tt: int(tt <= tf.m)

drdt_eSO2 = lambda r, tt: -k_SO2 * r*ud + Ae_SO2 * np.exp(-tau*tt*u.s) * tb(tt)
rho_peak_eSO2_num = scipy.integrate.odeint(drdt_eSO2, 0, t)
drdt_cSO2 = lambda r, tt: -k_SO2 * r*ud + Ac_SO2 * int(tt*u.s <= tf)
rho_peak_cSO2_num = scipy.integrate.odeint(drdt_cSO2, 0, t)
drdt_eash = lambda r, tt: -k_ash * r*ud + Ae_ash * np.exp(-tau*tt*u.s) * tb(tt)
rho_peak_eash_num = scipy.integrate.odeint(drdt_eash, 0, t)
drdt_cash = lambda r, tt: -k_ash * r*ud + Ac_ash * int(tt*u.s <= tf)
rho_peak_cash_num = scipy.integrate.odeint(drdt_cash, 0, t)
print('done')

# new analytic solution
tbound = np.minimum(t, tf)
rho_peak_eash = (Ae_ash * np.exp(-k_ash*t) * ( -1+np.exp(tbound*(k_ash-tau)))) / (k_ash-tau)
rho_peak_cash = (Ac_ash * np.exp(-k_ash*t) * ( -1+np.exp(tbound*k_ash))) / k_ash

rho_peak_eSO2 = (Ae_SO2 * np.exp(-k_SO2*t) * ( -1+np.exp(tbound*(k_SO2-tau)))) / (k_SO2-tau)
rho_peak_cSO2 = (Ac_SO2 * np.exp(-k_SO2*t) * ( -1+np.exp(tbound*k_SO2))) / k_SO2




# ========== plot ==========

data_crs = ccrs.PlateCarree()
fig = plt.figure(figsize=(8.2,6))
spec = fig.add_gridspec(2, 5)
#ax1 = fig.add_subplot(spec[0,0], projection=ccrs.AlbersEqualArea(lon0, lat0))
ax1 = fig.add_subplot(spec[1,0:2], projection=ccrs.PlateCarree(lon0))
ax2 = fig.add_subplot(spec[0,0:2])
ax3 = fig.add_subplot(spec[0,2:])
ax4 = fig.add_subplot(spec[1,2:])

zmid = np.searchsorted(z, z0)
latmid = np.searchsorted(lat, lat0)
rr = ((LON - lon0) * (a+z0)).to(u.km)

P0 = 1000*u.hPa
T0 = 250*u.K
H = const.Rd*T0/const.g
p = P0 * np.exp(-Z/H)
pmid = P0 * np.exp(-z[zmid]/H)

cmap = plt.cm.OrRd
#levels = [4, 3, 2, 2.5, 1, 0.5]
levels = [0, 0.5, 1, 1.5, 2, 3, 4]
poww = math.ceil(-np.log10(np.max(fe_SO2[:,latmid,:]).m))
vmin=-10**-poww

ax1.contour(LON[:,:,zmid]*r2d, LAT[:,:,zmid]*r2d, fe_SO2[:,:,zmid]*10**poww, levels=levels, transform=data_crs, cmap=cmap, vmin=vmin)
cs = ax2.contour(rr[:,latmid,:], p[:,latmid,:], fe_SO2[:,latmid,:]*10**poww, levels=levels, cmap=cmap, vmin=vmin)

degd = 4
extent = [lon0*r2d-degd, lon0*r2d+degd, lat0*r2d-degd, lat0*r2d+degd]
lonlim = np.arange(int(extent[0]), int(extent[1]), 2)
latlim = np.arange(int(extent[2]), int(extent[3]), 2)
ax1.set_extent(extent)
gl = ax1.gridlines(draw_labels=True, zorder=0)
gl.xlocator = mticker.FixedLocator(lonlim)
gl.ylocator = mticker.FixedLocator(latlim)
gl.right_labels = []
gl.top_labels = []
ax1.coastlines(resolution='50m', color='k', linestyle='-', alpha=1, zorder=9)
#ax1.set_xlabel(r'lat [deg]')
#ax1.set_ylabel(r'lon [deg]')
ax1.set_title('{:.0f} km, $\sim$40 hPa'.format(z[zmid].to(u.km).m))

ax2.set_xlabel(r'r [km]')
ax2.set_ylabel(r'p [hPa]')
ax2.set_xlim([-300, 300])
ax2.set_yscale('log')
ax2.invert_yaxis()
ax2.yaxis.set_major_formatter(ScalarFormatter())
ax2.clabel(cs, inline=1, fontsize=10)#, fmt='%1.0f')
ax2.set_title(r'$d\rho/dt \times 10^{{{pp}}}$ ($t$=0)'.format(pp=poww))

ax22 = ax2.twinx()
ax22.set_ylabel(r'Z [km]')
ax22.contour(rr[:,latmid,:], Z[:,latmid,:].to(u.km), fe_SO2[:,latmid,:], levels=15, colors='k', alpha=0)
#ax22.set_yscale('log')
#ax22.invert_yaxis()
#ax2.yaxis.set_major_formatter(ScalarFormatter())
#ptick = ax2.get_yticks()
#ax22.set_yticks(ptick)
ax22.set_ylim( (H * np.log(P0/(ax2.get_ylim()*u.hPa))).to(u.km).m )
#ztick = (H * np.log(P0/(ptick*u.hPa))).to(u.km).m
#ax22.set_yticklabels(['{:.0f}'.format(ztick[i]) for i in range(len(ztick))])

ax3.plot(t.m/3600, fpe_SO2.m, '-k', label=r'exponential $T(t)$')
ax3.plot(t.m/3600, fpc_SO2.m, '--k', label=r'constant $T(t)$')
ax3.set_ylabel(r'peak $f(t)$ [kg/m$^3$/s]', fontsize=11)
ax3.set_xticklabels([])
ax3.yaxis.tick_right()
ax3.yaxis.set_label_position("right")
ax3.legend(frameon=False)
yy = ax3.get_ylim()
yy = [yy[0], yy[1]]
tftf = [tf.to(u.hr).m, tf.to(u.hr).m]
ax3.plot(tftf, yy, ':k', lw=0.8)
ax3.set_ylim(yy)
#ax3.text(tf.to(u.hr).m-0.5, -0.7e-10, r'$t_f$')

ax4.plot(t.m/3600, rho_peak_eSO2, '-r', label='SO2')
ax4.plot(t.m/3600, rho_peak_cSO2, '--r')
ax4.plot(t.m/3600, rho_peak_eash, '-c', label='Ash')
ax4.plot(t.m/3600, rho_peak_cash, '--c')

ax4.plot(t.m/3600, rho_peak_eSO2_num, '-k', lw=0.6, label='numerical')
ax4.plot(t.m/3600, rho_peak_cSO2_num, '-k', lw=0.6)
ax4.plot(t.m/3600, rho_peak_eash_num, '-k', lw=0.6)
ax4.plot(t.m/3600, rho_peak_cash_num, '-k', lw=0.6)

ax4.set_xlabel(r'time [hr]')
ax4.set_ylabel(r'peak $\rho(t)$ [kg/m$^3$]', fontsize=11)
ax4.yaxis.tick_right()
ax4.yaxis.set_label_position("right")
ax4.legend(frameon=False)
yy = ax4.get_ylim()
yy = [yy[0], yy[1]]
ax4.plot(tftf, yy, ':k', lw=0.8)
ax4.set_ylim([0, yy[1]])

ticks = ax3.get_xticks().tolist()
ticklabs = ax3.get_xticklabels()
ticks.insert(np.searchsorted(ticks, tf.m/3600), tf.m/3600)
ticklabs.insert(np.searchsorted(ticks, tf.m/3600), '$t_f$')
ax4.set_xticks(ticks)
ax4.set_xticklabels(ticklabs)
ax3.set_xticks(ticks)

for ax in [ax2, ax3, ax4]:
    ax.xaxis.set_ticks_position('both')
    ax.xaxis.set_tick_params(direction='in')
    if(ax != ax2):
        ax.yaxis.set_ticks_position('both')
    ax.yaxis.set_tick_params(direction='in')
ax3.set_xlim(-0.5, np.max(t.m/3600))
ax4.set_xlim(-0.5, np.max(t.m/3600))

plt.tight_layout()
#plt.savefig('plume.png', dpi=300)
plt.show()



