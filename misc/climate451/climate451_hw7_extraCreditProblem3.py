import pdb
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from metpy.units import units
from metpy import constants as const

# --- constants
cp = const.Cp_d
g  = const.g
Rd = const.Rd
k  = Rd/cp

# --- read data
data  = xr.open_dataset('CAM_FV3L30_HS_3D.nc')
data  = data.assign_coords(lev=data.lev*100).isel(time=0)
T     = data.T * units.K
u     = data.U * units.m/units.s
v     = data.V * units.m/units.s
ps    = data.PS * units.Pa
p0    = data.P0 * units.Pa
lat   = data.lat * units.deg
# get pressure as 3d field
dims_2d = {'lat':len(T.lat), 'lon':len(T.lon)}
p       = data.lev.expand_dims(dims_2d).broadcast_like(data.T) * units.Pa

# --- compute potential temperature
theta = T * (p0/p)**(Rd/cp)

# --- compute TPE
Tdp = T.integrate('lev') * units.Pa
TPE = cp/g * Tdp

# --- compute KE
KE = 1/(2*g) * (u**2 + v**2)
KE = KE.integrate('lev') * units.Pa

# --- compute APE with zonal averaging
theta_zm       = theta.mean('lon')
theta_eddy     = theta - theta_zm
dthetazm_dp    = theta_zm.differentiate('lev') / units.Pa
theta_eddy2_zm = (theta_eddy**2).mean('lon')
integrand      = p**(k-1) * (-g * dthetazm_dp)**(-1) * theta_eddy2_zm
APE            = Rd*ps**(-k)/2 * integrand.integrate('lev') * units.Pa

# --- compute APE with global horizontal averaging
horz_mean = lambda x: x.weighted(np.cos(np.deg2rad(lat))).mean('lat', 'lon')
theta_avg       = horz_mean(theta)
theta_eddy_avg  = theta - theta_avg
dthetaavg_dp    = theta_avg.differentiate('lev') / units.Pa
theta_eddy2_avg = horz_mean(theta_eddy_avg**2)
p_avg           = horz_mean(p)
integrand       = p_avg**(k-1) * (-g * dthetaavg_dp)**(-1) * theta_eddy2_avg
APE2d           = Rd*ps**(-k)/2 * integrand.integrate('lev') * units.Pa

# --- compute APE with zonal averaging, using p0 rather than ps
APE2d0          = Rd*p0**(-k)/2 * integrand.integrate('lev') * units.Pa

# --- compute zonal means
TPE_zm    = TPE.mean('lon')
KE_zm     = KE.mean('lon')
APE_zm    = APE.mean('lon')
APE2d_zm  = APE2d.mean('lon')
APE2d0_zm = APE2d0.mean('lon')
print('APE 2D w/ P0 = {:.3e}'.format(APE2d0_zm))

# --- plot
for var,dat in {'TPE':TPE_zm, 'KE':KE_zm, 'APE':APE_zm, 'APE 2D':APE2d_zm}.items():
    fig = plt.figure(figsize=(6,3))
    ax = fig.add_subplot()
    ax.plot(data.lat, dat, {'TPE':'b','KE':'r','APE':'g', 'APE 2D':'g'}[var], lw=2)
    ax.set_xlim(-90, 90)
    ax.set_xlabel('lat', fontsize=14)
    ax.set_ylabel('zonal-mean {} [J/m\u00B2]'.format(var), fontsize=14)
    ax.grid()
    plt.tight_layout()
    plt.savefig('figs/{}.pdf'.format(var), dpi=100)

# --- compare TPE, KE, APE at 50N
TPE_50N = TPE_zm.sel(lat=50, method='nearest')
KE_50N  = KE_zm.sel(lat=50, method='nearest')
APE_50N = APE_zm.sel(lat=50, method='nearest')
APE2d_50N = APE2d_zm.sel(lat=50, method='nearest')
print('TPE/APE = {:.2e}/{:.2e} = {:.2f}'.format(TPE_50N, APE_50N, TPE_50N/APE_50N))
print('TPE/APE2d = {:.2e}/{:.2e} = {:.2f}'.format(TPE_50N, APE2d_50N, TPE_50N/APE2d_50N))
print('TPE/KE = {:.2e}/{:.2e} = {:.2f}'.format(TPE_50N, KE_50N, TPE_50N/KE_50N))

# ---- compare the zonal-mean and 2d-averaged evaluations of the APE
fig = plt.figure(figsize=(6,3))
ax = fig.add_subplot()
ax.plot(data.lat, APE2d_zm/APE_zm, 'k', lw=2)
ax.set_xlim(-90, 90)
ax.set_ylabel(r'$\frac{\text{APE w/ 2D averaging}}{\text{APE w/ zonal averaging}}$', fontsize=16)
ax.set_xlabel('lat', fontsize=14)
ax.grid()
plt.tight_layout()
ax.set_yscale('log')
plt.savefig('figs/APE_comparison.pdf', dpi=100)

plt.show()
