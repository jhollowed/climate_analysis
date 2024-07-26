import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import glob
import pdb
import climate_toolbox as ctb
import matplotlib as mpl
from matplotlib.colors import ListedColormap
from matplotlib.colors import BoundaryNorm

# the lat-lon input data that was ran through PyTEMDiags
dat = '/pscratch/sd/j/jhollo/E3SM/historical_data/TEM_test_data/v2.LR.WCYCL20TR.pmcpu.limvar.ens1.eam.h1.1998-11-23_TEM_VARIABLES_remap_180x360_aave.nc'
dat = xr.open_dataset(dat)

# the TEM data that resulted from running PyTEMDiags on this lat-lon data
jhtem = '/pscratch/sd/j/jhollo/PyTEMDiags/output/TEM_180x360_1.0deg_pcoord_L45_polesFalse_attrsTrue.nc'
jhtem = xr.open_dataset(jhtem)

# the TEM data that resulted from running CJs Fortran code on this lat-lon data
# original version of data supplied by CJ:
#cjtem = '/pscratch/sd/j/jhollo/PyTEMDiags/output/E3SMv2.LR.WCYCL20TR.pmcpu.limvar.ens1.eam.h1.1998-11-23_TEM_VARIABLES_remap_180x360_aave.tem_p_converted_z.nc'
# re-running with my modifications to output more variables:
cjtem='/pscratch/sd/j/jhollo/E3SM/historical_data/TEM_test_data/v2.LR.WCYCL20TR.pmcpu.limvar.ens1.eam.h1.1998-11-23_TEM_VARIABLES_remap_180x360_aave.tem_p.nc'
cjtem = xr.open_dataset(cjtem)

# reduce data and results...
jhtem = jhtem.mean('time')
cjtem  = cjtem.mean(['time', 'lon'])
lev    = np.intersect1d(jhtem.lev, cjtem.lev)
lat    = np.intersect1d(jhtem.lat, cjtem.lat)
jhtem = jhtem.sel(lev=lev)
jhtem = jhtem.sel(lat=lat)
cjtem  = cjtem.sel(lev=lev)
cjtem  = cjtem.sel(lat=lat)

# TMP: ------ truncate to NH
jhtem = jhtem.sel(lat=slice(0, 90))
cjtem = cjtem.sel(lat=slice(0, 90))
lat = jhtem.lat

dat    = dat.mean(['time', 'lon'])
dat    = dat.sel(lev=lev)
dat    = dat.sel(lat=lat)

# -----------------------------

# because the input was STRUCTURED, we can do verification. The arithmetic zonal mean on the
# input latlon data should be the same as the zonal means computed in PyTEMDiags for e.g. OMEGA

# top row of plots is jhtem, bottom row is cj_tem
# figure will include:
# - wtem
# - zonal mean omega
# - zonal mean v'theta'
# - psi
# - dpsicoslat_dlat

print('creating figure...')
fig = plt.figure(figsize=(22, 8))
ax0_jh = fig.add_subplot(3, 7, 1)
ax1_jh = fig.add_subplot(3, 7, 2)
ax2_jh = fig.add_subplot(3, 7, 3)
ax3_jh = fig.add_subplot(3, 7, 4)
ax4_jh = fig.add_subplot(3, 7, 5)
ax5_jh = fig.add_subplot(3, 7, 6)
ax6_jh = fig.add_subplot(3, 7, 7)
ax0_cj = fig.add_subplot(3, 7, 8)
ax1_cj = fig.add_subplot(3, 7, 9)
ax2_cj = fig.add_subplot(3, 7, 10)
ax3_cj = fig.add_subplot(3, 7, 11)
ax4_cj = fig.add_subplot(3, 7, 12)
ax5_cj = fig.add_subplot(3, 7, 13)
ax6_cj = fig.add_subplot(3, 7, 14)
ax0_diff = fig.add_subplot(3, 7, 15)
ax1_diff = fig.add_subplot(3, 7, 16)
ax2_diff = fig.add_subplot(3, 7, 17)
ax3_diff = fig.add_subplot(3, 7, 18)
ax4_diff = fig.add_subplot(3, 7, 19)
ax5_diff = fig.add_subplot(3, 7, 20)
ax6_diff = fig.add_subplot(3, 7, 21)


# ----------------- vertical residual velocity -------------------
print('working on vertical residual velocity...')

jh_wtem = jhtem['wtem'] * 1e3
cj_wtem = cjtem['wtem'] * 1e3
diff = (jh_wtem.T - cj_wtem) / np.abs(np.maximum(jh_wtem.T, cj_wtem))

levels = np.arange(-0.005, 0.0051, 0.001) * 1e3
cmap='rainbow'

cf = ax0_cj.contourf(lat, lev, cj_wtem, levels=levels, cmap=cmap, extend='both')
plt.colorbar(cf, ax=ax0_cj, location='right', label='wtem [mm/s]')
ax0_cj.set_title('CJ')
ax0_cj.set_xlabel('lat')
ax0_cj.set_ylabel('lev')
ax0_cj.set_yscale('log')
ax0_cj.invert_yaxis()

cf = ax0_jh.contourf(lat, lev, jh_wtem.T, levels=cf.levels, cmap=cmap, extend='both')
plt.colorbar(cf, ax=ax0_jh, location='right', label='wtem [mm/s]')
ax0_jh.set_title('PyTEMDiags')
ax0_jh.set_xlabel('lat')
ax0_jh.set_ylabel('lev')
ax0_jh.set_yscale('log')
ax0_jh.invert_yaxis()

levels = np.linspace(-1, 1, 11)
cf = ax0_diff.contourf(lat, lev, diff, levels=levels, cmap=cmap, extend='both', 
                       norm=mpl.colors.CenteredNorm())
plt.colorbar(cf, ax=ax0_diff, location='right', label='wtem [mm/s]')
ax0_diff.set_title('PyTEMDiags - CJ')
ax0_diff.set_xlabel('lat')
ax0_diff.set_ylabel('lev')
ax0_diff.set_yscale('log')
ax0_diff.invert_yaxis()


# ----------------- vertical residual pressure velocity -------------------
print('working on vertical residual presure velocity...')

jh_omtem = jhtem['omegatem'] * 1e3
cj_omtem = cjtem['wtem_jh'] * 1e3
diff = (jh_omtem.T - cj_omtem) / np.abs(np.maximum(jh_omtem.T, cj_omtem))
levels = np.linspace(-0.001, 0.001, 11) * 1e3
cmap='rainbow'

cf = ax1_cj.contourf(lat, lev, cj_omtem, levels=levels, cmap=cmap, extend='both')
plt.colorbar(cf, ax=ax1_cj, location='right', label='omega tem [1e-3 Pa/s]')
ax1_cj.set_title('CJ')
ax1_cj.set_xlabel('lat')
ax1_cj.set_ylabel('lev')
ax1_cj.set_yscale('log')
ax1_cj.invert_yaxis()

cf = ax1_jh.contourf(lat, lev, jh_omtem.T, levels=cf.levels, cmap=cmap, extend='both')
plt.colorbar(cf, ax=ax1_jh, location='right', label='omega tem [1e-3 Pa/s]')
ax1_jh.set_title('PyTEMDiags')
ax1_jh.set_xlabel('lat')
ax1_jh.set_ylabel('lev')
ax1_jh.set_yscale('log')
ax1_jh.invert_yaxis()

levels = np.linspace(-1, 1, 11)
cf = ax1_diff.contourf(lat, lev, diff, levels=levels, cmap=cmap, extend='both', norm=mpl.colors.CenteredNorm())
plt.colorbar(cf, ax=ax1_diff, location='right', label='omega tem difference [%]')
ax1_diff.set_title('PyTEMDiags - CJ')
ax1_diff.set_xlabel('lat')
ax1_diff.set_ylabel('lev')
ax1_diff.set_yscale('log')
ax1_diff.invert_yaxis()


# ----------------- zonal mean omega -------------------
print('working on zonal mean omega...')

cj_omb = cjtem['w_zm'] * 1e3
jh_omb = jhtem['wapb'] * 1e3
diff = (jh_omb.T - cj_omb) / np.maximum(np.abs(jh_omb.T), np.abs(cj_omb))
levels = np.arange(-0.2, 0.21, 0.02) * 2
cmap='rainbow'

cf = ax2_cj.contourf(lat, lev, cj_omb, levels=levels, cmap=cmap, extend='both')
plt.colorbar(cf, ax=ax2_cj, location='right', label='zm(omega) [1e-3 Pa/s]')
ax2_cj.set_title('CJ')
ax2_cj.set_xlabel('lat')
ax2_cj.set_ylabel('lev')
ax2_cj.set_yscale('log')
ax2_cj.invert_yaxis()

cf = ax2_jh.contourf(lat, lev, jh_omb.T, levels=cf.levels, cmap=cmap, extend='both')
plt.colorbar(cf, ax=ax2_jh, location='right', label='zm(omega) [1e-3 Pa/s]')
ax2_jh.set_title('PyTEMDiags')
ax2_jh.set_xlabel('lat')
ax2_jh.set_ylabel('lev')
ax2_jh.set_yscale('log')
ax2_jh.invert_yaxis()

levels = np.linspace(-1, 1, 11)
cf = ax2_diff.contourf(lat, lev, diff, levels=levels, cmap=cmap, extend='both', norm=mpl.colors.CenteredNorm())
plt.colorbar(cf, ax=ax2_diff, location='right', label='zm(omega) difference [%]')
ax2_diff.set_title('PyTEMDiags - CJ')
ax2_diff.set_xlabel('lat')
ax2_diff.set_ylabel('lev')
ax2_diff.set_yscale('log')
ax2_diff.invert_yaxis()


# ----------------- zonal v'theta' -------------------
print('working on v\'theta\'...')
cj_vptpb = cjtem['vtheta_zm']
jh_vptpb = jhtem['vptpb']
diff = (jh_vptpb.T - cj_vptpb) / np.maximum(np.abs(jh_vptpb.T), np.abs(cj_vptpb))
levels = 15
cmap='rainbow'

cf = ax3_cj.contourf(lat, lev, cj_vptpb, levels=levels, cmap=cmap, extend='both')
plt.colorbar(cf, ax=ax3_cj, location='right', label='zm(v\'theta\') [K m/s]')
ax3_cj.set_title('CJ')
ax3_cj.set_xlabel('lat')
ax3_cj.set_ylabel('lev')
ax3_cj.set_yscale('log')
ax3_cj.invert_yaxis()

cf = ax3_jh.contourf(lat, lev, jh_vptpb.T, levels=cf.levels, cmap=cmap, extend='both')
plt.colorbar(cf, ax=ax3_jh, location='right', label='zm(v\'theta\') [K m/s]')
ax3_jh.set_title('PyTEMDiags')
ax3_jh.set_xlabel('lat')
ax3_jh.set_ylabel('lev')
ax3_jh.set_yscale('log')
ax3_jh.invert_yaxis()

levels = np.linspace(-1, 1, 11)
cf = ax3_diff.contourf(lat, lev, diff, levels=levels, cmap=cmap, extend='both', norm=mpl.colors.CenteredNorm())
plt.colorbar(cf, ax=ax3_diff, location='right', label='zm(v\'theta\') difference [%]')
ax3_diff.set_title('PyTEMDiags - CJ')
ax3_diff.set_xlabel('lat')
ax3_diff.set_ylabel('lev')
ax3_diff.set_yscale('log')
ax3_diff.invert_yaxis()


# ----------------- d (zm(theta)) / dp -------------------
print('working on d (zm(theta)) / dp...')
cj_dtheta_dp = cjtem['dtheta_dp']
jh_dtheta_dp = jhtem['dthetab_dp']
diff = (jh_dtheta_dp.T - cj_dtheta_dp) / np.maximum(np.abs(jh_dtheta_dp.T), np.abs(cj_dtheta_dp))

cmap=plt.cm.rainbow
levels = [-0.01, -0.03, -0.1, -0.3, -1, -3, -10][::-1]
cmap = ListedColormap([cmap(i/len(levels)) for i in range(len(levels))])
norm = BoundaryNorm(levels, len(levels) - 1)

cf = ax4_cj.contourf(lat, lev, cj_dtheta_dp, levels=levels, cmap=cmap, extend='both', norm=norm)
plt.colorbar(cf, ax=ax4_cj, location='right', label='d (zm(theta)) / dp [K/Pa]', ticks=levels)
ax4_cj.set_title('CJ')
ax4_cj.set_xlabel('lat')
ax4_cj.set_ylabel('lev')
ax4_cj.set_yscale('log')
ax4_cj.invert_yaxis()

cf = ax4_jh.contourf(lat, lev, jh_dtheta_dp.T, levels=cf.levels, cmap=cmap, extend='both',norm=norm)
plt.colorbar(cf, ax=ax4_jh, location='right', label='d (zm(theta)) / dp [K/Pa]', ticks=levels)
ax4_jh.set_title('PyTEMDiags')
ax4_jh.set_xlabel('lat')
ax4_jh.set_ylabel('lev')
ax4_jh.set_yscale('log')
ax4_jh.invert_yaxis()

levels = np.linspace(-1, 1, 11)
cf = ax4_diff.contourf(lat, lev, diff, levels=levels, cmap=cmap, extend='both', norm=mpl.colors.CenteredNorm())
plt.colorbar(cf, ax=ax4_diff, location='right', label='d (zm(theta)) / dp difference [%]')
ax4_diff.set_title('PyTEMDiags - CJ')
ax4_diff.set_xlabel('lat')
ax4_diff.set_ylabel('lev')
ax4_diff.set_yscale('log')
ax4_diff.invert_yaxis()


# ----------------- psi -------------------
print('working on psi...')
cj_psi = cjtem['psi']
jh_psi = jhtem['psi']
diff = (jh_psi.T - cj_psi) / np.maximum(np.abs(jh_psi.T), np.abs(cj_psi))
levels = 15
cmap='rainbow'

cf = ax5_cj.contourf(lat, lev, cj_psi, levels=levels, cmap=cmap, extend='both')
plt.colorbar(cf, ax=ax5_cj, location='right', label='psi [Pa m/s]')
ax5_cj.set_title('CJ')
ax5_cj.set_xlabel('lat')
ax5_cj.set_ylabel('lev')
ax5_cj.set_yscale('log')
ax5_cj.invert_yaxis()

cf = ax5_jh.contourf(lat, lev, jh_psi.T, levels=cf.levels, cmap=cmap, extend='both')
plt.colorbar(cf, ax=ax5_jh, location='right', label='psi [Pa m/s]')
ax5_jh.set_title('PyTEMDiags')
ax5_jh.set_xlabel('lat')
ax5_jh.set_ylabel('lev')
ax5_jh.set_yscale('log')
ax5_jh.invert_yaxis()

levels = np.linspace(-1, 1, 11)
cf = ax5_diff.contourf(lat, lev, diff, levels=levels, cmap=cmap, extend='both', norm=mpl.colors.CenteredNorm())
plt.colorbar(cf, ax=ax5_diff, location='right', label='psi difference [%]')
ax5_diff.set_title('PyTEMDiags - CJ')
ax5_diff.set_xlabel('lat')
ax5_diff.set_ylabel('lev')
ax5_diff.set_yscale('log')
ax5_diff.invert_yaxis()


# ----------------- dpsicoslat_dlat -------------------
print('working on d(psicoslat)/dlat...')
cj_dpsicoslat_dlat = cjtem['dpsicoslat_dlat']
jh_dpsicoslat_dlat = jhtem['dpsicoslat_dlat']
diff = (jh_dpsicoslat_dlat.T - cj_dpsicoslat_dlat) / np.maximum(np.abs(jh_dpsicoslat_dlat.T), np.abs(cj_dpsicoslat_dlat))
levels = np.arange(-8000, 8001, 2000) / 4
cmap='rainbow'

cf = ax6_cj.contourf(lat, lev, cj_dpsicoslat_dlat, levels=levels, cmap=cmap, extend='both')
plt.colorbar(cf, ax=ax6_cj, location='right', label='d(psi cos(lat))/dlat [Pa/s]')
ax6_cj.set_title('CJ')
ax6_cj.set_xlabel('lat')
ax6_cj.set_ylabel('lev')
ax6_cj.set_yscale('log')
ax6_cj.invert_yaxis()

cf = ax6_jh.contourf(lat, lev, jh_dpsicoslat_dlat.T, levels=cf.levels, cmap=cmap, extend='both')
plt.colorbar(cf, ax=ax6_jh, location='right', label='d(psi cos(lat))/dlat [Pa/s]')
ax6_jh.set_title('PyTEMDiags')
ax6_jh.set_xlabel('lat')
ax6_jh.set_ylabel('lev')
ax6_jh.set_yscale('log')
ax6_jh.invert_yaxis()

levels = np.linspace(-1, 1, 11)
cf = ax6_diff.contourf(lat, lev, diff, levels=levels, cmap=cmap, extend='both', norm=mpl.colors.CenteredNorm())
plt.colorbar(cf, ax=ax6_diff, location='right', label='d(psi cos(lat))/dlat difference [%]')
ax6_diff.set_title('PyTEMDiags - CJ')
ax6_diff.set_xlabel('lat')
ax6_diff.set_ylabel('lev')
ax6_diff.set_yscale('log')
ax6_diff.invert_yaxis()


# ------------------------------------------------------


plt.tight_layout()
plt.savefig('compare.png', dpi=200)
plt.show()

# CONLUSIONS FOR NOW:
# OMEGA BAR IS FINE
# WTEM IS NOT FINE
# BECAUSE OMEGATEM IS NOT FINE
# BECAUSE d(ψ*cos(φ))/dφ IS NOT FINE
# BUT PSI MUST BE FINE BECAUSE VTEM IS FINE (MAYBE)


