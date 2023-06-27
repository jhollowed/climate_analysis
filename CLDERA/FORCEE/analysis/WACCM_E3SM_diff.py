import sys
import pdb
import glob
import scipy
import numpy as np
import xarray as xr
from matplotlib import colors
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

TICK_FS = 11
LABEL_FS = 14
plt.rc('font', size=TICK_FS)          # controls default text sizes
plt.rc('axes', titlesize=LABEL_FS)     # fontsize of the axes title
plt.rc('axes', labelsize=LABEL_FS)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=TICK_FS)    # fontsize of the tick labels
plt.rc('ytick', labelsize=TICK_FS)    # fontsize of the tick labels
plt.rc('legend', fontsize=TICK_FS)    # legend fontsize
plt.rc('figure', titlesize=LABEL_FS)  # fontsize of the figure title

# --------------------------------------------------------------

# ------ monthly zonally-meaned data
pdir         = '/pscratch/sd/j/jhollo/E3SM/historical/CLDERA_2ndHistorical_1950-2014_monthly'
E90_E3SM     = '{}/histoircal_h0_E90j_1850-2014.regrid.fv0.9x1.25_bilinear.nc'.format(pdir)
zm_E90_E3SM  = '{}/histoircal_h0_E90j_1850-2014.regrid.fv0.9x1.25_bilinear.zonalmean.nc'.format(pdir)
zm_ST80_E3SM = '{}/histoircal_h0_ST80_25j_1850-2014.regrid.fv0.9x1.25_bilinear.zonalmean.nc'.format(pdir)
zm_AOA_E3SM  = '{}/histoircal_h0_AOA_1850-2014.regrid.fv0.9x1.25_bilinear.zonalmean.nc'.format(pdir)
zm_PS_E3SM   = '{}/histoircal_h0_PS_1850-2014.regrid.fv0.9x1.25_bilinear.zonalmean.nc'.format(pdir)

pdir_WACCM   = '/pscratch/sd/j/jhollo/E3SM/historical/WACCM_CCMI_2022'
zm_WACCM     = '{}/f.e21.FWHIST.f09_f09_mg17.cesm2.1.4.REFD1.ens_mean.cam.h0zm.196001-201812.nc'.format(pdir_WACCM)
#zm_AOA_WACCM = '{}/AoA_waccm6_refd1.04_AOA1mf_1970-2019_ba_0_100.0_ck_0_50.0.nc'.format(pdir_WACCM)
#AOA_NAME='AOA'
zm_AOA_WACCM     = '{}/f.e21.FWHIST.f09_f09_mg17.cesm2.1.4.REFD1.H4.cam.h0zm.196001-201812.AOA1mf.nc'.format(pdir_WACCM)
AOA_NAME='AOA1mf'



# ------ read data
print('reading data...')
dat_e3sm   = xr.open_dataset(E90_E3SM) 
e90_e3sm   = xr.open_dataset(zm_E90_E3SM)['E90j']      # in ppb
st80_e3sm  = xr.open_dataset(zm_ST80_E3SM)['ST80_25j'] # in ppb
aoa_e3sm   = xr.open_dataset(zm_AOA_E3SM)['AOA']       # in years
ps_e3sm    = xr.open_dataset(zm_PS_E3SM)['PS']         # in years
p0_e3sm    = dat_e3sm['P0']
hy_e3sm    = [dat_e3sm['hyam'], dat_e3sm['hyai'], 
              dat_e3sm['hybm'], dat_e3sm['hybi']]

dat_waccm  = xr.open_dataset(zm_WACCM)
e90_waccm  = dat_waccm['E90'] * 1e9                    # in ppb
st80_waccm = dat_waccm['ST80_25'] * 1e9                # in ppb
aoa_waccm  = xr.open_dataset(zm_AOA_WACCM)[AOA_NAME]  # in years?
#aoa1_waccm = dat_waccm['AOA1']                         # in years?
#aoa2_waccm = dat_waccm['AOA2']                         # in years?
ps_waccm   = dat_waccm['PS']                           # in Pa
p0_waccm   = dat_waccm['P0']
hy_waccm   = [dat_waccm['hyam'], dat_waccm['hyai'], 
              dat_waccm['hybm'], dat_waccm['hybi']]
pdb.set_trace()

# ------ take time means
print('taking time means...')
e90_e3sm   = e90_e3sm.mean('time')
st80_e3sm  = st80_e3sm.mean('time')
aoa_e3sm   = aoa_e3sm.mean('time')
ps_e3sm    = ps_e3sm.mean('time')

e90_waccm  = e90_waccm.mean('time')
st80_waccm = st80_waccm.mean('time')
aoa_waccm  = aoa_waccm.mean('time')
aoa1_waccm = aoa1_waccm.mean('time')
aoa2_waccm = aoa2_waccm.mean('time')
ps_waccm   = ps_waccm.mean('time')

# ------ get coords
print('constructing P...')
lat_e3sm = dat_e3sm['lat']
lev_e3sm = dat_e3sm['lev']
LAT_e3sm, LEV_e3sm = np.meshgrid(lat_e3sm, lev_e3sm)
p_e3sm = (hy_e3sm[0] * p0_e3sm + hy_e3sm[2] * ps_e3sm) / 100
points_e3sm = []
for k in range(len(lev_e3sm)):
    for i in range(len(lat_e3sm)):
        points_e3sm.append([lat_e3sm.values[i], p_e3sm.values[k][i]])
points_e3sm = np.array(points_e3sm)

lat_waccm = dat_waccm['lat']
lev_waccm = dat_waccm['lev']
LAT_waccm, LEV_waccm = np.meshgrid(lat_waccm, lev_waccm)
p_waccm = (hy_waccm[0] * p0_waccm + hy_waccm[2] * ps_waccm) / 100
points_waccm = []
for k in range(len(lev_waccm)):
    for i in range(len(lat_waccm)):
        points_waccm.append([lat_waccm.values[i], p_waccm.values[k][i]])
points_waccm = np.array(points_waccm)

# ------ verify grids
verify_grid = False
if verify_grid:
    print('plotting grid verification...')
    fig = plt.figure(figsize=(12, 12))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)

    ax1.plot([0,0], [0,0], '-k', label='E3SM model levels', lw=0.7)
    ax1.plot(LAT_e3sm.T, LEV_e3sm.T, '-k', lw=0.5)
    ax1.plot([0,0], [0,0], ':k', label='WACCM model levels', lw=0.7)
    ax1.plot(LAT_waccm.T, LEV_waccm.T, '--k', lw=0.5)
    ax1.plot(points_e3sm.T[0], points_e3sm.T[1], '.r', 
             label='E3SM pressure positions', ms=0.85)
    ax1.plot(points_waccm.T[0], points_waccm.T[1], '.b', 
             label='WACCM pressure positions', ms=0.85)
    ax1.set_ylabel('p [hPa]')
    ax1.set_xlabel('lat [deg]')
    ax1.legend(loc='upper left', framealpha=1, fancybox=False)
    ax1.set_title('surface')

    ax2.plot(LAT_e3sm.T, LEV_e3sm.T, '-k', lw=0.7)
    ax2.plot(LAT_waccm.T, LEV_waccm.T, '--k', lw=0.7)
    ax2.plot(points_e3sm.T[0], points_e3sm.T[1], '.r', ms=0.85)
    ax2.plot(points_waccm.T[0], points_waccm.T[1], '.b', ms=0.85)
    ax2.set_ylabel('p [hPa]')
    ax2.set_xlabel('lat [deg]')
    ax2.set_yscale('log')
    ax2.set_title('stratosphere')
    ax2.yaxis.set_label_position("right")
    ax2.yaxis.tick_right()

    ax1.invert_yaxis()
    ax2.invert_yaxis()

    plt.tight_layout()
    plt.show()

# ------ interpolating E3SM data
print('interpolating E3SM data to WACCM grid...')
METHOD='cubic'
e90_e3sm  = scipy.interpolate.griddata(points_e3sm, values=np.ravel(e90_e3sm.values), 
                                        xi=points_waccm, method=METHOD)
st80_e3sm = scipy.interpolate.griddata(points_e3sm, values=np.ravel(st80_e3sm.values), 
                                        xi=points_waccm, method=METHOD)
aoa_e3sm  = scipy.interpolate.griddata(points_e3sm, values=np.ravel(aoa_e3sm.values), 
                                        xi=points_waccm, method=METHOD)

# ------ clean difference (remove nans outside of E3SM convex hull)
print('selecting data intersection...')
mask = np.isnan(np.ravel(e90_e3sm))
e90_e3sm   = np.ravel(e90_e3sm)[~mask]
st80_e3sm  = np.ravel(st80_e3sm)[~mask]
aoa_e3sm   = np.ravel(aoa_e3sm)[~mask]
e90_waccm  = np.ravel(e90_waccm)[~mask]
st80_waccm = np.ravel(st80_waccm)[~mask]
aoa_waccm  = np.ravel(aoa_waccm)[~mask] # TODO: also test aoa1, aoa2?
LAT = points_waccm.T[0]
LEV = points_waccm.T[1]
LAT = LAT[~mask]
LEV = LEV[~mask]

# ------ take diff
print('taking E3SM - WACCM difference...')
e90_diff  = e90_e3sm - e90_waccm
st80_diff = st80_e3sm - st80_waccm
aoa_diff  = aoa_e3sm - aoa_waccm


# ----------------------------------------------
# ------------------ PLOTTING ------------------

print('creating figure object...')
fig = plt.figure(figsize=(14, 15))
ax1_e3sm  = fig.add_subplot(331)
ax2_e3sm  = fig.add_subplot(332)
ax3_e3sm  = fig.add_subplot(333)
ax1_waccm = fig.add_subplot(334)
ax2_waccm = fig.add_subplot(335)
ax3_waccm = fig.add_subplot(336)
ax1_diff  = fig.add_subplot(337)
ax2_diff  = fig.add_subplot(338)
ax3_diff  = fig.add_subplot(339)

ax_e90   = [ax1_e3sm, ax1_waccm, ax1_diff]
e90_all  = [e90_e3sm, e90_waccm, np.abs(e90_diff)]
ax_st80  = [ax2_e3sm, ax2_waccm, ax2_diff]
st80_all = [st80_e3sm, st80_waccm, np.abs(st80_diff)]
ax_aoa   = [ax3_e3sm, ax3_waccm, ax3_diff]
aoa_all  = [aoa_e3sm, aoa_waccm, np.abs(aoa_diff)]
dat_labels = ['E3SM', 'WACCM', '|E3SM - WACCM|']


# ------------------- E90
print('plotting E90...')

for i in range(3):

    levels = [0.1, 1, 10, 20, 30, 40, 50, 60, 70, 
              80, 90, 100, 110, 120, 130, 140, 150]
    lablevels = [1, 10, 50, 100, 110, 130, 150]
    divnorm=colors.TwoSlopeNorm(vmin=0, vcenter=50, vmax=150)
    if(i == 2):
        # is difference
        levels = [0.1, 1, 2, 4, 6, 8, 10]
        lab_levels = [2, 4, 6, 8, 10, 12]
        divnorm = None

    ax = ax_e90[i]
    e90 = e90_all[i]
    dat_label = dat_labels[i]
    
    cf = ax.tricontourf(LAT, LEV, e90, levels=levels, cmap=plt.cm.RdYlBu_r, 
                         norm=divnorm, extend='both')
    ax.tricontour(LAT, LEV, e90, levels=levels, colors=['k'], linewidths=[0.4])          # draw contours
    cfl = ax.tricontour(LAT, LEV, e90, levels=lablevels, colors=['k'], linewidths=[0.4]) # contours to label
    if(i != 3):
        # is not difference
        ax.tricontour(LAT, LEV, e90, levels=[90], colors=['k'], linewidths=[2.5])            # bold 90 ppb
        cflm01 = ax.tricontour(LAT, LEV, e90, levels=[0.1], colors=['k'], linewidths=[0.4])  # label 0.1 ppb
        cflm005 = ax.tricontour(LAT, LEV, e90, levels=[0.05], colors=['k'], linewidths=[0.4])  # label 0.1 ppb
        cflm001 = ax.tricontour(LAT, LEV, e90, levels=[0.01], colors=['k'], linewidths=[0.4])  # label 0.1 ppb

    ax.clabel(cfl, inline=1, fontsize=TICK_FS, fmt='%1.0f')
    ax.clabel(cflm01, inline=1, fontsize=TICK_FS, fmt='%1.1f')
    ax.clabel(cflm005, inline=1, fontsize=TICK_FS, fmt='%1.2f')
    ax.clabel(cflm001, inline=1, fontsize=TICK_FS, fmt='%1.2f')
    if(i != 1):
        cb  = fig.colorbar(cf, ax=ax, location='top')
        if(i == 0): 
            cb.set_label('E90 [ppb]')
            cbt = [0, 50, 100, 150]
            cb.set_ticks(cbt)
            cb.set_ticklabels([str(t) for t in cbt])

    ax.set_ylim([10, 700])
    ax.set_yscale('log')
    ax.invert_yaxis()
    ax.set_ylabel('{}\np [hPa]'.format(dat_label))
    if(i==2):
        ax.set_xlabel('lat [deg]')
    ax.yaxis.set_major_formatter(ScalarFormatter())
    ax.yaxis.set_minor_formatter(ScalarFormatter())
    for label in ax.get_yticklabels(minor=True)[::2]:
        label.set_visible(False)


# ------------------- ST80
print('plotting ST80...')

for i in range(3):
    
    ax = ax_st80[i]
    st80 = st80_all[i]
    dat_label = dat_labels[i]
    
    levels=[20, 40, 60, 80, 100, 120, 140, 160, 180, 199]
    lablevels=[20, 60, 100, 140, 180, 199]
    #if(i == 2):
    #    # is difference
    #    levels = len(levels)
    #    lab_levels = levels
    
    cf = ax.tricontourf(LAT, LEV, st80, levels=levels, cmap=plt.cm.YlOrRd, extend='both')
    ax.tricontour(LAT, LEV, st80, levels=levels, colors=['k'], linewidths=[0.4])          # draw contours
    cfl = ax.tricontour(LAT, LEV, st80, levels=lablevels, colors=['k'], linewidths=[0.4]) # contours to label

    ax.clabel(cfl, inline=1, fontsize=TICK_FS, fmt='%1.0f')
    if(i != 1):
        cb  = fig.colorbar(cf, ax=ax, location='top')
        if(i == 0): cb.set_label('ST80_25 [ppb]')

    ax.set_ylim([50, 400])
    ax.set_yscale('log')
    ax.invert_yaxis()
    #ax.set_ylabel('lev [hPa]')
    if(i==2):
        ax.set_xlabel('lat [deg]')
    ax.yaxis.set_major_formatter(ScalarFormatter())
    ax.yaxis.set_minor_formatter(ScalarFormatter())


# ------------------- AOA
print('plotting AOA...')

for i in range(3):
    
    ax = ax_aoa[i]
    aoa = aoa_all[i]
    dat_label = dat_labels[i]
    
    levels=[0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5]
    lablevels=[1, 2, 3, 4, 5]
    if(i == 2):
        # is difference
        levels = [0.25, 0.5, 0.75, 1.0, 1.25, 1.5]
        lab_levels = [0.5, 1.0, 1.5, 2.0]
    
    cf = ax.tricontourf(LAT, LEV, aoa, levels=levels, cmap=plt.cm.YlGnBu_r, extend='both')
    ax.tricontour(LAT, LEV, aoa, levels=levels, colors=['k'], linewidths=[0.4])          # draw contours
    cfl = ax.tricontour(LAT, LEV, aoa, levels=lablevels, colors=['k'], linewidths=[0.4]) # contours to label

    ax.clabel(cfl, inline=1, fontsize=TICK_FS, fmt='%1.0f')
    if(i != 1):
        cb  = fig.colorbar(cf, ax=ax, location='top')
        if(i == 0): cb.set_label('AOA [years]')

    ax.set_ylim([1, 100])
    ax.set_yscale('log')
    ax.yaxis.set_label_position("right")
    ax.yaxis.tick_right()
    ax.invert_yaxis()
    ax.set_ylabel('lev [hPa]')
    if(i==2):
        ax.set_xlabel('lat [deg]')
    ax.yaxis.set_major_formatter(ScalarFormatter())
    ax.yaxis.set_minor_formatter(ScalarFormatter())
    for label in ax.get_yticklabels(minor=True)[1::2]:
        label.set_visible(False)

plt.tight_layout(pad=3.0)
plt.show()
