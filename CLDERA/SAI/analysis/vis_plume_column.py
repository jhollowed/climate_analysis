import pdb
import numpy as np
from scipy import special
from scipy import integrate
import matplotlib.pyplot as plt
import metpy.constants as const
from metpy.units import units as u
import artist_utils as claut
import xarray as xr
from matplotlib import colors

params = {'text.usetex': False, 'mathtext.fontset': 'stixsans'}
plt.rcParams.update(params)

gmt = './cmaps/GMT_no_green.rgb'
cmap_gmt = claut.ncar_rgb_to_cmap(gmt)

# ---------- parameters ----------
lat0 = 15.15 * u.deg
lon0 = 120.35 * u.deg
mu = (22.59 * u.km).to(u.m)
sigma = (4 * u.km).to(u.m)
alpha = -2
zmin = (17 * u.km).to(u.m)
zmax = (30 * u.km).to(u.m)
tf = (9 * u.hr).to(u.s)
kso2 = 1/((25 * u.day).to(u.s))
ksulf = 1/((360*u.day).to(u.s))
kash = 1/((1*u.day).to(u.s))
Mso2 = 1.7e10 * u.kg
Mash = 5e10 * u.kg
dTstrat = 0.3 * u.K/u.day
dTsurf = 0.3 * u.K/u.day
w = 2.04
#w = 1

labels = ['SO2', 'Ash', 'Sulfate', 'Sulfate_nodecay']
kk = [kso2, kash, ksulf]
MM = [Mso2, Mash]

# vert distribution function
f = lambda z: (1/(np.sqrt(2*np.pi) * sigma)) * np.exp(-(z - mu)**2/(2*sigma**2)) *\
              (1 - special.erf(alpha * ((mu-z)/(np.sqrt(2)*sigma)).m))

#AA = 1.45416586e08 * u.kg*u.m/u.s
#VV = 0.000148156582 * 1/u.m
#t = np.linspace(0, 48*60*60, 1000) * u.s
#tmin = np.array([min(tt, tf).m for tt in t]) * t[0].u
#TESTA = (AA*VV/kk[0] * np.exp(-kk[0]*t) * (np.exp(kk[0] * tmin) - 1))
#TESTB = np.exp(-kk[0] * (t - t[500]))*TESTA[500]  +  (AA*VV/kk[0] * np.exp(-kk[0]*t) * (np.exp(kk[0] * tmin) - np.exp(kk[0] * min(tf, t[500]))))
#plt.ion()
#plt.plot(t, TESTA)
#plt.plot(t, TESTB)
#pdb.set_trace()v


def compute_tracer_evolution(nlev = 100, zk = None, dt_inject = None, dt_decay = None, tmax = 365*5*u.day, timestep_correction=False, analytic_only=False):
    '''
    Performs analytic and numerical solutions to SAI tracer tendency equations

    Parameters
    ----------
    nlev : int
        number of vertical levels from 12 to 35 km
    zk : Quantity array
        array of vertical positions in a length unit. If provided, nlev is ignored
    dt_inject : Quantity
        timestep during first 24 hours (injection period). Must have a unit
    dt_decay : Quantity
        timestep during decay (remaining 5 years). Must have a unit
    tmax : Quantity
        ending time of the decay period. Must have a unit. Defaults to 5 years
    timestep_correction : bool
        whether or not to enable timestep correction. If False, the correction is not performed 
        (numerical solution is done with the analytic tendencies). In either case, the method
        employed is a first-order explicit method (forward euler).
    analytic_only : bool
        whether or not to skip the numerical estimation and return only the analytic result. If
        True, the return for the numerical result is None
    '''
    # ---------- vertical distribution ----------
    if(zk is None):
        zk = (np.linspace(12, 35, nlev) * u.km).to(u.m)
    else:
        zk = zk.to(u.km)
    V = f(zk)
    V[np.logical_or(zk < zmin, zk > zmax)] = 0 * V[0].u
    sVk = np.sum(V)

    # ---------- normalization ----------
    A  = [MM[j] / (tf * sVk) for j in range(len(MM))]

    # ---------- tendencies ----------
    if(not timestep_correction):
        # analytic tendencies
        #dmdt = [lambda mk, tt, dt: -kso2*mk + A[0] * V * int(tt <= tf),
        #        lambda mk, tt, dt: -kash*mk + A[1] * V * int(tt <= tf),
        #        lambda mk, tt, mso2: -ksulf * mk + w*kso2*mso2, 
        #        lambda mk, tt, mso2: kso2*mso2]
        dmdt = [lambda mk, tt, dt: -kso2*mk + A[0] * V * int((tt+60*u.s) < tf),
                lambda mk, tt, dt: -kash*mk + A[1] * V * int((tt+60*u.s) < tf),
                lambda mk, tt, mso2: -ksulf * mk + w*kso2*mso2, 
                lambda mk, tt, mso2: kso2*mso2]
    else: 
        # correction "tendency" Δm/Δt for so2, ash, sulfur unchanged (same analytic tendency as above)
        dmdt = [lambda mk, tt, dt: ((np.exp(-kso2*dt)-1)*mk + A[0]*V/kso2 * np.exp(-kso2*(tt+dt)) * (np.exp(kso2*(min(tf,tt+dt))) - np.exp(kso2*min(tf,tt))))/dt,
                lambda mk, tt, dt: ((np.exp(-kash*dt)-1)*mk + A[1]*V/kash * np.exp(-kash*(tt+dt)) * (np.exp(kash*(min(tf,tt+dt))) - np.exp(kash*min(tf,tt))))/dt,
                lambda mk, tt, mso2: -ksulf * mk + w*kso2*mso2, 
                lambda mk, tt, mso2: kso2*mso2]

    # ---------- solve ----------

    
    # ----- numerical solutions

    # --- 2nd order RK (midpoint)
    #time_step = lambda dydt, mn, tn, dt: mn + dt * dydt(mn + 0.5*dt*dydt(mn, tn), tn + 0.5*dt)
    #time_step_sulf = lambda dydt, mn, tn, dt, mso2n: mn + dt * dydt(mn + 0.5*dt*dydt(mn, tn, mso2n), tn + 0.5*dt, mso2n)
    # --- 1st order Euler
    euler_step = lambda dydt, mn, tn, dt: mn + dt * dydt(mn, tn, dt)
    euler_step_sulf = lambda dydt, mn, tn, dt, mso2n: mn + dt * dydt(mn, tn, mso2n)

    # ----- numerical injection
    
    if(dt_inject is None):
        # default; do first day with fine resolution (~1.5 minute)
        dt_inject = (24/1000)*u.hr

    start_time = 0*u.s
    end_time = (1*u.day + dt_inject.to(u.day)).to(u.s)
    dt_inject = dt_inject.to(u.s)
    t_inject = np.arange(start_time.m, end_time.m, dt_inject.m) * u.s
    m_numerical_inject = np.zeros((4, len(t_inject), len(V))) * u.kg
    
    if(not analytic_only):
        print('solving numerical injection')
        for n in range(0, len(t_inject) - 1):
            print('{:.2f}%'.format(100*n/(len(t_inject)-1)), end='\r')
            m_numerical_inject[0, n+1, :] = euler_step(dmdt[0], m_numerical_inject[0, n], t_inject[n], dt_inject)
            m_numerical_inject[1, n+1, :] = euler_step(dmdt[1], m_numerical_inject[1, n], t_inject[n], dt_inject)
            m_numerical_inject[2, n+1, :] = euler_step_sulf(dmdt[2], m_numerical_inject[2, n], t_inject[n], dt_inject, m_numerical_inject[0, n])
            m_numerical_inject[3, n+1, :] = euler_step_sulf(dmdt[3], m_numerical_inject[3, n], t_inject[n], dt_inject, m_numerical_inject[0, n])

    # ----- numerical decay

    if(dt_decay is None):
        # default; decay with coarse resolution (~10 hr for tmax=5yr)
        dt_decay = (tmax/2000) # in days by default

    start_time = t_inject[-1]
    end_time = (tmax).to(u.s) + dt_decay.to(u.s)
    dt_decay = dt_decay.to(u.s)
    t_decay = np.arange(start_time.m, end_time.m, dt_decay.m) * u.s
    
    m_numerical_decay = np.zeros((4, len(t_decay), len(V))) * u.kg 
    m_numerical_decay[:,0,:] = m_numerical_inject[:,-1,:]

    if(not analytic_only):
        print('solving numerical decay')
        for n in range(0, len(t_decay) - 1):
            print('{:.2f}%'.format(100*n/(len(t_decay)-1)), end='\r')
            m_numerical_decay[0, n+1, :] = euler_step(dmdt[0], m_numerical_decay[0, n], t_decay[n], dt_decay)
            m_numerical_decay[1, n+1, :] = euler_step(dmdt[1], m_numerical_decay[1, n], t_decay[n], dt_decay)
            m_numerical_decay[2, n+1, :] = euler_step_sulf(dmdt[2], m_numerical_decay[2, n], t_decay[n], dt_decay, m_numerical_decay[0, n])
            m_numerical_decay[3, n+1, :] = euler_step_sulf(dmdt[3], m_numerical_decay[3, n], t_decay[n], dt_decay, m_numerical_decay[0, n])

    # ---- concat numerical inject, decay
    t = np.hstack([t_inject, t_decay])
    m_numerical = np.zeros((4, len(t), len(V))) * u.kg
    m_numerical[0] = np.concatenate([m_numerical_inject[0], m_numerical_decay[0]]) 
    m_numerical[1] = np.concatenate([m_numerical_inject[1], m_numerical_decay[1]]) 
    m_numerical[2] = np.concatenate([m_numerical_inject[2], m_numerical_decay[2]])
    m_numerical[3] = np.concatenate([m_numerical_inject[3], m_numerical_decay[3]])

    # ----- analytic solution
    m_analytic = np.zeros((4, len(t), len(V))) * u.kg

    tmin = np.array([min(tt, tf).m for tt in t]) * t[0].u
    m_analytic[0,:,:] = (A[0]/kk[0] * np.atleast_2d(V).T * np.exp(-kk[0]*t) * (np.exp(kk[0] * tmin) - 1)).T
    m_analytic[1,:,:] = (A[1]/kk[1] * np.atleast_2d(V).T * np.exp(-kk[1]*t) * (np.exp(kk[1] * tmin) - 1)).T
    m_analytic[2,:,:] = (w*A[0]/((kk[2] - kk[0])*kk[2]) * np.atleast_2d(V).T * \
                        np.exp(-kk[2]*t) * (\
                        kk[0]*(1 - np.exp(kk[2] * tmin)) - kk[2]*np.exp((kk[2]-kk[0])*t) * (1 - np.exp(kk[0] * tmin))\
                        )).T
    # VV only valid for w=1, ksulf=0
    m_analytic[3,:,:] = (A[0]/kk[0] * np.atleast_2d(V).T * np.exp(-kk[0]*t) * (1 - np.exp(kk[0]*tmin) + np.exp(kk[0]*t)*kk[0]*tmin)).T

    return m_analytic, m_numerical, zk, V, t
 




# --------- select plots to render 
main_figure = False
z_figure = False
t_figure = False
eva_figure = False
heating_figure = True

# ----- global plotting properties
label_fs =14
tick_fs = 12
grid_alph = 0.66

# ---------- solve!  ----------
if(main_figure or eva_figure):
    print('----- solving for main figure')
    m_analytic, m_numerical, zk, V, t = compute_tracer_evolution(100)
if(z_figure):
    print('----- solving for z figure')
    m_analytic_50, _, zk_50, V_50, t = compute_tracer_evolution(50, analytic_only=True)
    m_analytic_10, _, zk_10, V_10, _ = compute_tracer_evolution(10, analytic_only=True)
    m_analytic_5, _, zk_5, V_5, _ = compute_tracer_evolution(5, analytic_only=True)
if(t_figure):
    print('----- solving for t figure')
    # below we save the result to file so that we need not recompute it for plotting changes (takes a while)
    mfilename = './tmp/t_figure_sols_unit_kg_w{:.2f}.npz'.format(w)
    try:
        if(True): raise FileNotFoundError # switch to force re-computation and writing to file, even if exists
        mfile = np.load(mfilename)
        m_numerical_dt1, m_numerical_dt2, m_numerical_dt3, m_numerical_dt1_cor, m_numerical_dt2_cor, m_numerical_dt3_cor =\
                mfile['dt1']*u.kg, mfile['dt2']*u.kg, mfile['dt3']*u.kg, mfile['dt1c']*u.kg, mfile['dt2c']*u.kg, mfile['dt3c']*u.kg
        analytic_only = True
        print('read data from file')
    except FileNotFoundError:
        analytic_only = False

    nlev = 100
    tmax = 90 * u.day
    dt1 = 1800 * u.s
    dt2 = dt1/2
    dt3 = dt1*2

    # no correction
    m_analytic_dt1, _m_numerical_dt1, zk_dt1, V_dt1, t_dt1 = \
            compute_tracer_evolution(nlev, dt_inject=dt1, dt_decay=dt1, tmax=tmax, timestep_correction=False, analytic_only=analytic_only)
    m_analytic_dt2, _m_numerical_dt2, zk_dt2, V_dt2, t_dt2 = \
            compute_tracer_evolution(nlev, dt_inject=dt2, dt_decay=dt2, tmax=tmax, timestep_correction=False, analytic_only=analytic_only)
    m_analytic_dt3, _m_numerical_dt3, zk_dt3, V_dt3, t_dt3 = \
            compute_tracer_evolution(nlev, dt_inject=dt3, dt_decay=dt3, tmax=tmax, timestep_correction=False, analytic_only=analytic_only)
    # correction
    _, _m_numerical_dt1_cor, _, _, _ = \
            compute_tracer_evolution(nlev, dt_inject=dt1, dt_decay=dt1, tmax=tmax, timestep_correction=True, analytic_only=analytic_only)
    _, _m_numerical_dt2_cor, _, _, _ = \
            compute_tracer_evolution(nlev, dt_inject=dt2, dt_decay=dt2, tmax=tmax, timestep_correction=True, analytic_only=analytic_only)
    _, _m_numerical_dt3_cor, _, _, _ = \
            compute_tracer_evolution(nlev, dt_inject=dt3, dt_decay=dt3, tmax=tmax, timestep_correction=True, analytic_only=analytic_only)
    
    if(not analytic_only):
        print('saving data to file')
        m_numerical_dt1, m_numerical_dt2, m_numerical_dt3, m_numerical_dt1_cor, m_numerical_dt2_cor, m_numerical_dt3_cor = \
        _m_numerical_dt1, _m_numerical_dt2, _m_numerical_dt3, _m_numerical_dt1_cor, _m_numerical_dt2_cor, _m_numerical_dt3_cor
        np.savez(mfilename, dt1=m_numerical_dt1.m, dt2=m_numerical_dt2.m, dt3=m_numerical_dt3.m, 
                            dt1c=m_numerical_dt1_cor.m, dt2c=m_numerical_dt2_cor.m, dt3c=m_numerical_dt3_cor.m)
if(heating_figure):
    zk = xr.open_dataset('data/E3SM_Z3_HSyr6_horzAvg.nc')['Z3']
    zk = sorted(zk.values) * u(zk.units)
    m_analytic_heat, _, zk_heat, _, t_heat = compute_tracer_evolution(zk=zk, analytic_only=True)

    
# ==================================================================
# ====================== MAIN FIGURE ===============================
if(main_figure):
    print('----- plotting main figure')

    # ---- get densities in kg/m^3
    cell_volume = (200*u.km)**2 * np.diff(zk)[0]
    m_dens = m_analytic.to(u.kg) / cell_volume.to(u.m**3)
    m_dens[m_dens.m < 1e-12] = 1e-12 * (u.kg/u.m**3)
    m_dens = np.log10(m_dens.m)

    # ---- get total masses
    m_tot = np.sum(m_analytic.to('Mt'), axis=2)
    m_tot_num = np.sum(m_numerical.to('Mt'), axis=2)

    # ---------- vis ----------
    fig = plt.figure(figsize=(9,10))
    gs = fig.add_gridspec(nrows=4, ncols=4)
    ax0 = fig.add_subplot(gs[:,0])
    ax1 = fig.add_subplot(gs[0,1:])
    ax2 = fig.add_subplot(gs[1,1:])
    ax3 = fig.add_subplot(gs[2,1:])
    ax4 = fig.add_subplot(gs[3,1:])

    # ---- plotting properties
    clev = np.linspace(-7, -4, 10)
    clev_thin_contours = clev[clev % 1 != 0]
    clev_contours = np.linspace(-7, -4, 4)
    cm = plt.cm.RdYlBu_r
    cm = cmap_gmt
    cm = plt.cm.YlOrRd
    cbar_tick_labels = np.array(['{:.0f}'.format(lev) for lev in clev])
    cbar_tick_labels[clev % 1 != 0] = ''

    # ---- global ax settings
    for ax in [ax0, ax1, ax2, ax3, ax4]:
        if(ax != ax0):
            ax.yaxis.set_label_position("right")
            ax.yaxis.tick_right()
        ax.grid(ls='--', alpha=grid_alph)
        ax.set_ylabel('$z$ [km]', fontsize=label_fs)
        ax.tick_params(axis='both', which='major', labelsize=tick_fs)
        claut.format_ticks(ax)

    # ----- plot V(z) profile
    #  get higher res zk
    zk_fine = (np.linspace(12, 35, 1000) * u.km).to(u.m)
    V_fine = f(zk_fine)
    V_fine[np.logical_or(zk_fine < zmin, zk_fine > zmax)] = 0 * V_fine[0].u
    ax0.plot(V_fine.to(u.km**(-1)), zk_fine.to(u.km), '-k', lw=2.5)
    ax0.set_xlabel('$V(z)$ [km$^{-1}$]', fontsize=label_fs)
    ax0.set_ylim([12, 30])

    # ----- contour plots of denisty evolution
    cf1 = ax1.contourf(t.to(u.day), zk.to(u.km), m_dens[1].T, cmap=cm, levels=clev, extend='both', vmin=-7.5)
    ax1.contour(t.to(u.day), zk.to(u.km), m_dens[1].T, levels=clev_thin_contours, colors=['k']*len(clev_contours), linewidths=0.5)
    ax1.set_xlim([0, 14])
    ax1.set_xlabel('$t$ [days]', fontsize=label_fs)
    claut.add_annotation_box(ax1, 'ash', loc='upper right', fs=label_fs)

    ax2.contourf(t.to(u.month), zk.to(u.km), m_dens[0].T, cmap=cm, levels=clev, extend='both', vmin=-7.5)
    ax2.contour(t.to(u.day), zk.to(u.km), m_dens[0].T, levels=clev_thin_contours, colors=['k']*len(clev_contours), linewidths=0.5)
    ax2.set_xlim([0, 12])
    ax2.set_xlabel('$t$ [months]', fontsize=label_fs)
    claut.add_annotation_box(ax2, 'SO$_2$', loc='upper right', fs=label_fs)

    ax3.contourf(t.to(u.year), zk.to(u.km), m_dens[2].T, cmap=cm, levels=clev, extend='both', vmin=-7.5)
    ax3.contour(t.to(u.day), zk.to(u.km), m_dens[2].T, levels=clev_thin_contours, colors=['k']*len(clev_contours), linewidths=0.5)
    ax3.set_xlim([0, 5])
    ax3.set_xlabel('$t$ [years]', fontsize=label_fs)
    claut.add_annotation_box(ax3, 'sulfate', loc='upper right', fs=label_fs)

    # ----- color bar
    cbar = fig.colorbar(cf1, ax=ax1, orientation='horizontal', location='top')
    cbar.set_label('log$_{10}$( density [kg/m$^3$] )', fontsize=label_fs)
    cbar.set_ticklabels(cbar_tick_labels)
    cbar.ax.tick_params(labelsize=tick_fs)

    # ----- line plots of total mass
    #ax4.plot(t.to(u.month), m_tot[1], 'c', label='ash')
    ax4.plot(t.to(u.month), m_tot[0], 'c', label='SO$_2$', lw=4)
    ax4.plot(t.to(u.month), m_tot[2], 'orange', label='sulfate', lw=4)
    ax4.plot(t.to(u.month), m_tot[3], ':', color='orange', label='sulfate validation', lw=3)
    ax4.plot(t.to(u.month), m_tot_num[0], 'k', label='numerical solutions', lw=0.85)
    ax4.plot(t.to(u.month), m_tot_num[2], 'k', lw=0.85)
    ax4.plot(t.to(u.month), m_tot_num[3], 'k', lw=0.85)
    ax4.set_xlabel('$t$ [months]', fontsize=label_fs)
    ax4.set_ylabel('total mass [Mt]', fontsize=label_fs)
    ax4.legend(fancybox=False, fontsize=tick_fs, loc='lower right')
    ax4.set_xlim([-0.5, 12])

    plt.tight_layout()

    # manually add contour labels after tight layout
    cc3 = ax3.contour(t.to(u.day), zk.to(u.km), m_dens[2].T, levels=clev_contours, colors=['k']*len(clev_contours))
    ax3.clabel(cc3, cc3.levels, inline=True, fmt='%.0f', fontsize=tick_fs)
    cc1 = ax1.contour(t.to(u.day), zk.to(u.km), m_dens[1].T, levels=clev_contours, colors=['k']*len(clev_contours))
    ax1.clabel(cc1, cc1.levels, inline=True, fmt='%.0f', fontsize=tick_fs, manual=((0.78, 21.7), (2.9, 22.3), (4.71, 23.2), (6.65, 23.4)))
    cc2 = ax2.contour(t.to(u.day), zk.to(u.km), m_dens[0].T, levels=clev_contours, colors=['k']*len(clev_contours))
    ax2.clabel(cc2, cc2.levels, inline=True, fmt='%.0f', fontsize=tick_fs, manual=((0.94, 21.2), (2.67, 22.7), (4.81, 23.3)))

    plt.savefig('figs/sai_column_main_fig.png', dpi=300)



# ==================================================================
# ====================== VARIOUS NLEV  =============================
if(z_figure):
    print('----- plotting z figure')
    
    # ---------- vis ----------
    fig = plt.figure(figsize=(8,5))
    gs = fig.add_gridspec(nrows=1, ncols=2)
    axV = fig.add_subplot(gs[0,0])
    axM = fig.add_subplot(gs[0,1:])

    # ---- get total masses
    m_tot_50 = np.sum(m_analytic_50.to('Mt'), axis=2)
    m_tot_10 = np.sum(m_analytic_10.to('Mt'), axis=2)
    m_tot_5 = np.sum(m_analytic_5.to('Mt'), axis=2)

    # ---- global ax settings
    for ax in [axV, axM]:
        ax.grid(ls='--', alpha=grid_alph)
        ax.tick_params(axis='both', which='major', labelsize=tick_fs)
        claut.format_ticks(ax)

    # ----- plot V(z) profiles
    axV.plot(V_50.to(u.km**(-1)), zk_50.to(u.km), '-ok', lw=2, label='nlev=50')
    axV.plot(V_10.to(u.km**(-1)), zk_10.to(u.km), '-or', lw=2, label='nlev=10')
    axV.plot(V_5.to(u.km**(-1)), zk_5.to(u.km), '-om', lw=2, label='nlev=5')
    axV.legend(fancybox=False, fontsize=tick_fs, loc='upper right')
    axV.set_xlabel('$V(z)$ [km$^{-1}$]', fontsize=label_fs)
    axV.set_ylabel('$z$ [km]', fontsize=label_fs)
    axV.set_ylim([12, 30])

    # ----- line plots of total mass
    axM.plot(t.to(u.month), (m_tot_10[0]-m_tot_50[0]) * 1e14, '-r', lw=2, label='M(nlev=10) - M(nlev=50)')
    axM.plot(t.to(u.month), (m_tot_5[0]-m_tot_50[0]) * 1e14, '-m', lw=2, label='M(nlev=5) - M(nlev=50)')

    axM.set_xlabel('$t$ [months]', fontsize=label_fs)
    axM.set_ylabel('total SO$_2$ mass difference [$10^{-14}$ Mt]', fontsize=label_fs)
    axM.legend(fancybox=False, fontsize=tick_fs, loc='upper right')
    axM.set_xlim([-0.1, 4])
    axM.yaxis.set_label_position("right")
    axM.yaxis.tick_right()
    ticklabs = np.array(['{:.1f}'.format(t) for t in axM.get_yticks()])
    ticklabs[axM.get_yticks() % 0.5 != 0] = ''
    axM.set_yticklabels(ticklabs)

    plt.tight_layout()
    plt.savefig('figs/sai_column_nlev_fig.png', dpi=300)


# ==================================================================
# ====================== VARIOUS TIMESTEP  =========================
if(t_figure):
    print('----- plotting t figure')
    
    # ---------- vis ----------
    fig = plt.figure(figsize=(13,6.4))
    gs = fig.add_gridspec(nrows=3, ncols=2)
    axM = fig.add_subplot(gs[0:2,0])
    axF = fig.add_subplot(gs[2,0])
    axT = fig.add_subplot(gs[:,1])

    # ---- get total masses
    m_tot_dt1 = np.sum(m_analytic_dt1.to('Mt'), axis=2)
    m_tot_num_dt1 = np.sum(m_numerical_dt1.to('Mt'), axis=2)
    m_tot_num_dt1_cor = np.sum(m_numerical_dt1_cor.to('Mt'), axis=2)
    
    m_tot_dt2 = np.sum(m_analytic_dt2.to('Mt'), axis=2)
    m_tot_num_dt2 = np.sum(m_numerical_dt2.to('Mt'), axis=2)
    m_tot_num_dt2_cor = np.sum(m_numerical_dt2_cor.to('Mt'), axis=2)
    
    m_tot_dt3 = np.sum(m_analytic_dt3.to('Mt'), axis=2)
    m_tot_num_dt3 = np.sum(m_numerical_dt3.to('Mt'), axis=2)
    m_tot_num_dt3_cor = np.sum(m_numerical_dt3_cor.to('Mt'), axis=2)

    m_tot_num = [m_tot_num_dt1, m_tot_num_dt2, m_tot_num_dt3]
    m_tot_num_cor = [m_tot_num_dt1_cor, m_tot_num_dt2_cor, m_tot_num_dt3_cor]
    m_tot = [m_tot_dt1, m_tot_dt2, m_tot_dt3]

    t_dt = [t_dt1.to(u.day), t_dt2.to(u.day), t_dt3.to(u.day)]

    # fractional errors
    err_dt1 = (m_tot_num_dt1 - m_tot_dt1) / m_tot_dt1 * 100
    err_dt1_cor = (m_tot_num_dt1_cor - m_tot_dt1) / m_tot_dt1 * 100
    err_dt2 = (m_tot_num_dt2 - m_tot_dt2) / m_tot_dt2 * 100
    err_dt2_cor = (m_tot_num_dt2_cor - m_tot_dt2) / m_tot_dt2 * 100
    err_dt3 = (m_tot_num_dt3 - m_tot_dt3) / m_tot_dt3 * 100
    err_dt3_cor = (m_tot_num_dt3_cor - m_tot_dt3) / m_tot_dt3 * 100

    err = [err_dt1.m, err_dt2.m, err_dt3.m]
    err_cor = [err_dt1_cor.m, err_dt2_cor.m, err_dt3_cor.m]
    
    # ---- global ax settings
    tti = -1*u.day
    ttf = 90*u.day
    so2_lw = 1.2
    sulf_lw = 2.5
    for ax in [axM, axF, axT]:
        ax.grid(ls='--', alpha=grid_alph)
        ax.tick_params(axis='both', which='major', labelsize=tick_fs)
        claut.format_ticks(ax)
    colors = plt.cm.plasma((0.8, 0.5, 0.2))
    labels = ['$\Delta t$', '$\Delta t\>\2$', '$2\Delta t$']

    # ----- line plots of total mass
    for i in range(3):
        axM.plot(t_dt[i], m_tot_num[i][0], '-', lw=so2_lw, color=colors[i])
        axM.plot(t_dt[i], m_tot_num_cor[i][0], ':', lw=so2_lw, color=colors[i])
        # sulfate
        if(i == 0): sfx = ' = 1800 s'
        else: sfx = ''
        axM.plot(t_dt[i], m_tot_num[i][2], '-', lw=sulf_lw, color=colors[i], label='$\Delta t${}'.format(sfx))
        axM.plot(t_dt[i], m_tot_num_cor[i][2], ':', lw=sulf_lw, color=colors[i], label='$\Delta t$, corrected')
    # analytic from the finest solution
    axM.plot(t_dt[1], m_tot[1][0], '-k', lw=2, label='analytic')
    axM.plot(t_dt[1], m_tot[1][2], '-k', lw=2)
    
    #axM.set_ylim([0, 20])
    axM.set_xlim([tti, ttf])
    
    axM.set_ylabel('total SO$_2$ mass [Mt]', fontsize=label_fs)
    axM.set_xlabel('$t$ [hours]', fontsize=label_fs)
    leg = axM.legend(fancybox=False, fontsize=tick_fs, ncol=2, loc='center right')
    bb = leg.get_bbox_to_anchor().inverse_transformed(axM.transAxes)
    # Change to location of the legend.
    yOffset = -0.08
    bb.y0 += yOffset
    bb.y1 += yOffset
    leg.set_bbox_to_anchor(bb, transform = axM.transAxes)
    
    # ----- text labels of line groups
    axM.text(20, 8, 'SO$_2$', rotation=-20, fontsize=label_fs)
    axM.text(30, 27, 'sulfate', rotation=25, fontsize=label_fs)
    
    # ----- line plots of total mass fractional error
    for i in range(3):
        axF.plot(t_dt[i], err[i][0], '-', lw=so2_lw, color=colors[i])
        axF.plot(t_dt[i], err_cor[i][0], ':', lw=so2_lw, color=colors[i])
        # sulfate
        axF.plot(t_dt[i], err[i][2], '-', lw=sulf_lw, color=colors[i])
        axF.plot(t_dt[i], err_cor[i][2], ':', lw=sulf_lw, color=colors[i])
    
    axF.set_yticks([-10, -5, 0, 5, 10, 15])
    axF.set_ylim([-7, 13])
    axF.set_xlim([tti/6, ttf/9])
    
    axF.set_xlabel('$t$ [hours]', fontsize=label_fs)
    axF.set_ylabel('total SO$_2$\nmass error [%]', fontsize=label_fs)
    
    # ----- line plots of peak mass fractional error vs dt
    dtt = np.array([dt1.m, dt2.m, dt3.m])*u.s
    peak_mass = [np.max(m_tot[i], axis=1) for i in range(3)]
    peak_mass_num = [np.max(m_tot_num[i], axis=1) for i in range(3)]
    peak_mass_num_cor = [np.max(m_tot_num_cor[i], axis=1) for i in range(3)]
    peak_mass_err = np.array([np.abs(peak_mass[i] - peak_mass_num[i])/peak_mass[i] * 100 for i in range(3)]).T
    peak_mass_err_cor = np.array([np.abs(peak_mass[i] - peak_mass_num_cor[i])/peak_mass[i] * 100 for i in range(3)]).T

    # maximum absolute value of relative err
    #max_err_so2 = [np.nanmax(ee[0] for ee in err]
    #max_err_cor_so2 = [np.nanmax(ee[0]) for ee in err_cor]
    
    # absolute value of relative err at 90 hr
    #max_err_so2 = [ee[0][-1] for ee in err]
    #max_err_cor_so2 = [ee[0][-1] for ee in err_cor]
    axT.plot(dtt, peak_mass_err[2], '-o', lw=sulf_lw+2, color='grey', label='sulfate', fillstyle='none', markeredgecolor='r')
    axT.plot(dtt, peak_mass_err_cor[2], ':o', lw=sulf_lw+2, color='grey', label='sulfate, corrected', fillstyle='none', markeredgecolor='r')
    axT.plot(dtt, peak_mass_err[0], '-d', lw=so2_lw, color='k', label='SO$_2$', fillstyle='none', markeredgecolor='r')
    axT.plot(dtt, peak_mass_err_cor[0], ':d', lw=so2_lw, color='k', label='SO$_2$, corrected', fillstyle='none', markeredgecolor='r')
    axT.plot(dtt, peak_mass_err[1], '-s', lw=sulf_lw, color='k', label='ash', fillstyle='none', markeredgecolor='r')
    axT.plot(dtt, peak_mass_err_cor[1], ':s', lw=sulf_lw, color='k', label='ash, corrected', fillstyle='none', markeredgecolor='r')
    
    axT.set_xticks(np.arange(0, 1800*3, 1800/2))
    axT.set_xlim([500, 4000])
    axT.set_xlabel('$\Delta t$ [seconds]', fontsize=label_fs)
    axT.set_ylabel('error in peak tracer mass [%]', fontsize=label_fs)
    axT.yaxis.set_label_position("right")
    axT.yaxis.tick_right()
    axT.legend(fancybox=False, fontsize=tick_fs, ncol=1, loc='upper left')


    plt.tight_layout()
    plt.savefig('figs/sai_column_dt_fig.png', dpi=300)


# ==================================================================
# ====================== EVA COMPARISON  ===========================
if(eva_figure):
    print('----- plotting EVA figure')
    
    # ---------- vis ----------
    fig = plt.figure(figsize=(6,5))
    gs = fig.add_gridspec(nrows=1, ncols=1)
    ax = fig.add_subplot(gs[0,0])

    # --- molar masses, eva params
    mm_s = 32.07
    mm_so2 = 64.066
    mm_h2so4 = 98.1
    facid = 0.75
    mso20 = 9*u.Mt
    tprod = 180*u.day
    tprod_mod = 1/kso2
    tloss = 330*u.day
    tloss_mod = 1/ksulf

    # ---- get total masses
    m_tot = np.sum(m_analytic.to('Mt'), axis=2)
    m_tot_eva = (mm_h2so4/facid/mm_so2) * (mm_so2/mm_s) * mso20 * (1 - np.exp(-t/tprod)) * np.exp(-t/tloss)
    m_tot_eva_mod = (mm_h2so4/facid/mm_so2) * (mm_so2/mm_s) * mso20 * (1 - np.exp(-t/tprod_mod)) * np.exp(-t/tloss)

    # ---- global ax settings
    ax.grid(ls='--', alpha=0.3, color='k')
    ax.tick_params(axis='both', which='major', labelsize=tick_fs)
    claut.format_ticks(ax)

    # ----- line plots of total mass
    ax.plot(t.to(u.year), m_tot[2], '-r', lw=2, label='this work')
    ax.fill_between(t.to(u.year), m_tot[2]*1.25, m_tot[2]*0.75, color='r', alpha=0.5)
    ax.plot(t.to(u.year), m_tot_eva, '-k', lw=2, label='EVA v1.0, $1/k_\\mathrm{SO2} = 180$ days')
    ax.fill_between(t.to(u.year), m_tot_eva*1.25, m_tot_eva*0.75, color='k', alpha=0.3)
    ax.plot(t.to(u.year), m_tot_eva_mod, '--k', lw=2, label='EVA v1.0, $1/k_\\mathrm{SO2} = 30$ days')
    ax.fill_between(t.to(u.year), m_tot_eva_mod*1.25, m_tot_eva_mod*0.75, color='k', alpha=0.3)
    
    ax.set_xlim([-0.1, 3])
    ax.set_xlabel('$t$ [years]', fontsize=label_fs)
    ax.set_ylabel('total sulfate mass [Mt]', fontsize=label_fs)
    ax.legend(fancybox=False, fontsize=tick_fs, loc='upper right')
    #axM.set_xlim([-0.1, 4])

    plt.tight_layout()
    plt.savefig('figs/sai_column_eva_fig.png', dpi=300)



# ==================================================================
# ====================== HEATING VIS ===============================
if(heating_figure):

    # ---------- vis ----------
    fig = plt.figure(figsize=(10,5))
    gs = fig.add_gridspec(nrows=1, ncols=1)
    ax1 = fig.add_subplot(gs[0,0])
    axes = [ax1]
    #ax2 = fig.add_subplot(gs[1,0])
    #axes = [ax1, ax2]

    # --- heating parameters
    dTstrat = 0.35 * u.K/u.day
    qstar = 2e-4 * u.kg/u.m**3
    dTsurf = -0.012 * u.K/u.day
    taustar = 5.883e10 * u.kg
    maxaod_z = 5 * u.km
    maxaod_idx = np.searchsorted(zk_heat, maxaod_z)
    # four entries for so2, ash, sulfate, and sulfate_contol
    # turn off sulfate_control contribution via toggle at index 3 V
    b = np.array([1, 1, 1, 0]) * 1/u.kg
    strat_heat_toggle = [1, 1, 1, 0]
    aod_toggle = [1, 1, 1, 0]
    norm_contribution = np.array([1, 0, 1, 0])
    
    # ---- get densities in kg/m^3
    cell_volume = (200*u.km)**2 * np.abs(np.diff(zk))
    m_dens = np.zeros(m_analytic_heat.shape) * u.kg/(u.m**3)
    for k in range(len(zk_heat)-1):
        m_dens[:,:,k] = m_analytic_heat.to(u.kg)[:,:,k] / cell_volume.to(u.m**3)[k]
    #zk_mid = (zk_heat[:-1] + np.diff(zk)/2)

    # --- copute strat. heating
    # update heating normalization, ignoring ash
    qstar = np.max(np.sum(m_dens * norm_contribution[:,np.newaxis,np.newaxis], axis=0))
    print('qstar = {:.2e}'.format(qstar))
    strat_heating = (m_dens / qstar) * dTstrat
    for j in range(len(strat_heating)):
        strat_heating[j] = strat_heating[j] * strat_heat_toggle[j]
    
    # --- compute AOD
    aod_cooling = np.zeros(m_analytic_heat.shape) * u.K/u.day
    aod_cooling2 = np.zeros(m_analytic_heat.shape) * u.K/u.day
    tau = np.flip(np.flip(m_analytic_heat.to('kg'), 2).cumsum(axis=2), axis=2) * b[:, np.newaxis, np.newaxis] 
    for j in range(len(tau)):
        tau[j] = tau[j] * aod_toggle[j]
    # update cooling normalization, ignoring ash
    taustar = np.max(np.sum(tau*norm_contribution[:,np.newaxis,np.newaxis], axis=0))
    print('taustar = {:.2e}'.format(taustar))
    aod_cooling[:,:,0:maxaod_idx+1] = tau[:,:,0:maxaod_idx+1]/taustar * dTsurf 
    aod_cooling2 = tau/taustar * dTsurf

    # --- total heating rate
    dTdt1 = np.sum(strat_heating + aod_cooling, axis=0)
    dTdt2 = np.sum(strat_heating + aod_cooling2, axis=0)
    #dTdt = [dTdt1, dTdt2]
    dTdt = [dTdt1]
    
    for i in range(len(axes)):
        axi = axes[i]

        # ---- global ax settings
        axi.grid(ls=':', alpha=0.3, color='k')
        axi.tick_params(axis='both', which='major', labelsize=tick_fs)
        claut.format_ticks(axi)
        cm = plt.cm.RdYlBu_r
        clev = np.linspace(-0.012, 0, 7)
        clev = np.hstack([clev, np.linspace(0.05, 0.35, 7)])
        clev_lab = [-0.01, 0, 0.01, 0.05, 0.1, 0.2, 0.3]
        divnorm=colors.TwoSlopeNorm(vcenter=0.)
        bg = np.array([245,251,218])/255
        axi.set_facecolor(bg)
        
        # ----- contour plots of heating rate
        cf = axi.contourf(t_heat.to(u.month), zk_heat.to(u.km), dTdt[i].T, cmap=cm, levels=clev, extend='both', norm=divnorm)
        cfl = axi.contour(t_heat.to(u.month), zk_heat.to(u.km), dTdt[i].T, levels=clev_lab, colors = 'k', linewidths=1)
        
        # ----- color bar
        cbar = fig.colorbar(cf, ax=axi, orientation='vertical', location='right')
        cbar.set_label('dT/dt [K/day]', fontsize=label_fs)
        cbar.ax.tick_params(labelsize=tick_fs)
    
        plt.tight_layout()
        
        # ----- contour labels
        xx,yy = 4.31, 31.7
        axi.clabel(cfl, cfl.levels, inline=True, fmt='%.2f', fontsize=tick_fs, inline_spacing=0, 
                   manual = ((2.33, 5.0), (4.13, 30), (4.13, 26.0), (6.56, 23.9), (8.77, 22.1), (5.23, 20.8), (2.44, 20.7)) )
        
        axi.set_xlabel('$t$ [months]', fontsize=label_fs)
        axi.set_ylabel('$z$ [km]', fontsize=label_fs)
        axi.set_ylim([0, 35])
        axi.set_xlim([-1, 12])

    plt.savefig('figs/sai_column_heating_fig.png', dpi=300)


plt.show()





