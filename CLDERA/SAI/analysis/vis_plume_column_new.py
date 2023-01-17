import pdb
import metpy
import numpy as np
from ambiance import Atmosphere as stdatm
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
mu = (17 * u.km).to(u.m)
sigma = (1.5 * u.km).to(u.m)
t0 = (30*u.day).to(u.s)
tf = t0 + (24 * u.hr).to(u.s)
dt = tf-t0
kso2 = 1/((25 * u.day).to(u.s))
ksulf = 1/((360*u.day).to(u.s))
kash = 1/((1*u.day).to(u.s))
Mso2 = 1.7e10 * u.kg
Mash = 5e10 * u.kg
w = 2.04

labels = ['SO2', 'Ash', 'Sulfate', 'Sulfate_nodecay']
kk = [kso2, kash, ksulf]
MM = [Mso2, Mash]

# vert distribution function
f = lambda z: np.exp(-0.5 * (z-mu)**2/sigma**2) * 1/u.km

def compute_tracer_evolution(nlev = 'e3sm', zk = None, dt_inject = None, dt_decay = None, tmax = 365*5*u.day, analytic_only=False):
    '''
    Performs analytic and numerical solutions to SAI tracer tendency equations

    Parameters
    ----------
    nlev : int
        number of vertical levels from 10 to 26 km
    zk : Quantity array
        array of vertical positions in a length unit. If provided, nlev is ignored
    dt_inject : Quantity
        timestep during first 24 hours (injection period). Must have a unit
    dt_decay : Quantity
        timestep during decay (remaining 5 years). Must have a unit
    tmax : Quantity
        ending time of the decay period. Must have a unit. Defaults to 5 years
    analytic_only : bool
        whether or not to skip the numerical estimation and return only the analytic result. If
        True, the return for the numerical result is None
    '''
    # ---------- vertical distribution ----------
    z0 = 12 * u.km
    zT = 24 * u.km
    if(zk is None):
        if(nlev == 'e3sm'):
           zk = np.array([12066.939, 12551.622, 13032.191, 13509.571, 13984.683, 14458.534, 14931.965, 15405.396, 15879.088, 16353.361, 
                          16828.559, 17318.906, 17882.68, 18564.715, 19366.602, 20290.275, 21337.922, 22513.977, 23836.834, 25323.637, 
                          26980.406, 28813.516, 30832.33, 33047.85]) * u.m
        else:
            zk = (np.linspace(z0, zT, nlev)).to(u.m)
    else:
        zk = np.sort(zk.to(u.m))
    V = f(zk)
    sVk = np.sum(V)

    # ---------- normalization ----------
    A  = [MM[j] / (dt * sVk) for j in range(len(MM))]

    # ---------- tendencies ----------
    # analytic tendencies
    dmdt = [lambda mk, tt, dt: -kso2*mk + A[0] * V * int((tt+60*u.s) < tf),
            lambda mk, tt, dt: -kash*mk + A[1] * V * int((tt+60*u.s) < tf),
            lambda mk, tt, mso2: -ksulf * mk + w*kso2*mso2, 
            lambda mk, tt, mso2: kso2*mso2]

    # ---------- solve ----------

    # ----- numerical solutions
    # --- 1st order Euler
    euler_step = lambda dydt, mn, tn, dt: mn + dt * dydt(mn, tn, dt)
    euler_step_sulf = lambda dydt, mn, tn, dt, mso2n: mn + dt * dydt(mn, tn, mso2n)

    # ----- numerical injection 
    if(dt_inject is None):
        # default; do first day with fine resolution (~1.5 minute)
        #dt_inject = (24/1000)*u.hr
        dt_inject = (24/200)*u.hr

    start_time = t0.to(u.s)
    end_time = start_time + (1*u.day + dt_inject.to(u.day)).to(u.s)
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
        #dt_decay = (tmax/2000) # in days by default
        dt_decay = (tmax/200) # in days by default

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

    tmin = np.array([min(tt-t0, dt).m for tt in t]) * t[0].u
    m_analytic[0,:,:] = (A[0]/kk[0] * np.atleast_2d(V).T * np.exp(-kk[0]*(t-t0)) * (np.exp(kk[0] * tmin) - 1)).T
    m_analytic[1,:,:] = (A[1]/kk[1] * np.atleast_2d(V).T * np.exp(-kk[1]*(t-t0)) * (np.exp(kk[1] * tmin) - 1)).T
    m_analytic[2,:,:] = (w*A[0]/((kk[2] - kk[0])*kk[2]) * np.atleast_2d(V).T * \
                        np.exp(-kk[2]*(t-t0)) * (\
                        kk[0]*(1 - np.exp(kk[2] * tmin)) - kk[2]*np.exp((kk[2]-kk[0])*(t-t0)) * (1 - np.exp(kk[0] * tmin))\
                        )).T
    # VV only valid for w=1, ksulf=0
    m_analytic[3,:,:] = (A[0]/kk[0] * np.atleast_2d(V).T * np.exp(-kk[0]*(t-t0)) * (1 - np.exp(kk[0]*tmin) + np.exp(kk[0]*(t-t0))*kk[0]*tmin)).T
    
    return m_analytic, m_numerical, zk, V, t
 


# --------- select plots to render 
main_figure = False
z_figure = False
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
    m_analytic_100, _, zk_500, V_500, t = compute_tracer_evolution(500, analytic_only=True)
    m_analytic_8, _, zk_8, V_8, _ = compute_tracer_evolution(8, analytic_only=True)
    m_analytic_e3sm, _, zk_e3sm, V_e3sm, _ = compute_tracer_evolution('e3sm', analytic_only=True)
    A = 17*1e9*u.kg/(9*3600*u.s * np.sum(V_e3sm))
    dmdt = A * V_e3sm
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
    for ax in [ax1, ax2, ax3]:
        ax.set_facecolor((1,1,204/255))

    # ----- plot V(z) profile
    #  get higher res zk
    zk_fine = (np.linspace(10, 26, 1000) * u.km).to(u.m)
    V_fine = f(zk_fine)
    ax0.plot(V_fine.to(u.km**(-1)), zk_fine.to(u.km), '-k', lw=2.5)
    ax0.set_xlabel('$V(z)$ [km$^{-1}$]', fontsize=label_fs)
    ax0.set_ylim([10, 26])

    # ----- contour plots of denisty evolution
    #cf1 = ax1.contourf(t.to(u.day), zk.to(u.km), m_dens[1].T, cmap=cm, levels=clev, extend='both', vmin=-7.5)
    cf1 = ax1.contourf(t.to(u.day), zk.to(u.km), m_dens[1].T, cmap=cm, levels=clev, extend='both')
    ax1.contour(t.to(u.day), zk.to(u.km), m_dens[1].T, levels=clev_thin_contours, colors=['k']*len(clev_contours), linewidths=0.5)
    ax1.set_xlim([t0.to(u.day).m-1.33, t0.to(u.day).m+14])
    ax1.set_xlabel('$t$ [days]', fontsize=label_fs)
    claut.add_annotation_box(ax1, 'ash', loc='upper right', fs=label_fs)

    ax2.contourf(t.to(u.month), zk.to(u.km), m_dens[0].T, cmap=cm, levels=clev, extend='both', vmin=-7.5)
    ax2.contour(t.to(u.day), zk.to(u.km), m_dens[0].T, levels=clev_thin_contours, colors=['k']*len(clev_contours), linewidths=0.5)
    ax2.set_xlim([0, 12])
    ax2.set_xlabel('$t$ [months]', fontsize=label_fs)
    claut.add_annotation_box(ax2, 'SO$_2$', loc='upper right', fs=label_fs)

    ax3.contourf(t.to(u.year), zk.to(u.km), m_dens[2].T, cmap=cm, levels=clev, extend='both', vmin=-7.5)
    ax3.contour(t.to(u.day), zk.to(u.km), m_dens[2].T, levels=clev_thin_contours, colors=['k']*len(clev_contours), linewidths=0.5)
    ax3.set_xlim([-0.35, 5])
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
    ax4.legend(fancybox=False, fontsize=tick_fs, loc='upper right', ncol=2, frameon=False)
    ax4.set_xlim([0, 24])

    plt.tight_layout()

    # manually add contour labels after tight layout
    cc3 = ax3.contour(t.to(u.day), zk.to(u.km), m_dens[2].T, levels=clev_contours, colors=['k']*len(clev_contours))
    ax3.clabel(cc3, cc3.levels, inline=True, fmt='%.0f', fontsize=tick_fs)
    cc1 = ax1.contour(t.to(u.day), zk.to(u.km), m_dens[1].T, levels=clev_contours, colors=['k']*len(clev_contours))
    ax1.clabel(cc1, cc1.levels, inline=True, fmt='%.0f', fontsize=tick_fs, manual=((31, 15.2), (33.5, 18.3), (35.4, 18.9), (37.3, 19.4)))
    cc2 = ax2.contour(t.to(u.day), zk.to(u.km), m_dens[0].T, levels=clev_contours, colors=['k']*len(clev_contours))
    ax2.clabel(cc2, cc2.levels, inline=True, fmt='%.0f', fontsize=tick_fs, manual=((2.45, 18.53), (4, 18.9), (5.64, 19.3)))

    plt.savefig('figs/sai_column_main_fig_new.png', dpi=300)



# ==================================================================
# ====================== VARIOUS NLEV  =============================
if(z_figure):
    print('----- plotting z figure')
    
    # ---------- vis ----------
    fig = plt.figure(figsize=(8.2,5.5))
    gs = fig.add_gridspec(nrows=1, ncols=2)
    axV = fig.add_subplot(gs[0,0])
    axM = fig.add_subplot(gs[0,1:])

    # ---- get total masses
    m_tot_100 = np.sum(m_analytic_100.to('Mt'), axis=2)
    m_tot_8 = np.sum(m_analytic_8.to('Mt'), axis=2)
    m_tot_e3sm = np.sum(m_analytic_e3sm.to('Mt'), axis=2)

    # ---- global ax settings
    for ax in [axV, axM]:
        ax.grid(ls='--', alpha=grid_alph)
        ax.tick_params(axis='both', which='major', labelsize=tick_fs)
        claut.format_ticks(ax)

    # ----- plot V(z) profiles
    axV.plot(V_500.to(u.km**(-1)), zk_500.to(u.km), '-k', lw=3, label='analytic', alpha=1)
    axV.plot(V_e3sm[zk_e3sm.m<24000].to(u.km**(-1)), zk_e3sm[zk_e3sm.m<24000].to(u.km), '-or', lw=1.6, label='nlev={} (E3SMv2)'.format(len(zk_e3sm[zk_e3sm.m<24000])))
    axV.plot(V_8.to(u.km**(-1)), zk_8.to(u.km), '-sb', lw=1.6, label='nlev=8')
    axV.legend(fancybox=False, fontsize=tick_fs, loc='upper right')
    axV.set_xlabel('$V(z)$ [km$^{-1}$]', fontsize=label_fs)
    axV.set_ylabel('$z$ [km]', fontsize=label_fs)
    axV.set_ylim([12, 24])

    # ----- line plots of total mass
    axM.plot(t.to(u.month), (m_tot_8[0]-m_tot_100[0]) * 1e14, '-b', lw=2, 
                            label='M(nlev=8) $-$ M(analytic)')
    axM.plot(t.to(u.month), (m_tot_e3sm[0]-m_tot_100[0]) * 1e14, '-r', lw=2, 
                             label='M(nlev=24) $-$ M(analytic)')
    
    axM.set_xlabel('$t$ [months]', fontsize=label_fs)
    axM.set_ylabel('total SO$_2$ mass difference [$10^{-14}$ Mt]', fontsize=label_fs)
    axM.legend(fancybox=False, fontsize=tick_fs, loc='upper right')
    axM.set_xlim([0.9, 5])
    axM.yaxis.set_label_position("right")
    axM.yaxis.tick_right()
    #ticklabs = np.array(['{:.1f}'.format(t) for t in axM.get_yticks()])
    #ticklabs[axM.get_yticks() % 0.5 != 0] = ''
    #axM.set_yticklabels(ticklabs)

    plt.tight_layout()
    plt.savefig('figs/sai_column_nlev_fig_new.png', dpi=300)


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
    
    # whether or not to apply alternate parameter tunings that avoid aod saturation
    retune = True

    # ---------- vis ----------
    fig = plt.figure(figsize=(7,5))
    axes = fig.subplots(2,1,sharex=True,gridspec_kw={'height_ratios':[2,1],'hspace':0.1})
    ax1 = axes[0]
    ax2 = axes[1]

    # --- heating parameters
    sb= 5.670374419e-8 * u.W/u.m**2/u.K**4
    cp = const.Cp_d
    cell_area = (200*u.km)**2
    maxaod_z = 0.2 * u.km
    maxaod_idx = np.searchsorted(zk_heat, maxaod_z)
    zeta = 4.7e-4
    lat = 0 *np.pi/180
    Isw = 558.54*u.W/u.m**2 * np.cos(lat)
    Ilw = sb*(315*u.K - 60*u.K * np.sin(lat)**2)**4
    # four entries for so2, ash, sulfate, and sulfate_contol
    # turn off sulfate_control contribution via toggle at index 3 V
    
    if(retune):
        blw = np.array([1, 1, 1, 0]) * 0.0042 * u.m**2/u.kg # for single column
        bsw = np.array([1, 1, 1, 0]) * 0.3 * u.m**2/u.kg # for single column
    else:
        blw = np.array([1, 1, 1, 0]) * 0.062 * u.m**2/u.kg # for E3SM
        bsw = np.array([1, 1, 1, 0]) * 400 * u.m**2/u.kg # for E3SM
    dz = np.diff(np.hstack([[0*u.m],zk_heat]))
    
    # ---- get densities in kg/m^3
    cell_volume = (200*u.km)**2 * np.abs(dz)
    m_dens = np.zeros(m_analytic_heat.shape) * u.kg/(u.m**3)
    for k in range(len(zk_heat)-1):
        m_dens[:,:,k] = m_analytic_heat.to(u.kg)[:,:,k] / cell_volume.to(u.m**3)[k]
    # --- get air masses, assuming exponential atm
    #zk_mid = (zk_heat[:-1] + np.diff(zk)/2)
    dens_air = stdatm(zk_heat.to(u.m).m).density * u.kg/u.m**3
    m_air = (dens_air * cell_area * dz).to('kg')
    

    # ------- copute strat. heating
    # update heating normalization, ignoring ash
    dIlw = np.zeros(m_dens.shape[1:]) * u.W/u.m**2
    # for testing approximations...
    dIlw_simp = np.zeros(m_dens.shape[1:]) * u.W/u.m**2
    dIlw_aprx = np.zeros(m_dens.shape[1:]) * u.W/u.m**2
    dIlw_simp_aprx = np.zeros(m_dens.shape[1:]) * u.W/u.m**2
    
    for k in range(len(zk_heat)-1):
        
        test_impl = 1
        if(test_impl== 1):
            cell_od = np.zeros(m_dens.shape[1])    
            cell_od_sum = np.zeros(m_dens.shape[1])    
            for j in range(len(blw)):
                cell_od = cell_od + blw[j]*m_dens[j,:,k]*dz[k]
                for kk in range(k):
                    cell_od_sum = cell_od_sum + blw[j]*m_dens[j,:,kk]*dz[kk]
        print('CELL OD CEHCK 1: {}'.format(cell_od_sum[100]))
        
        test_impl = 2
        if(test_impl== 2):
            cell_od2 = np.zeros(m_dens.shape[1])    
            if(k == 0): cell_od_sum2 = np.zeros(m_dens.shape[1])    
            
            cell_od2 = blw[0] * m_dens[0,:,k] * dz[k] + \
                      blw[1] * m_dens[1,:,k] * dz[k] + \
                      blw[2] * m_dens[2,:,k] * dz[k]
            cell_od_sum2 = cell_od_sum2 + \
                      blw[0] * m_dens[0,:,k-1] * dz[k-1] + \
                      blw[1] * m_dens[1,:,k-1] * dz[k-1] + \
                      blw[2] * m_dens[2,:,k-1] * dz[k-1]
        print('CELL OD CEHCK 2: {}'.format(cell_od_sum2[100]))
        
        #if(0):
            # --sanity check for all equal blw VV
            #cell_od = cell_od + blw[0] * np.sum(m_dens[:,:,k], axis=0) * dz[k]
            #for kk in range(k):
            #    cell_od_sum = cell_od_sum + blw[0] * np.sum(m_dens[:,:,kk], axis=0) * dz[kk]

        # THIS INDEXING NEEDS TO REVERSE IN VERTICAL FOR E3SM
        # unapproximated (exp) form
        dIlw[:,k] = Ilw * np.exp(-cell_od_sum) * (1 - np.exp(-cell_od))
        # unapproximated simple form
        dIlw_simp[:,k] = Ilw  * (1 - np.exp(-cell_od))
        # approximated form
        dIlw_aprx[:,k] = Ilw * (1 - cell_od_sum) * (cell_od)
        # approximated simple form
        dIlw_simp_aprx[:,k] = Ilw * cell_od
    
    print('\nmax lw diff from unapprox, simple: {:.2e}'.format(np.nanmax(np.abs(dIlw_simp-dIlw))))
    print('max lw diff from unapprox, approx: {:.2e}'.format(np.nanmax(np.abs(dIlw_aprx-dIlw))))
    print('max lw diff from unapprox, simple approx: {:.2e}'.format(np.nanmax(np.abs(dIlw_simp_aprx-dIlw))))
    # sum across species and apply species flags
    strat_heating = np.zeros(np.shape(m_dens)[1:]) * u.K/u.day
    strat_heating[:,:] = strat_heating[:,:] + (cell_area*dIlw[:,:]/m_air[:] / cp)\
                                               .to(u.K/u.day)
    
    # for testing approximations...
    strat_heating_aprx = np.zeros(np.shape(m_dens)[1:]) * u.K/u.day
    strat_heating_simp = np.zeros(np.shape(m_dens)[1:]) * u.K/u.day
    strat_heating_simp_aprx = np.zeros(np.shape(m_dens)[1:]) * u.K/u.day
    
    # for testing approximations...
    strat_heating_aprx[:,:] = strat_heating_aprx[:,:] + (cell_area*dIlw_aprx[:,:]/m_air[:] / cp)\
                                               .to(u.K/u.day)
    strat_heating_simp[:,:] = strat_heating_simp[:,:] + (cell_area*dIlw_simp[:,:]/m_air[:] / cp)\
                                               .to(u.K/u.day)
    strat_heating_simp_aprx[:,:] = strat_heating_simp_aprx[:,:] + (cell_area*dIlw_simp_aprx[:,:]/m_air[:] / cp)\
                                               .to(u.K/u.day)
    print('max dT diff from unapprox, simple: {:.2e}'.format(np.nanmax(np.abs(strat_heating_simp-strat_heating))))
    print('max dT diff from unapprox, approx: {:.2e}'.format(np.nanmax(np.abs(strat_heating_aprx-strat_heating))))
    print('max dT diff from unapprox, simple approx: {:.2e}'.format(np.nanmax(np.abs(strat_heating_simp_aprx-strat_heating))))

    
    # ------- compute AOD
    nlev = len(zk_heat[zk_heat<maxaod_z])
    col_mass = np.flip(np.flip(m_analytic_heat.to('kg'), 2).cumsum(axis=2), axis=2)
    aod = np.zeros(np.shape(strat_heating)[0])
    dIsw = np.zeros(np.shape(strat_heating)[0]) * u.W/u.m**2
    aod_cooling = np.zeros(np.shape(strat_heating)[0:]) * u.K/u.day
    # scale flux density to cooling rate
    # THIS INDEXING NEEDS TO REVERSE IN VERTICAL FOR E3SM
    # sum across species and apply species flags
    for j in range(len(col_mass)):
        aod[:] = aod[:] + ((bsw[j] * col_mass[j,:,0])/cell_area).to_base_units()
    dIsw[:] = Isw*(np.exp(-aod[:]) - 1)
    for n in range(nlev):
        aod_cooling[:,n] = (zeta/cp * dIsw[:] * cell_area/np.sum(m_air[:nlev])).to(u.K/u.day)
   
    # --------------------------
    print('\nDONE COMPUTING')
    print('max strat heating: {}'.format(np.max(strat_heating)))
    print('max mix.ratio: {:.3E}'.format(np.max(m_dens / dens_air)))
    print('min surf cooling: {}'.format(np.min(aod_cooling)))
    print('max burden: {:.3E}'.format(np.max(col_mass)))
    print('max aod: {}\n'.format(np.max(aod)))

    # --- total heating rate
    dTdt = strat_heating + aod_cooling
    
    # ---- global ax settings
    ax1.grid(ls=':', alpha=0.3, color='k')
    ax2.grid(ls=':', alpha=0.3, color='k')
    ax1.tick_params(which='major', labelsize=tick_fs)
    ax2.tick_params(which='major', labelsize=tick_fs)
    claut.format_ticks(ax1)
    claut.format_ticks(ax2)
    cm = plt.cm.RdYlBu_r
    
    if(retune):
        clev1 = np.linspace(0.05, 0.35, 7)
        clev2 = np.linspace(-0.02, 0, 5)
        clev = np.hstack([clev2, clev1])
        clev_lab = [-0.02, -0.01, 0.1, 0.2, 0.3]
        clab_fmt = '.2f'
        clab_pos1 = ((4.95,14.75), (8.75, 16.05), (5.21,17.44))
        clab_pos2 = ((9.08, 203),)
    else:
        clev1 = np.linspace(0.5, 4, 8)
        clev2 = np.linspace(-0.1, 0, 5)
        clev = np.hstack([clev2, clev1])
        clev_lab = [-0.1, -0.05, 1, 2, 3]
        clab_fmt = '.2f'
        clab_pos1 = ((5.07,14.38), (4.66, 16.83))
        clab_pos2 = ((9.08, 203),)
    
    divnorm=colors.TwoSlopeNorm(vcenter=0.)
    bg = np.array([245,251,218])/255
    ax1.set_facecolor(bg)
    ax2.set_facecolor(bg)
    
    # ----- contour plots of heating rate
    cf1 = ax1.contourf(t_heat.to(u.month), zk_heat.to(u.km), dTdt.T, cmap=cm, levels=clev, extend='both', norm=divnorm)
    cfl1 = ax1.contour(t_heat.to(u.month), zk_heat.to(u.km), dTdt.T, levels=clev, colors = 'k', linewidths=0.75)
    cf2 = ax2.contourf(t_heat.to(u.month), zk_heat.to(u.m), dTdt.T, cmap=cm, levels=clev, extend='both', norm=divnorm)
    cfl2 = ax2.contour(t_heat.to(u.month), zk_heat.to(u.m), dTdt.T, levels=clev, colors = 'k', linewidths=0.75)
    
    # ----- color bar
    cbar = fig.colorbar(cf1, ax=axes.ravel().tolist(), orientation='vertical', location='right')
    cbar.set_label('dT/dt [K/day]', fontsize=label_fs)
    cbar.ax.tick_params(labelsize=tick_fs)

    
    # ----- contour labels
    ax1.clabel(cfl1, inline=True, fmt='%{}'.format(clab_fmt), fontsize=tick_fs, 
               inline_spacing=0, manual=iter(clab_pos1))
    ax2.clabel(cfl2, inline=True, fmt='%{}'.format(clab_fmt), fontsize=tick_fs, 
               inline_spacing=0, manual=iter(clab_pos2))
    
    ax2.set_xlabel('$t$ [months]', fontsize=label_fs)
    ax2.set_ylabel('$z$ [m]', fontsize=label_fs)
    ax1.set_ylabel('$z$ [km]', fontsize=label_fs)
    ax1.set_ylim([10, 23])
    ax1.set_xlim([0, 12])
    ax2.set_ylim([0, 400])
    ax2.set_xlim([0, 12])
    claut.add_annotation_box(ax1, 'stratosphere', loc='lower right', fs=tick_fs)
    claut.add_annotation_box(ax2, 'surface', loc='upper right', fs=tick_fs)

    fig.tight_layout()

    plt.savefig('figs/sai_column_heating_fig_new_{}.png'.format(['a','b'][int(not retune)]), 
                dpi=300, bbox_inches='tight')


plt.show()





