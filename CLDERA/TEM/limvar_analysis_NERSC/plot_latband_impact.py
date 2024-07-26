import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import sys
import pdb
from matplotlib import colors
import matplotlib as mpl
from matplotlib.ticker import ScalarFormatter
import matplotlib.patches as mpatches
import warnings
import variable_plotting_settings
from matplotlib.ticker import FuncFormatter
import math


var  = sys.argv[1]       # variable to plot
overlay_panel = bool(int(sys.argv[2])) # whether or not to render fourth panel with impact overlay on cf
try:
    q    = sys.argv[3]   # either 'aoa', 'e90', 'none', or not passed
except IndexError:
    q = None

loc = '/pscratch/sd/j/jhollo/E3SM/historical_data/limvar/analysis'

bands = ['SHpole', 'SHmid', 'tropics', 'NHmid', 'NHpole']

for bi in range(len(bands)):
    band = bands[bi]

    print('reading data')
    data   = xr.open_dataset('{}/data_ensmean_{}.nc'.format(loc, band))
    cf     = xr.open_dataset('{}/cf_ensmean_{}.nc'.format(loc, band))
    impact = xr.open_dataset('{}/impact_ensmean_{}.nc'.format(loc, band))
    pval   = xr.open_dataset('{}/pval_{}.nc'.format(loc, band))
    print('available data vars: {}'.format(list(data.data_vars)))

    tem_data      = xr.open_dataset('{}/tem_data_ensmean_{}.nc'.format(loc, band))
    tem_cf        = xr.open_dataset('{}/tem_cf_ensmean_{}.nc'.format(loc, band))
    tem_impact    = xr.open_dataset('{}/tem_impact_ensmean_{}.nc'.format(loc, band))
    tem_pval      = xr.open_dataset('{}/tem_pval_{}.nc'.format(loc, band))
    budget_data   = xr.open_dataset('{}/budget_data_ensmean_{}.nc'.format(loc, band))
    budget_cf     = xr.open_dataset('{}/budget_cf_ensmean_{}.nc'.format(loc, band))
    budget_impact = xr.open_dataset('{}/budget_impact_ensmean_{}.nc'.format(loc, band))
    budget_pval   = xr.open_dataset('{}/budget_pval_{}.nc'.format(loc, band))
    print('available TEM vars: {}'.format(list(tem_data.data_vars)))
    print('available TEM budget vars: {}'.format(list(budget_data.data_vars)))

    aoa_tem_data      = xr.open_dataset('{}/tem_data_ensmean_TRACER-AOA_{}.nc'.format(loc, band))
    aoa_tem_cf        = xr.open_dataset('{}/tem_cf_ensmean_TRACER-AOA_{}.nc'.format(loc, band))
    aoa_tem_impact    = xr.open_dataset('{}/tem_impact_ensmean_TRACER-AOA_{}.nc'.format(loc, band))
    aoa_tem_pval      = xr.open_dataset('{}/tem_pval_TRACER-AOA_{}.nc'.format(loc, band))
    aoa_budget_data   = xr.open_dataset('{}/budget_data_ensmean_TRACER-AOA_{}.nc'.format(loc, band))
    aoa_budget_cf     = xr.open_dataset('{}/budget_cf_ensmean_TRACER-AOA_{}.nc'.format(loc, band))
    aoa_budget_impact = xr.open_dataset('{}/budget_impact_ensmean_TRACER-AOA_{}.nc'.format(loc, band))
    aoa_budget_pval   = xr.open_dataset('{}/budget_pval_TRACER-AOA_{}.nc'.format(loc, band))
    print('available AOA TEM vars: {}'.format(list(aoa_tem_data.data_vars)))
    print('available AOA budget vars: {}'.format(list(aoa_budget_data.data_vars)))

    e90_tem_data      = xr.open_dataset('{}/tem_data_ensmean_TRACER-E90j_{}.nc'.format(loc, band))
    e90_tem_cf        = xr.open_dataset('{}/tem_cf_ensmean_TRACER-E90j_{}.nc'.format(loc, band))
    e90_tem_impact    = xr.open_dataset('{}/tem_impact_ensmean_TRACER-E90j_{}.nc'.format(loc, band))
    e90_tem_pval      = xr.open_dataset('{}/tem_pval_TRACER-E90j_{}.nc'.format(loc, band))
    e90_budget_data   = xr.open_dataset('{}/budget_data_ensmean_TRACER-E90j_{}.nc'.format(loc, band))
    e90_budget_cf     = xr.open_dataset('{}/budget_cf_ensmean_TRACER-E90j_{}.nc'.format(loc, band))
    e90_budget_impact = xr.open_dataset('{}/budget_impact_ensmean_TRACER-E90j_{}.nc'.format(loc, band))
    e90_budget_pval   = xr.open_dataset('{}/budget_pval_TRACER-E90j_{}.nc'.format(loc, band))
    print('available E90 TEM vars: {}'.format(list(e90_tem_data.data_vars)))
    print('available E90 budget vars: {}'.format(list(e90_budget_data.data_vars)))


    # find which dataset the variable belongs to, read
    print('extracting variable...')
    try:
        data, cf, impact, pval = [data[var], cf[var], impact[var], pval[var]]
    except KeyError: pass
    try:
        data, cf, impact, pval = [tem_data[var], tem_cf[var], tem_impact[var], tem_pval[var]]
    except KeyError: pass
    try:
        data, cf, impact, pval = [budget_data[var], budget_cf[var], budget_impact[var], budget_pval[var]]
    except KeyError: pass
    try:
        if(q == 'aoa'):
            data, cf, impact, pval = [aoa_tem_data[var], aoa_tem_cf[var], 
                                      aoa_tem_impact[var], aoa_tem_pval[var]]
            var = '{}_aoa'.format(var)
        elif(q == 'e90'):
            data, cf, impact, pval = [e90_tem_data[var], e90_tem_cf[var], 
                                      e90_tem_impact[var], e90_tem_pval[var]]
            var = '{}_e90'.format(var)
    except KeyError: pass
    try:
        if(q == 'aoa'):
            data, cf, impact, pval = [aoa_budget_data[var], aoa_budget_cf[var], 
                                      aoa_budget_impact[var], aoa_budget_pval[var]]
            var = '{}_aoa'.format(var)
        elif(q == 'e90'):
            data, cf, impact, pval = [e90_budget_data[var], e90_budget_cf[var], 
                                      e90_budget_impact[var], e90_budget_pval[var]]
            var = '{}_e90'.format(var)
    except KeyError: pass

    # --------------------------------------------------------------------------

    plev = data.plev
    time = data.time

    # --- do slicing
    ti, tf = 0, int(1.5*360) # plot for first year and a half
    tsl    = slice(ti, tf)
    data   = data.isel(time=tsl)
    cf     = cf.isel(time=tsl)
    impact = impact.isel(time=tsl)
    pval   = pval.isel(time=tsl)
    time   = data.time
    month = [tt.month+((tt.year-1991)*12)-6+tt.day/30 for tt in time.values]

    # --- do vertical slicing
    pmin, pmax = 1, 250
    levslice = slice(pmin, pmax)
    plev     = plev.sel(plev = levslice)
    data     = data.sel(plev = levslice)
    cf       = cf.sel(plev = levslice)
    impact   = impact.sel(plev = levslice)
    pval     = pval.sel(plev = levslice)

    # ---- get plotting settings for this variable
    opt = variable_plotting_settings.lat_p_plots
    data_levels    = opt[var]['data_lev']
    impact_levels  = opt[var]['impact_lev']
    data_norm      = opt[var]['data_norm']
    impact_norm    = opt[var]['impact_norm']
    scaling        = opt[var]['scaling']
    impact_scaling = opt[var]['impact_scaling']
    units          = opt[var]['units']
    impact_units   = opt[var]['impact_units']
    cmap           = opt[var]['cmap']
    fmt            = opt[var]['fmt']

    # ---- get default settings
    if(var == 'E90j'): var = 'E90'
    if(scaling is None):          scaling = 1
    if(impact_scaling is None):   impact_scaling = 1
    if(cmap is None):             cmap = 'Spectral_r'
    if(units is not None):        varstr = '{} [{}]'.format(var, units)
    else:                         varstr = var
    if(impact_units is not None): impactstr = '{} impact [{}]'.format(var, impact_units)
    else:                         impactstr = var

    titlestr = varstr
    if(len(varstr) > 10): 
        varstr = varstr.split(' [')[0] + '\n[' + varstr.split('[')[-1]
    if(len(impactstr) > 10): 
        impactstr = impactstr.split(' [')[0] + '\n[' + impactstr.split('[')[-1]

    # ---- scale variable
    if(scaling == 'log'):
        data   = np.log10(data)
        cf     = np.log10(cf)
    else:
        data   = data * scaling
        cf     = cf * scaling
    if(impact_scaling == 'log'): impact = np.log10(impact)
    else: impact = impact * impact_scaling

    # --- make colormap
    cmap_str = cmap
    if(data_levels is None): nlev = 12
    else: nlev = len(data_levels)
    cmap = plt.get_cmap(cmap)
    newcolors = np.vstack((
        cmap([0.0]),  # Extended color for lower bound
        cmap(np.linspace(1/nlev, 1-(1/nlev), cmap.N)),
        cmap([1.0])   # Extended color for upper bound
    ))
    cmap = colors.ListedColormap(newcolors)

    # --- make impact colormap
    if(impact_levels is None): nlev = 12
    else: nlev = len(impact_levels)
    impact_cmap = plt.get_cmap('RdBu_r')
    newcolors = np.vstack((
        impact_cmap([0.0]),  # Extended color for lower bound
        impact_cmap(np.linspace(1/nlev, 1-(1/nlev), impact_cmap.N)),
        impact_cmap([1.0])   # Extended color for upper bound
    ))
    impact_cmap = colors.ListedColormap(newcolors)

    # --- set parameters
    mpl.rcParams['hatch.linewidth'] = 0.25
    var_lw = 0.75
    pthresh  = 0.05 # p-value threshold for significance

    # --- font sizes
    suptitlefs = 15
    titlefs = 11
    labelfs = 9

    # --- function for formatting level labels on the colorbar
    def cbarfmt(v, pos):
        if(fmt is None):
            # Format the float with up to 2 decimal places
            vfmt = "{:.2f}".format(v)
            # Remove trailing zeros and the decimal point if not needed
            vfmt = vfmt.rstrip('0').rstrip('.')
            # if after this formatting, the number is zero, must have been smaller than 0.01
            if(float(vfmt) == 0 and v != 0): 
                vfmt = '{:.2e}'.format(v).replace('e-0','e-').replace('e+0','e+')
                # Remove trailing zeros and the decimal point if not needed
                val, exp = vfmt.split('e')
                val = val.rstrip('0').rstrip('.')
                vfmt = val + 'e' + exp
            #done
            return vfmt
        else:
            return(fmt(v))

    # --- function for automatically finding contour levels
    # if levels is None (default), automatically find levels
    # if levels provided, just construct the norm object and return
    def get_levels(x, pp=0.5, levels=None, norm=None):
        cmap_out = None
        if(levels is None):
            low, high = np.percentile(np.ravel(x), pp), np.percentile(np.ravel(x), 100-pp)
            if(low < 0 and high > 0):
                neglevels = np.linspace(low, 0, 6)
                poslevels = np.linspace(0, high, 6)[1:]
                levels = np.hstack([neglevels, poslevels])
                norm = colors.TwoSlopeNorm(vmin=min(levels), vcenter=0, vmax=max(levels))
            else:
                levels = np.linspace(np.percentile(np.ravel(x), pp), np.percentile(np.ravel(x), 100-pp), 12)
                norm = None
        elif(levels is not None and norm is None):
            low, high = min(levels), max(levels)
            if(low < 0 and high > 0):
                norm = colors.TwoSlopeNorm(vmin=min(levels), vcenter=0, vmax=max(levels))
            else:
                norm = None
        elif(levels is not None and norm is not None):
            if(norm == 'uneven'):
                norm = colors.BoundaryNorm(levels, cmap.N, extend='both')
            if(norm == 'log'):
                norm = colors.LogNorm(vmin=levels.min(), vmax=levels.max())
            if(norm == 'symlog'):
                norm = colors.SymLogNorm(linthresh=min(levels[levels>0]), linscale=1)
        return levels, norm

    # ---- configure levels
    data_levels, data_norm = get_levels(data, levels=data_levels, norm=data_norm)
    impact_levels, impact_norm = get_levels(impact, levels=impact_levels, norm=impact_norm)

    # --- make figure
    if(overlay_panel): num_plt = 4
    else:              num_plt = 3
    fig, ax = plt.subplots(num_plt, 1, figsize=(5, 20*num_plt/4), layout='constrained')

    # -------------------------------------------------------

    print('plotting band {}'.format(band))
    
    # ---- forced run
    x1 = data
    c1 = ax[0].contourf(month, plev, x1, cmap=cmap, norm=data_norm, levels=data_levels, 
                          extend='both')
    ax[0].contour(month, plev, x1, colors='k', levels=data_levels, alpha=0.5, linewidths=var_lw)
    # if zero-contour exists, make bold
    if(0 in data_levels):
        ax[0].contour(month, plev, x1, colors='k', levels=[0], alpha=0.5, linewidths=var_lw*1.5)


    # ---- counterfactual
    x2 = cf
    c2 = ax[1].contourf(month, plev, x2, cmap=cmap, norm=c1.norm, levels=c1.levels, extend='both')
    ax[1].contour(month, plev, x2, colors='k', levels=c1.levels, alpha=0.5, linewidths=var_lw)
    # if zero-contour exists, make bold
    if(0 in c1.levels):
        ax[1].contour(month, plev, x2, colors='k', levels=[0], alpha=0.5, linewidths=var_lw*2)
    
    # ---- impact
    x3 = impact
    c3 = ax[2].contourf(month, plev, x3, cmap=impact_cmap, norm=impact_norm, 
                          levels=impact_levels, extend='both')

    # ---- imapct significance
    x4 = pval
    c4 = ax[2].contour(month, plev, x4, colors='k', levels=[pthresh], linewidths=var_lw*1.25)
    c5 = ax[2].contourf(month, plev, x4, levels=[pthresh, x4.max()], 
                           hatches=['////'], extend='right', colors='none', alpha=0)

    if(overlay_panel):
        # ---- plot sign agreement of significant impact with counterfactual
        ax[3].contourf(month, plev, x2, cmap=cmap, norm=c2.norm, levels=c2.levels, extend='both')
     
        pmask = x4 < pthresh
        x5    = x3.where(pmask, other=0) * np.sign(x2)
        x5neg = x5.where(x5<0, other=0)
        x5pos = x5.where(x5>0, other=0)
        if(x5neg.values.min() < 0):
            c6 = ax[3].contourf(month, plev, x5neg, colors='none', levels=[x5neg.values.min(), -1e-20], 
                            hatches=['**'], extend='left')
        if(x5pos.values.max() > 0):
            c7 = ax[3].contourf(month, plev, x5pos, colors='none', levels=[1e-20, x5pos.values.max()], 
                            hatches=['OO'], extend='right')
        #c8 = ax[3].contour(month, plev, x5neg, colors='blue', levels=[c3.levels[c3.levels < 0][-1]], 
        #                linewidths=var_lw*1.25)
        #c9 = ax[3].contour(month, plev, x5pos, colors='red', levels=[c3.levels[c3.levels > 0][0]], 
        #                linewidths=var_lw*1.25)

    # --- legend
    with warnings.catch_warnings(action="ignore"):
        cc4 = plt.plot([0,0], [0,0], '-k', lw=var_lw*1.25)[0]
        cc5 = [mpatches.Patch(facecolor='w', hatch=pc.get_hatch()) for pc in c5.collections][0]
        if(overlay_panel):
            cc6 = [mpatches.Patch(facecolor='w', hatch=pc.get_hatch()) for pc in c6.collections][0]
            cc7 = [mpatches.Patch(facecolor='w', hatch=pc.get_hatch()) for pc in c7.collections][0]
    if(overlay_panel):
        legend_lines = [cc4, cc5, cc6, cc7]
        labels = ['pval = 0.05', 'pval > 0.05', 'impact decelerates', 'impact accelerates']
    else:
        legend_lines = [cc4, cc5]
        labels = ['pval = 0.05', 'pval > 0.05']
    fig.legend(legend_lines, labels, bbox_to_anchor=[0.55, 0], 
               loc='lower center', ncol=[2,4][month is None], fontsize=titlefs-1)

    # ---- colorbars
    cb1 = fig.colorbar(c2, ax=[ax[0]], location='left', aspect=7, pad=0.1, 
                       format=FuncFormatter(cbarfmt))
    cb1.set_label('10 Tg\n{}'.format(varstr), fontsize=titlefs)
    cb2 = fig.colorbar(c2, ax=[ax[1]], location='left', aspect=7, pad=0.1, 
                       format=FuncFormatter(cbarfmt)) 
    cb2.set_label('0 Tg\n{}'.format(varstr), fontsize=titlefs)
    cb3 = fig.colorbar(c3, ax=[ax[2]], location='left', aspect=7, pad=0.1, 
                       format=FuncFormatter(cbarfmt))
    cb3.set_label('(10 Tg - 0 Tg)\n{}'.format(impactstr), fontsize=titlefs)
    if(overlay_panel):
        cb4 = fig.colorbar(c2, ax=[ax[3]], location='left', aspect=7, pad=0.1, 
                           format=FuncFormatter(cbarfmt))
        cb4.set_label('0 Tg\n{}'.format(varstr), fontsize=titlefs)
    # if norm = uneven, then the contours an non-uniformly spaced, and every one of them
    # needs to be labeled. To make room for this, make font smaller
    if(isinstance(data_norm, colors.BoundaryNorm)):
        cb1.set_ticks(data_levels)
        cb1.set_ticklabels([cbarfmt(vv, None) for vv in cb1.get_ticks()])
        cb1.ax.tick_params(labelsize=cb1.ax.get_yticklabels()[0].get_size() * 0.75)
        cb2.set_ticks(data_levels)
        cb2.set_ticklabels([cbarfmt(vv, None) for vv in cb2.get_ticks()])
        cb2.ax.tick_params(labelsize=cb2.ax.get_yticklabels()[0].get_size() * 0.75)
        if(overlay_panel):
            cb4.set_ticks(data_levels)
            cb4.set_ticklabels([cbarfmt(vv, None) for vv in cb4.get_ticks()])
            cb4.ax.tick_params(labelsize=cb4.ax.get_yticklabels()[0].get_size() * 0.75)
    if(isinstance(impact_norm, colors.BoundaryNorm)):
        cb3.set_ticks(impact_levels)
        cb3.set_ticklabels([cbarfmt(vv, None) for vv in cb3.get_ticks()])
        cb3.ax.tick_params(labelsize=cb3.ax.get_yticklabels()[0].get_size() * 0.75)
    
    # --- format
    for axi in ax:
        axi.invert_yaxis()
        axi.set_yscale('log')
        axi.minorticks_off()
        yticks = np.array([0.1, 0.3, 1, 3, 10, 30, 100, 300, 700, 1000])
        yticks = yticks[np.logical_and(yticks >= plev.values.min(), yticks <= plev.values.max())]
        axi.set_yticks(yticks)
        axi.set_ylabel('pressure [hPa]')
        axi.yaxis.set_major_formatter(ScalarFormatter())
        axi.tick_params(labelright=False, labelleft=True, left=True, right=True)
        axi.set_xlim([0, np.max(month)])
    
    ax[0].tick_params(labeltop=True, labelbottom=False, bottom=True, top=True)
    ax[1].tick_params(labeltop=False, labelbottom=False, bottom=True, top=True)
    ax[2].tick_params(labeltop=False, labelbottom=False, bottom=True, top=True)
    ax[num_plt-1].tick_params(labeltop=False, labelbottom=True, bottom=True, top=True)
    ax[num_plt-1].set_xlabel('time [months]\n\n\n')

    ax[0].set_title('{} {}'.format(varstr, band), fontsize=titlefs)

overlay_str = ['', '_overlay'][overlay_panel is None]
plt.savefig('figs/{}_{}{}.png'.format(var, band, overlay_str), dpi=200)
plt.show()
