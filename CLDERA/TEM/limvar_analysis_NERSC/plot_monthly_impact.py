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
import plotting_utils as putil


var  = sys.argv[1]       # variable to plot
year = int(sys.argv[2])  # year to plot
overlay_panel = bool(int(sys.argv[3])) # whether or not to render fourth panel with impact overlay on cf
try:                     # if month not provided, plot entire calendar
    month = int(sys.argv[4])
except (IndexError, ValueError):
    month = None
try:
    q    = sys.argv[5]   # either 'aoa', 'e90', 'none', or not passed
except IndexError:
    q = None


# ----- read the data
data, cf, impact, pval, coherence = putil.read_variable(var)

# --------------------------------------------------------------------------

lat = data.lat
plev = data.plev
time = data.time
years = np.array([t.year for t in time.values])

# --- do time slicing
if(year == 1991): ti, tf = 0, 6
elif(year == 1992): ti, tf = 7, 7+12
elif(year == 1993): ti, tf, = 19, 19+12
elif(year == 1994): ti, tf, = 30, len(years)-1
tsl       = slice(ti, tf)
data      = data.isel(time=tsl)
cf        = cf.isel(time=tsl)
impact    = impact.isel(time=tsl)
pval      = pval.isel(time=tsl)
coherence = coherence.isel(time=tsl)
time      = data.time

# extract single month if provided
if(month is not None):
    tidx      = [t.month for t in time.values].index(month)
    tidx      = slice(tidx, tidx+1)
    data      = data.isel(time=tidx)
    cf        = cf.isel(time=tidx)
    impact    = impact.isel(time=tidx)
    pval      = pval.isel(time=tidx)
    coherence = coherence.isel(time=tidx)
    time      = data.time

# --- do vertical slicing
pmin, pmax = 1, 400
levslice  = slice(pmin, pmax)
plev      = plev.sel(plev = levslice)
data      = data.sel(plev = levslice)
cf        = cf.sel(plev = levslice)
impact    = impact.sel(plev = levslice)
pval      = pval.sel(plev = levslice)
coherence = coherence.sel(plev = levslice)

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

# ---- plotting settings for tropopause
trop_lw = 3.2 # suppressed for now
trop_data_color='pink'
trop_cf_color='pink'
trop_ls = '-'
trop_label = 'tropopause'

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
pval_levels = [0.025, 0.05] # contours to plot in pvalue
coherence_levels = [0.799] # contours to plot in pvalue
coherence_color = 'yellow'

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
fig, ax = plt.subplots(num_plt, len(time), figsize=(4*len(time), 7*num_plt/4), layout='constrained')
fig.suptitle(titlestr)
if(len(time)==1):
    # compatability of loops below when running for single month
    ax = np.atleast_2d(ax).T

# -------------------------------------------------------

for j in range(len(time)):

    print('plotting month {}...'.format(j+1))

    # get tropopause data
    data_troppj = data_tropp.isel(time=j)
    cf_troppj = cf_tropp.isel(time=j)
    
    # ---- forced run
    x1 = data.isel(time=j)
    c1 = ax[0,j].contourf(lat, plev, x1.T, cmap=cmap, norm=data_norm, levels=data_levels, 
                          extend='both')
    ax[0,j].contour(lat, plev, x1.T, colors='k', levels=data_levels, alpha=0.5, linewidths=var_lw)
    # if zero-contour exists, make bold
    if(0 in data_levels):
        ax[0,j].contour(lat, plev, x1.T, colors='k', levels=[0], alpha=0.5, linewidths=var_lw*1.5)
    # overlay tropopause
    ax[0, j].plot(lat, data_troppj, ls=trop_ls, color=trop_data_color, lw=trop_lw)


    # ---- counterfactual
    x2 = cf.isel(time=j)
    c2 = ax[1,j].contourf(lat, plev, x2.T, cmap=cmap, norm=c1.norm, levels=c1.levels, extend='both')
    ax[1,j].contour(lat, plev, x2.T, colors='k', levels=c1.levels, alpha=0.5, linewidths=var_lw)
    # if zero-contour exists, make bold
    if(0 in c1.levels):
        ax[1,j].contour(lat, plev, x2.T, colors='k', levels=[0], alpha=0.5, linewidths=var_lw*2)
    # overlay tropopause
    ax[1, j].plot(lat, cf_troppj, ls=trop_ls, color=trop_cf_color, lw=trop_lw)
    
    # ---- impact
    x3 = impact.isel(time=j)
    c3 = ax[2,j].contourf(lat, plev, x3.T, cmap=impact_cmap, norm=impact_norm, 
                          levels=impact_levels, extend='both')
    # overlay tropopause
    ax[2, j].plot(lat, cf_troppj, ls=trop_ls, color=trop_cf_color, lw=trop_lw)
    #ax[2, j].plot(lat, data_troppj, ls=trop_ls, color=trop_data_color, lw=trop_lw)

    # ---- imapct significance
    x4 = pval.isel(time=j)
    c4 = ax[2,j].contour(lat, plev, x4.T, colors='k', levels=pval_levels, linewidths=var_lw*1.25)
    #c5 = ax[2,j].contourf(lat, plev, x4.T, levels=[pthresh, x4.max()], 
    #                       hatches=['////'], extend='right', colors='none', alpha=0) 
    # ---- imapct coherence
    x4_c = coherence.isel(time=j)
    c4_c = ax[2,j].contour(lat, plev, x4_c.T, colors = coherence_color, 
                           levels=coherence_levels, linewidths=var_lw*1.25)
    c5 = ax[2,j].contourf(lat, plev, x4_c.T, levels=[0, coherence_levels[0]], 
                           hatches=['////'], extend='right', colors='none', alpha=0) 

    if(overlay_panel):
        # ---- plot sign agVreement of significant impact with counterfactual
        ax[3,j].contourf(lat, plev, x2.T, cmap=cmap, norm=c2.norm, levels=c2.levels, extend='both')
        
        # overlay tropopause
        ax[3, j].plot(lat, cf_troppj, ls=trop_ls, color=trop_cf_color, lw=trop_lw)
        #ax[3, j].plot(lat, data_troppj, ls=trop_ls, color=trop_data_color, lw=trop_lw)
     
        pmask = x4 < pthresh
        x5    = x3.where(pmask, other=0) * np.sign(x2)
        x5neg = x5.where(x5<0, other=0)
        x5pos = x5.where(x5>0, other=0)
        if(x5neg.values.min() < 0):
            c6 = ax[3,j].contourf(lat, plev, x5neg.T, colors='none', 
                                  levels=[x5neg.values.min(), -1e-20], hatches=['**'], extend='left')
        else:  
            c6 = None
        if(x5pos.values.max() > 0):
            c7 = ax[3,j].contourf(lat, plev, x5pos.T, colors='none', 
                                  levels=[1e-20, x5pos.values.max()], hatches=['OO'], extend='right')
        else:
            c7 = None
        #c8 = ax[3,j].contour(lat, plev, x5neg.T, colors='blue', levels=[c3.levels[c3.levels < 0][-1]], 
        #                linewidths=var_lw*1.25)
        #c9 = ax[3,j].contour(lat, plev, x5pos.T, colors='red', levels=[c3.levels[c3.levels > 0][0]], 
        #                linewidths=var_lw*1.25)

    # --- legend
    if(j == 0):
        with warnings.catch_warnings(action="ignore"):
            cc4   = plt.plot([0,0], [0,0], '-k', lw=var_lw*1.25)[0]
            cc4_2 = plt.plot([0,0], [0,0], '-', color=coherence_color, lw=var_lw*1.25)[0]
            cc5   = [mpatches.Patch(facecolor='w', hatch=pc.get_hatch()) for pc in c5.collections][0]
            if(overlay_panel):
                if(c6 is not None):
                    cc6 = [mpatches.Patch(facecolor='w', hatch=pc.get_hatch()) \
                                                   for pc in c6.collections][0]
                if(c7 is not None):
                    cc7 = [mpatches.Patch(facecolor='w', hatch=pc.get_hatch()) \
                                                   for pc in c7.collections][0]
        #legend_lines = [cc4, cc5, cc4_2]
        #labels = ['pval = 0.05', 'pval > 0.05', 'coherence = 0.8']
        legend_lines = [cc4, cc5]
        labels = ['pval = 0.05', 'coherence < 0.8']
        if(overlay_panel):
            if(c6 is not None):
                legend_lines.extend([cc6])
                labels.extend('impact decelerates')
            if(c7 is not None):
                legend_lines.extend([cc7])
                labels.extend('impact accelerates')
        fig.legend(legend_lines, labels, bbox_to_anchor=[0.55, 0], 
                   loc='lower center', ncol=[2,4][month is None], fontsize=titlefs-1)

    # ---- colorbars
    if(j == 0):
        cb1 = fig.colorbar(c2, ax=[ax[0,j]], location='left', aspect=7, pad=0.1, 
                           format=FuncFormatter(cbarfmt))
        cb1.set_label('10 Tg\n{}'.format(varstr), fontsize=titlefs)
        cb2 = fig.colorbar(c2, ax=[ax[1,j]], location='left', aspect=7, pad=0.1, 
                           format=FuncFormatter(cbarfmt)) 
        cb2.set_label('0 Tg\n{}'.format(varstr), fontsize=titlefs)
        cb3 = fig.colorbar(c3, ax=[ax[2,j]], location='left', aspect=7, pad=0.1, 
                           format=FuncFormatter(cbarfmt))
        cb3.set_label('(10 Tg - 0 Tg)\n{}'.format(impactstr), fontsize=titlefs)
        if(overlay_panel):
            cb4 = fig.colorbar(c2, ax=[ax[3,j]], location='left', aspect=7, pad=0.1, 
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
    for axi in ax[:,j]:
        axi.invert_yaxis()
        axi.set_yscale('log')
        axi.minorticks_off()
        axi.set_xticks([-90, -60, -30, 0, 30, 60, 90])
        yticks = np.array([0.1, 0.3, 1, 3, 10, 30, 100, 300, 700, 1000])
        yticks = yticks[np.logical_and(yticks >= plev.values.min(), yticks <= plev.values.max())]
        axi.set_yticks(yticks)
        if(j == 0):
            axi.set_ylabel('pressure [hPa]')
            axi.yaxis.set_major_formatter(ScalarFormatter())
            axi.tick_params(labelright=False, labelleft=True, left=True, right=True)
        elif(j == len(time)-1):
            axi.set_ylabel('pressure[hPa]')
            axi.yaxis.set_major_formatter(ScalarFormatter())
            axi.yaxis.set_label_position("right")
            axi.tick_params(labelright=True, labelleft=False, left=True, right=True)
        else:
            axi.tick_params(labelright=False, labelleft=False, left=True, right=True)
    
    ax[0,j].tick_params(labeltop=True, labelbottom=False, bottom=True, top=True)
    ax[1,j].tick_params(labeltop=False, labelbottom=False, bottom=True, top=True)
    ax[2,j].tick_params(labeltop=False, labelbottom=False, bottom=True, top=True)
    ax[num_plt-1,j].tick_params(labeltop=False, labelbottom=True, bottom=True, top=True)
    ax[num_plt-1,j].set_xlabel('lat\n\n\n')

    ax[0,j].set_title('{} {}'.format(time.values[j].strftime('%b'), time.values[j].year), 
                              fontsize=titlefs)

# ----- save figure to file
print('saving figure...')
month_str = ['_{}'.format(month), ''][month is None]
overlay_str = ['_overlay', ''][overlay_panel is None]
plt.savefig('figs/coherence/{}_{}{}{}.png'.format(var, year, month_str, overlay_str), dpi=200)
#plt.show()
