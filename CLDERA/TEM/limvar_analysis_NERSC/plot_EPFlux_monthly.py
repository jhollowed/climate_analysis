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
try:
    var2 = sys.argv[6]   # an optonal second variable to over-plot in contours
except IndexError:
    var2 = None
try:
    q2 = sys.argv[7]   # either 'aoa', 'e90', or 'none' for the optional second variable 
                       # to over-plot in contours
except IndexError:
    q2 = None
try:
    vec1 = sys.argv[8]   # optional variables to plot as vector field
    vec2 = sys.argv[9]
except IndexError:
    vec1, vec2 = None, None


# --------------------------------------------------------------------------


# ----- read the data
data, cf, impact, pval, coherence                          = putil.read_variable(var, q)
data2, cf2, impact2, pval2, coherence2                     = putil.read_variable(var2, q)
data_vec1, cf_vec1, impact_vec1, pval_vec1, coherence_vec1 = putil.read_variable(vec1, q)
data_vec2, cf_vec2, impact_vec2, pval_vec2, coherence_vec2 = putil.read_variable(vec2, q)

# ---- do pressure, time slicing
slice_args={'pmin':1,'pmax':400,'year':year,'month':month}
all_data = [data, cf, impact, pval, coherence]
data, cf, impact, pval, coherence = do_slicing(all_data, slice_args**)
all_data = [data2, cf2, impact2, pval2, coherence2]
data2, cf2, impact2, pval2, coherence2 = do_slicing(all_data, slice_args**)
all_data = [data_vec1, cf_vec1, impact_vec1, pval_vec1, coherence_vec1]
data_vec1, cf_vec1, impact_vec1, pval_vec1, coherence_vec1 = do_slicing(all_data, slice_args**)
all_data = [data_vec1, cf_vec1, impact_vec1, pval_vec1, coherence_vec1]
data_vec2, cf_vec2, impact_vec2, pval_vec2, coherence_vec2 = do_slicing(all_data, slice_args**)

# ---- get dims
time = data.time
lat = data.lat
plev = data.plev

def get_settings(var_in):
    # ---- get plotting settings for this variable
    opt = variable_plotting_settings.lat_p_plots
    return opt[var_in]['data_lev'], opt[var_in]['impact_lev'], opt[var_in]['data_norm'],\
           opt[var_in]['impact_norm'], opt[var_in]['scaling'], opt[var_in]['impact_scaling'], \
           opt[var_in]['units'], opt[var_in]['impact_units'], opt[var_in]['cmap'], opt[var_in]['fmt']

data_levels, impact_levels, data_norm, impact_norm, scaling, impact_scaling, units, impact_units, cmap, fmt = get_settings(var)
data2_levels, impact2_levels, data2_norm, impact2_norm, scaling2, impact2_scaling, units2, impact2_units, _, fmt2 = get_settings(var2)
data_vec1_levels, impact_vec1_levels, data_vec1_norm, impact_vec1_norm, scaling_vec1, impact_vec1_scaling, units_vec1, impact_vec1_units, _, fmt_vec1 = get_settings(vec1)
data_vec2_levels, impact_vec2_levels, data_vec2_norm, impact_vec2_norm, scaling_vec2, impact_vec2_scaling, units_vec2, impact_vec2_units, _, fmt_vec2 = get_settings(vec2)

def get_defaults(var_in, scaling_in, impact_scaling_in, cmap_in, units_in, impact_units_in):
    # ---- get default settings
    if(var_in == 'E90j'): var = 'E90'
    if(scaling_in is None):          scaling_in = 1
    if(impact_scaling_in is None):   impact_scaling_in = 1
    if(cmap_in is None):             cmap_in = 'Spectral_r'
    if(units_in is not None):        varstr_in = '{} [{}]'.format(var_in, units_in)
    else:                            varstr_in = var_in
    if(impact_units_in is not None): impactstr_in = '{} impact [{}]'.format(var_in, impact_units_in)
    else:                            impactstr_in = var_in
    if(len(varstr_in) > 10): 
        varstri_in = varstr_in.split(' [')[0] + '\n[' + varstr_in.split('[')[-1]
    if(len(impactstr_in) > 10): 
        impactstr_in = impactstr_in.split(' [')[0] + '\n[' + impactstr_in.split('[')[-1]
    return scaling_in, impact_scaling_in, cmap_in, varstr_in, impactstr_in

scaling, impact_scaling, cmap, varstr, impactstr = get_defaults(var, scaling, impact_scaling, cmap, units, impact_units) 
scaling2, impact2_scaling, _, varstr2, impactstr2 = get_defaults(var2, scaling2, impact2_scaling, None, units2, impact2_units) 
scaling_vec1, impact_vec1_scaling, _, varstr_vec1, impactstr_vec1 = get_defaults(vec1, scaling_vec1, impact_vec1_scaling, None, units_vec1, impact_vec1_units) 
scaling_vec2, impact_vec2_scaling, _, varstr_vec2, impactstr_vec2 = get_defaults(vec2, scaling_vec2, impact_vec2_scaling, None, units_vec2, impact_vec2_units) 

titlestr = varstr

# ---- plotting settings for tropopause
trop_lw = 3.2 # suppressed for now
trop_data_color='pink'
trop_cf_color='pink'
trop_ls = '-'
trop_label = 'tropopause'

def scale_vars(data_in, cf_in, impact_in, scaling_in, impact_scaling_in):
    # ---- scale variable
    if(scaling_in == 'log'):
        data_in   = np.log10(data_in)
        cf_in     = np.log10(cf_in)
    else:
        data_in   = data_in * scaling_in
        cf_in     = cf_in * scaling_in
    if(impact_scaling_in == 'log'): impact_in = np.log10(impact_in)
    else: impact_in = impact_in * impact_scaling_in
    return data_in, cf_in, impact_in

data, cf, impact = scale_vars(data, cf, impact, scaling, impact_scaling)
data2, cf2, impact2 = scale_vars(data2, cf2, impact2, scaling2, impact2_scaling)
data_vec1, cf_vec1, impact_vec1 = scale_vars(data_vec1, cf_vec1, impact_vec1, scaling_vec1, impact_vec1_scaling)
data_vec2, cf_vec2, impact_vec2 = scale_vars(data_vec2, cf_vec2, impact_vec2, scaling_vec2, impact_vec2_scaling)

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
var2_lw = 1.5
pthresh  = 0.05 # p-value threshold for significance
pval_levels = [0.025] # contours to plot in pvalue
coherence_levels = [0.799] # contours to plot in pvalue
coherence_color = 'yellow'
var2_color = 'grey'
vec_color = 'black'
sig_color = 'yellow'

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
data2_levels, data2_norm = get_levels(data2, levels=data2_levels, norm=data2_norm)
impact2_levels, impact2_norm = get_levels(impact2, levels=impact2_levels, norm=impact2_norm)
data_vec1_levels, data_vec1_norm = get_levels(data_vec1, levels=data_vec1_levels, norm=data_vec1_norm)
impact_vec1_levels, impact_vec1_norm = get_levels(impact_vec1, levels=impact_vec1_levels, norm=impact_vec1_norm)
data_vec2_levels, data_vec2_norm = get_levels(data_vec2, levels=data_vec2_levels, norm=data_vec2_norm)
impact_vec2_levels, impact_vec2_norm = get_levels(impact_vec2, levels=impact_vec2_levels, norm=impact_vec2_norm)

# ---- scale vector data
#sf = np.sqrt(data_vec1.plev) # scale factor
sf = 10
data_vecmag = (data_vec1**2 + data_vec2**2)**(1/2)
#data_vecmag = 1
data_vec1 = data_vec1/data_vecmag * sf
data_vec2 = data_vec2/data_vecmag * sf
cf_vecmag = (cf_vec1**2 + cf_vec2**2)**(1/2)
#cf_vecmag = 1
cf_vec1 = cf_vec1/cf_vecmag * sf
cf_vec2 = cf_vec2/cf_vecmag * sf
impact_vecmag = (impact_vec1**2 + impact_vec2**2)**(1/2)
impact_vec1 = impact_vec1/impact_vecmag * sf
impact_vec2 = impact_vec2/impact_vecmag * sf

# scale vertical component by 10x
sfv = 10 # vertical component scale factor
data_vec2 = data_vec2*sfv
cf_vec2 = cf_vec2*sfv
impact_vec2 = impact_vec2*sfv

# remove vectors where not significant
mask = np.logical_or(pval_vec1 < 0.025, pval_vec2 < 0.025)
impact_vec1 = impact_vec1.where(mask)
impact_vec2 = impact_vec2.where(mask)

# zero out var2 impact where not significant
mask = pval2 < 0.025
impact2 = impact2 * mask

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
    #ax[0,j].contour(lat, plev, x1.T, colors='k', levels=data_levels, alpha=0.5, linewidths=var_lw)
    # if zero-contour exists, make bold
    #if(0 in data_levels):
    #    ax[0,j].contour(lat, plev, x1.T, colors='k', levels=[0], alpha=0.5, linewidths=var_lw*1.5)
    # overlay tropopause
    ax[0, j].plot(lat, data_troppj, ls=trop_ls, color=trop_data_color, lw=trop_lw)
    # overlay var2
    x1 = data2.isel(time=j)
    ax[0,j].contour(lat, plev, x1.T, colors=var2_color, levels = data2_levels, linewidths=var2_lw)
    # overlay vector field
    x1 = [data_vec1.isel(time=j), data_vec2.isel(time=j)]
    ax[0,j].quiver(lat, plev, x1[0].T, x1[1].T, color=vec_color)


    # ---- counterfactual
    x2 = cf.isel(time=j)
    c2 = ax[1,j].contourf(lat, plev, x2.T, cmap=cmap, norm=c1.norm, levels=c1.levels, extend='both')
    #ax[1,j].contour(lat, plev, x2.T, colors='k', levels=c1.levels, alpha=0.5, linewidths=var_lw)
    # if zero-contour exists, make bold
    #if(0 in c1.levels):
    #    ax[1,j].contour(lat, plev, x2.T, colors='k', levels=[0], alpha=0.5, linewidths=var_lw*2)
    # overlay tropopause
    ax[1, j].plot(lat, cf_troppj, ls=trop_ls, color=trop_cf_color, lw=trop_lw) 
    # overlay var2
    x2 = cf2.isel(time=j)
    ax[1,j].contour(lat, plev, x2.T, colors=var2_color, levels = data2_levels, linewidths=var2_lw)
    # overlay vector field
    x2 = [cf_vec1.isel(time=j), cf_vec2.isel(time=j)]
    ax[1,j].quiver(lat, plev, x2[0].T, x2[1].T, color=vec_color)
    
    # ---- impact
    x3 = impact.isel(time=j)
    c3 = ax[2,j].contourf(lat, plev, x3.T, cmap=impact_cmap, norm=impact_norm, 
                          levels=impact_levels, extend='both')
    # overlay tropopause
    ax[2, j].plot(lat, cf_troppj, ls=trop_ls, color=trop_cf_color, lw=trop_lw)
    #ax[2, j].plot(lat, data_troppj, ls=trop_ls, color=trop_data_color, lw=trop_lw)
    # overlay var2
    x3 = impact2.isel(time=j)
    ax[2,j].contour(lat, plev, x3.T, colors=var2_color, levels = impact2_levels, linewidths=var2_lw)
    # overlay vector field
    x2 = [impact_vec1.isel(time=j), impact_vec2.isel(time=j)]
    ax[2,j].quiver(lat, plev, x2[0].T, x2[1].T, color=vec_color)

    # ---- imapct significance
    x4 = pval.isel(time=j)
    c4 = ax[2,j].contour(lat, plev, x4.T, colors=sig_color, levels=pval_levels, linewidths=var_lw*1.25)
    #c5 = ax[2,j].contourf(lat, plev, x4.T, levels=[pthresh, x4.max()], 
    #                       hatches=['////'], extend='right', colors='none', alpha=0) 
    # ---- imapct coherence
    x4_c = coherence.isel(time=j)
    #c4_c = ax[2,j].contour(lat, plev, x4_c.T, colors = coherence_color, 
    #                       levels=coherence_levels, linewidths=var_lw*1.25)
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
#plt.savefig('figs/{}_{}{}{}.png'.format(var, year, month_str, overlay_str), dpi=200)
plt.show()
