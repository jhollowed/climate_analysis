import os
import sys
import pdb
import glob
import warnings
import numpy as np
import xarray as xr
import artist_utils as aut
import cartopy.crs as ccrs
import climate_toolbox as ctb
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from cartopy.util import add_cyclic_point
from matplotlib.offsetbox import AnchoredText

cmap = plt.cm.rainbow


# ==========================================================================================


def vertical_slice(x, y, var_dict, ax, plot_zscale=True, inverty=True, logy=True, center_x=None,
                   xlabel=None, ylabel=None, title=None, xlim=None, ylim=None, gridlines=False, 
                   gridlinesArgs=None, cyclic=True, slice_at=None, slice_at_loc='lower left'):
    '''
    Plot the 2D vertical slice of a variable

    Parameters
    ----------
    x : 1D array, or string
        The horizontal coordinate. Internally, this will be used with the vertical
        cordainte y to construct grid via np.meshgrid. If array, this should be the 
        coordinate values. If string, this should be the name of a field present in
        var which gives the coordinate (for e.g. a xarray DataArray).
    y : 1D array, or string
        The vertical coordinate. Internally, this will be used with the horizontal
        cordainte y to construct grid via np.meshgrid. If array, this should be the 
        coordinate values. If string, this should be the name of a field present in
        var which gives the coordinate (for e.g. a xarray DataArray).
    var_dict : dict, or list of dicts
        A dictionary for each variable to be plotted, containing the following items:
        -- var : 2D array object
            The 2D data to plot (xarray DataArray, numpy array, etc.) 
        -- plotType : string, optional
            String specifying a plotting type, matching the name of an attribute of the
            matplotlib axis class (e.g. contour, contourf, ...). 
            Defaults to 'contourf'.
        -- plotArgs : dict, optional
            Args sent to the plotType plotting call, as a dict
            Deafults to matplotlib defaults in most cases. Some specific plotTypes inherit
            custom defaults for keyword args that aren't specified by this item; see function
            source for definitions
        -- colorFormatter : string, optional
            String specifying a function to format colors on the plot, matching the name of
            either an attribute of the axis class, or the parent figure class. 
            Defailts to one of the current associations, or else None, in which no color 
            formatting will be performed:
                plt.contour --> plt.clabel
                plt.contourf --> fig.colorbar
        -- colorArgs : dict, optional
            Args sent to colorFormatter, as a dict. 
            Defaults to matplotlib defaults in most cases; see source for other defaults set
            locally. If colorFormatter is None, this is ignored and a warning will be raised
            if included in the dict.
        Any optional items not included in the dict will be inserted and set to their stated
        default values. If var_dict is a list of dicts, each variable will be overlain on the axis.
    ax : pyplot axis object
        The axis on which to render the plot.
    plot_zscale : bool, optional
        Whether or not to include a second y axis of values converting the original y-axis 
        values from pressure to height, assuming an isothermal atmosphere. Assumed that the
        original y-axis in pressure. Defaults to true.
    inverty : bool, optional
        Whether or not to invert the y-axis (often used for a vertical pressure coordinate). 
        Defaults to true.
    logy : string, optional
        Whether or not to set the yscale to 'log'. Defaults to True.
    center_x: float
        x-coordiante on which to center the data. It will be assumed that the x-data is
        periodic, and defined in degrees. 
    xlim : list of length 2
        xlimits to impose on figure.
        Defaults to None, in which case the default chosen by the plotting call is used.
    ylim : list of length 2
        ylimits to impose on figure.
        Defaults to None, in which case the default chosen by the plotting call is used.
    xlabel : string
        The label fot the x-axis. Fontsize will default to 12.
        Defaults to 'x'. Set to an empty string '' to disable. 
    ylabel : string
        The label for the y-axis. Fontsize will default to 12.
        Default to 'y'. Set to an empty string '' to disable.
    title : string
        The title for the axis. Fontsize will default to 14. 
        Default is None, in which case the title is blank.
    gridlines : bool, optional
        Whether or not to turn on gridlines with a default customized configuration.
        Default is False
    gridlinesArgs : dict, optional 
        Args sent to cartopy gridlines, as a dict.
    cyclic : bool
        Whether or not to add a cyclic point to the data, to avoid gaps in any contour plots.
        Defaults to True.
    slice_at : string
        String to print in a text box in a plot corner, giving the position of the slice 
        (or stating a mean, etc). Defaults to 'SLICE AT ?' as a placeholder. Set to an empty
        string '' to disable the text box.
    '''
    
    # -------- prepare --------
    fig = ax.get_figure()

    # -------- define default plot args --------
    default_args = {
                    'contour'  :  {'levels':12, 'colors':'k', 'extend':'both'},
                    'contourf' :  {'levels':12, 'cmap':'rainbow','extend':'both'},
                    'clabel'   :  {'inline':True, 'fmt':'%.2f', 'fontsize':9},
                    'colorbar' :  {'ax':ax, 'location':'right', 'orientation':'vertical',
                                   'extend':'both', 'extendrect':False, 'format':'%.0f'},
                    'gridlines':  {'draw_labels':True, 'dms':True, 'x_inline':False, 'y_inline':False, 
                                   'color':'k', 'lw':0.3, 'alpha':0.5, 'xformatter':aut.LON_DEG_FORMATTER}
                   }
    color_formatters = {
                        'contour'  : 'clabel',
                        'contourf' : 'colorbar'
                       }
    if(xlabel is None): xlabel = 'x'
    if(ylabel is None): ylabel = 'y'
    if(slice_at is None): slice_at = 'SLICE AT ?'

    
    # -------- check inputs, add missing dict fields --------
    valid_keys = ['var', 'plotType', 'plotArgs', 'colorArgs', 'colorFormatter']
    if isinstance(var_dict, dict): var_dict = [var_dict]
    
    for i in range(len(var_dict)):
        d = var_dict[i]

        # ignore unrecognized keys
        for key in d.keys():
            if(key not in valid_keys):
                warnings.warn('var_dict key {} not recognized; ignoring'.format(key))

        # check types, create missing items
        assert 'var' in d.keys(), '2D field variable must be given in var_dict with key \'var\''
        if 'plotType' not in d.keys():
            d['plotType'] = 'contourf'
        else:
            assert isinstance(d['plotType'], str), \
                   'plotType must be given as a string, not {}'.format(type(d['plotType']))
        pl = d['plotType']
        
        if 'plotArgs' not in d.keys():
            d['plotArgs'] = {}
        else:
            assert isinstance(d['plotArgs'], dict), \
                   'plotArgs must be given as a dict, not {}'.format(type(d['plotArgs']))
        
        if 'colorArgs' not in d.keys():
            d['colorArgs'] = {}
        else:
            assert isinstance(d['colorArgs'], dict), \
                   'colorArgs must be given as a dict, not {}'.format(type(d['colorArgs']))
        
        if 'colorFormatter' not in d.keys():
            if(pl in color_formatters.keys()):
                d['colorFormatter'] = color_formatters[pl]
            else:
                d['colorFormatter'] = None
        else:
            assert isinstance(d['colorFormatter'], str) or d['colorFormatter'] is None, \
                   'colorArgs must be given as a dict, not {}'.format(type(d['colorArgs']))
        
        # merge specified kwargs with defaults
        if(pl in default_args.keys()):
            d['plotArgs'] = {**default_args[pl], **d['plotArgs']}
        if(d['colorFormatter'] in default_args.keys()):
            d['colorArgs'] = {**default_args[d['colorFormatter']], **d['colorArgs']}
        if gridlinesArgs is not None:
            gridlinesArgs = {**default_args['gridlines'], **gridlinesArgs}
        else:
            gridlinesArgs = default_args['gridlines']
  
    # -------- plot variables --------
    for i in range(len(var_dict)):
        if(cyclic):
            d['var'], xcyc = add_cyclic_point(d['var'], coord=x, axis=1)
    x = xcyc
    
    if(center_x is not None):
        # recenter on center_x in degrees, assuming periodicity in x
        xcen = x - center_x  # shifted coords only used to find rolling index
        shift_right, shift_left = False, False
        if(np.max(xcen) >= 180): shift_right = True
        if(np.min(xcen) <= -180): shift_left = True
        assert not(shift_left & shift_right), 'FAILED on centering at center_x;\
                                               data nonunique? x not in degrees? bug here?'
        if(shift_right): xroll = np.searchsorted(xcen, 180, side='right')
        if(shift_left): xroll = np.searchsorted(xcen, -180, side='left')
        x = np.roll(x, -xroll)                       # center x on x_center via a matrix "roll"
        x[x > (center_x + 180)] -= 360
        xlim = [center_x - 180, center_x + 180]

    X, Y = np.meshgrid(x, y)
    
    plots = np.empty(len(var_dict), dtype=object)
    for i in range(len(var_dict)):
        d = var_dict[i]
        if(center_x is not None):
            d['var'] = np.roll(d['var'], -xroll, axis=1) # roll the data for center_x
        plotter = getattr(ax, d['plotType'])
        plots[i] = plotter(X, Y, d['var'], **d['plotArgs'])
        # bold zero contour if exists
        if d['plotType'] == 'contour':
            try:
                if(not isinstance(plots[i].levels, list)):
                    zero = plots[i].levels.tolist().index(0)
                else:
                    zero = plots[i].levels.index(0)
                bold = plots[i].collections[zero].get_linewidth() * 1.5
                plots[i].collections[zero].set_linewidth(bold)
            except ValueError:
                pass
        
        
    # -------- format colors --------
    cf = np.empty(len(var_dict), dtype=object)
    for i in range(len(var_dict)):
        d = var_dict[i]
        if d['colorFormatter'] is not None:
            try:
                colorFormatter = getattr(ax, d['colorFormatter'])
            except AttributeError:
                try:
                    colorFormatter = getattr(fig, d['colorFormatter'])
                except AttributeError:
                    raise AttributeError('Neither object {} or {} has attribute {}'.format(
                                          type(ax), type(fig), d['colorFormatter']))
            cf[i] = colorFormatter(plots[i], **d['colorArgs'])
    
    # -------- format figure -------- 
    if(inverty): ax.invert_yaxis()
    if(logy): ax.set_yscale('log')
    if(xlim is not None): ax.set_xlim(xlim)
    if(ylim is not None): ax.set_ylim(ylim)
    if(xlabel != ''): ax.set_xlabel(xlabel, fontsize=12)
    if(ylabel != ''): ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(title, fontsize=14)
    
    # x, y tick labels formats assuming pressure vs. degrees
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y,pos: \
                                 ('{{:.{:1d}f}}'.format(int(np.maximum(-np.log10(y),0)))).format(y)))
    ax.xaxis.set_major_formatter(aut.LON_DEG_FORMATTER)
    
    if(gridlines):
        ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False, 
                     color='k', lw=0.3, alpha=0.75)
    if(plot_zscale):
        ylimz = ctb.ptoz(ax.get_ylim()).m/1000
        axz = ax.twinx()
        plotterz = getattr(axz, d['plotType'])
        plotterz(X, ctb.ptoz(Y).m, d['var'], **d['plotArgs'], alpha=0)  # pressure must be in hPa
        axz.set_ylim(ylimz)
        axz.set_ylabel(r'Z [km]')
    if(slice_at != ''): 
        text_box = AnchoredText(slice_at, frameon=True, loc=slice_at_loc, pad=0.5)
        text_box.set_zorder(100)
        plt.setp(text_box.patch, facecolor='white', alpha=1)
        ax.add_artist(text_box)
    return cf
    
    
# -------------------------------------------------------------


def horizontal_slice(x, y, var_dict, ax, projection=ccrs.Robinson(),  
                     xlabel=None, ylabel=None, title=None, xlim=None, ylim=None, 
                     gridlines=True, gridlinesArgs=None, coastlines=True, coastlinesArgs=None, 
                     cyclic=True, slice_at=None):
    '''
    Plot the 2D horizontal slice of a variable

    Parameters
    ----------
    x : 1D array, or string
        The horizontal coordinate. Internally, this will be used with the vertical
        cordainte y to construct grid via np.meshgrid. If array, this should be the 
        coordinate values. If string, this should be the name of a field present in
        var which gives the coordinate (for e.g. a xarray DataArray).
    y : 1D array, or string
        The vertical coordinate. Internally, this will be used with the horizontal
        cordainte y to construct grid via np.meshgrid. If array, this should be the 
        coordinate values. If string, this should be the name of a field present in
        var which gives the coordinate (for e.g. a xarray DataArray).
    var_dict : dict, or list of dicts
        A dictionary for each variable to be plotted, containing the following items:
        -- var : 2D array object
            The 2D data to plot (xarray DataArray, numpy array, etc.) 
        -- plotType : string, optional
            String specifying a plotting type, matching the name of an attribute of the
            matplotlib axis class (e.g. contour, contourf, ...). 
            Defaults to 'contourf'.
        -- plotArgs : dict, optional
            Args sent to the plotType plotting call, as a dict
            Defaults to matplotlib defaults in most cases. Some specific plotTypes inherit
            custom defaults for keyword args that aren't specified by this item; see function
            source for definitions
        -- colorFormatter : string, optional
            String specifying a function to format colors on the plot, matching the name of
            either an attribute of the axis class, or the parent figure class. 
            Defailts to one of the current associations, or else None, in which no color 
            formatting will be performed:
                plt.contour --> plt.clabel
                plt.contourf --> fig.colorbar
        -- colorArgs : dict, optional
            Args sent to colorFormatter, as a dict. 
            Defaults to matplotlib defaults in most cases; see source for other defaults set
            locally. If colorFormatter is None, this is ignored and a warning will be raised
            if included in the dict.
        Any optional items not included in the dict will be inserted and set to their stated
        default values. If var_dict is a list of dicts, each variable will be overlain on the axis.
    ax : pyplot axis object
        The axis on which to render the plot
    projection : cartopy ccrs object, optional
        Projection to apply to the slice, as a catropy.ccrs object instance.
        Default is ccrs.Robinson.
    ylim : list of length 2, optional
        ylimits to impose on figure.
        Defaults to None, in which case the default chosen by the plotting call is used.
    xlim : list of length 2, optional
        xlimits to impose on figure.
        Defaults to None, in which case the default chosen by the plotting call is used.
    xlabel : string
        The label fot the x-axis. Fontsize will default to 12.
        Defaults to 'x'. Set to an empty string '' to disable. 
    ylabel : string
        The label for the y-axis. Fontsize will default to 12.
        Default to 'y'. Set to an empty string '' to disable.
    title : string, optional
        The title for the axis. Fontsize will default to 14. 
        Default is None, in which case the title is blank.
    gridlines : bool, optional
        Whether or not to turn on gridlines with a default customized configuration.
        Default is True
    gridlinesArgs : dict, optional 
        Args sent to cartopy gridlines, as a dict.
    coastlines : bool, optional 
        Whether or not to plot coastlines. Defaults to True
    coastlinesArgs : dict, optional 
        Args sent to cartopy coastlines, as a dict.
    cyclic : bool
        Whether or not to add a cyclic point to the data, to avoid gaps in any contour plots.
        Defaults to True.
    slice_at : string
        String to print in a text box in a plot corner, giving the position of the slice 
        (or stating a mean, etc). Defaults to 'SLICE AT' as a placeholder. Set to an empty
        string '' to disable the text box.
    '''
    
    # -------- prepare --------
    fig = ax.get_figure()
    transform = ccrs.PlateCarree()

    
    # -------- define default plot args --------
    default_args = {
            'contour'   :  {'levels':12, 'colors':'k', 'extend':'both', 'transform':transform},
            'contourf'  :  {'levels':12, 'cmap':'rainbow', 'extend':'both', 'transform':transform},
            'clabel'    :  {'inline':True, 'fmt':'%.2f', 'fontsize':9, 'transform':transform},
                           'colorbar' :  {'ax':ax, 'orientation':'vertical', 'extend':'both', 
                           'extendrect':False, 'format':'%.0f', 'transform':transform},
            'gridlines' :  {'draw_labels':True, 'dms':True, 'x_inline':False, 'y_inline':False, 
                            'color':'k', 'lw':0.5, 'alpha':0.5, 'linestyle':':', 'crs':transform,
                            'xformatter':aut.LON_DEG_FORMATTER},
            'coastlines':  {'resolution':'110m', 'color':'k', 'linestyle':'-', 'alpha':0.75}
                   }
    color_formatters = {
                        'contour'  : 'clabel',
                        'contourf' : 'colorbar',
                       }
    if(xlabel is None): xlabel = 'x'
    if(ylabel is None): ylabel = 'y'
    if(slice_at is None): slice_at = 'SLICE AT ?'
    
    # -------- check inputs, add missing dict fields --------
    valid_keys = ['var', 'plotType', 'plotArgs', 'colorArgs', 'colorFormatter']
    if isinstance(var_dict, dict): var_dict = [var_dict]
    
    for i in range(len(var_dict)):
        d = var_dict[i]

        # ignore unrecognized keys
        for key in d.keys():
            if(key not in valid_keys):
                warnings.warn('var_dict key {} not recognized; ignoring'.format(key))

        # check types, create missing items
        assert 'var' in d.keys(), '2D field variable must be given in var_dict with key \'var\''
        if 'plotType' not in d.keys():
            d['plotType'] = 'contourf'
        else:
            assert isinstance(d['plotType'], str), \
                   'plotType must be given as a string, not {}'.format(type(d['plotType']))
        pl = d['plotType']
        
        if 'plotArgs' not in d.keys():
            d['plotArgs'] = {}
        else:
            assert isinstance(d['plotArgs'], dict), \
                   'plotArgs must be given as a dict, not {}'.format(type(d['plotArgs']))
        
        if 'colorArgs' not in d.keys():
            d['colorArgs'] = {}
        else:
            assert isinstance(d['colorArgs'], dict), \
                   'colorArgs must be given as a dict, not {}'.format(type(d['colorArgs']))
        
        if 'colorFormatter' not in d.keys():
            if(pl in color_formatters.keys()):
                d['colorFormatter'] = color_formatters[pl]
            else:
                d['colorFormatter'] = None
        else:
            assert isinstance(d['colorFormatter'], str) or d['colorFormatter'] is None, \
                   'colorArgs must be given as a dict, not {}'.format(type(d['colorArgs']))
        
        # merge specified kwargs with defaults
        if(pl in default_args.keys()):
            d['plotArgs'] = {**default_args[pl], **d['plotArgs']}
        if(d['colorFormatter'] in default_args.keys()):
            d['colorArgs'] = {**default_args[d['colorFormatter']], **d['colorArgs']}
        if gridlinesArgs is not None:
            gridlinesArgs = {**default_args['gridlines'], **gridlinesArgs}
        else:
            gridlinesArgs = default_args['gridlines']
        if coastlinesArgs is not None:
            coastlinesArgs = {**default_args['coastlines'], **coastlinesArgs}
        else:
            coastlinesArgs = default_args['coastlines']

    # -------- plot variables --------
     
    if(cyclic):
        d['var'], x = add_cyclic_point(d['var'], coord=x, axis=1)
    
    X, Y = np.meshgrid(x, y)
    
    plots = np.empty(len(var_dict), dtype=object)
    for i in range(len(var_dict)):
        d = var_dict[i]
        plotter = getattr(ax, d['plotType'])
        plots[i] = plotter(X, Y, d['var'], **d['plotArgs'])
        # bold zero contour if exists
        if d['plotType'] == 'contour':
            try: 
                if(not isinstance(plots[i].levels, list)):
                    zero = plots[i].levels.tolist().index(0)
                else:
                    zero = plots[i].levels.index(0)
                bold = plots[i].collections[zero].get_linewidth() * 1.5
                plots[i].collections[zero].set_linewidth(bold)
            except ValueError:
                pass
        
    # -------- format colors --------
    for i in range(len(var_dict)):
        d = var_dict[i]
        if d['colorFormatter'] is not None:
            try:
                colorFormatter = getattr(ax, d['colorFormatter'])
            except AttributeError:
                try:
                    colorFormatter = getattr(fig, d['colorFormatter'])
                except AttributeError:
                    raise AttributeError('Neither object {} or {} has attribute {}'.format(
                                          type(ax), type(fig), d['colorFormatter'])) 
            colorFormatter(plots[i], **d['colorArgs'])
    
    # -------- format figure --------
    if(xlim is not None): ax.set_xlim(xlim)
    if(ylim is not None): ax.set_ylim(ylim) 
    if(xlabel != ''): ax.set_xlabel(xlabel, fontsize=12)
    if(ylabel != ''): ax.set_ylabel(ylabel, fontsize=12)
    if(gridlines):
       gl = ax.gridlines(**gridlinesArgs)
       gl.xlabels_top = False
       gl.ylabels_right = False
    if(coastlines):
        ax.coastlines(**coastlinesArgs) 
    if(slice_at != ''): 
        text_box = AnchoredText(slice_at, frameon=True, loc='lower left', pad=0.5)
        text_box.set_zorder(100)
        plt.setp(text_box.patch, facecolor='white', alpha=1)
        ax.add_artist(text_box)
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(title, fontsize=14)
