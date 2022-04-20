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
import cartopy.feature as cfeature

cmap = plt.cm.rainbow


# ==========================================================================================


def vertical_slice(x, y, var_dict, ax=None, savefig=None, figsize=None, inverty=True, 
                   plot_zscale=True, yscale='log', xscale='linear', xlabel=None, ylabel=None, 
                   title=None, xlim=None, ylim=None, gridlines=False, gridlinesArgs=None, 
                   center_x=None):
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
    ax : pyplot axis object, optional
        The axis on which to render the plot.
        Defaults to None, in which case a new single-axis figure will be generated and
        returned.
    figsize : tuple
        tuple giving the figure size width, height in inches (same as matplotlib figsize arg to Figure).
        Default is None, in which case the figsize is defaulted by matplotlib.
    savefig : string, optional
        String giving the destination to write out a png of the plot.
        Defaults to None, in which case no plot if saved to file.
    inverty : bool, optional
        Whether or not to invert the y-axis (often used for a vertical pressure coordinate). 
        Defaults to true.
    plot_zscale : bool, optional
        Whether or not to include a second y axis of values converting the original y-axis 
        values from pressure to height, assuming an isothermal atmosphere. Assumed that the
        original y-axis in pressure. Defaults to true.
    yscale : string, optional
        String giving the yscale, matching a valid option for pyplot.yscale(). 
        Defautls to 'log',
    xscale : string, optional
        String giving the zscale, matching a valid option for pyplot.yscale(). 
        Defautls to 'linear'.
    ylim : list of length 2
        ylimits to impose on figure.
        Defaults to None, in which case the default chosen by the plotting call is used.
    xlim : list of length 2
        xlimits to impose on figure.
        Defaults to None, in which case the default chosen by the plotting call is used.
    xlabel : string
        The label fot the x-axis. Fontsize will default to 12.
        Default is None, in which case the label is blank.
    ylabel : string
        The label for the y-axis. Fontsize will default to 12.
        Default is None, in which case the label is blank.
    title : string
        The title for the axis. Fontsize will default to 14. 
        Default is None, in which case the title is blank.
    gridlines : bool, optional
        Whether or not to turn on gridlines with a default customized configuration.
        Default is False
    gridlinesArgs : dict, optional 
        Args sent to cartopy gridlines, as a dict.
    center_x: float
        x-coordiante on which to center the data. It will be assumed that the x-data is
        periodic, and defined in degrees. 
        
    Returns
    -------
    matplotlib Figure object
        Returns a figure if ax=None.
    '''
    
    # -------- prepare --------
    if ax is None:
        if figsize is None: fig = plt.figure()
        else: fig = plt.figure(figsize)
        ax = fig.add_subplot(111)
        return_fig = True
    else:
        fig = ax.get_figure()
        return_fig = False


    # -------- define default plot args --------
    default_args = {
                    'contour'  :  {'levels':12, 'colors':'k', 'extend':'both'},
                    'contourf' :  {'levels':12, 'cmap':'rainbow','extend':'both'},
                    'clabel'   :  {'inline':True, 'fmt':'%.2f', 'fontsize':9},
                    'colorbar' :  {'ax':ax, 'location':'right', 'orientation':'vertical',
                                   'extend':'both', 'extendrect':True, 'format':'%.0f'},
                    'gridlines':  {'draw_labels':True, 'dms':True, 'x_inline':False, 'y_inline':False, 
                                   'color':'k', 'lw':0.3, 'alpha':0.5}
                   }
    color_formatters = {
                        'contour'  : 'clabel',
                        'contourf' : 'colorbar'
                       }
    
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
    
    if(center_x is not None):
        # recenter on center_x in degrees assuming periodicity in x
        #xmax = center_x + 180
        #xroll = np.searchsorted(x, xmax)
        #xorig = x
        #x = x - center_x
        #xmin = np.min(x)

        #shift_x_origin   = lambda x: np.hstack([x[xmax:], x[xmax:] - 360])
        #unshift_x_origin = lambda x: np.hstack([ (x+center_x)[x>180], [x>180] + 360])
        
        #x = shift_x_origin(x)
        #x = np.roll(x, -xroll)
        pass
    else:
        #xlabs = None
        pass
    
    X, Y = np.meshgrid(x, y)
    
    plots = np.empty(len(var_dict), dtype=object)
    for i in range(len(var_dict)):
        d = var_dict[i]
        plotter = getattr(ax, d['plotType'])
        plots[i] = plotter(X, Y, d['var'], **d['plotArgs'])
        
        
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
    if(inverty): ax.invert_yaxis()
    if(xlim is not None): ax.set_xlim(xlim)
    if(ylim is not None): ax.set_ylim(ylim)
    ax.set_yscale(yscale)
    ax.set_xscale(xscale)
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(title, fontsize=14)
    if(gridlines):
        ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False, 
                     color='k', lw=0.3, alpha=0.75)
    if(plot_zscale):
        ylimz = ctb.ptoz(ax.get_ylim()).m/1000
        axz = ax.twinx()
        plotterz = getattr(axz, d['plotType'])
        # pressure must be in hPa
        plotterz(X, ctb.ptoz(Y).m, d['var'], **d['plotArgs'], alpha=0)
        axz.set_ylim(ylimz)
        axz.set_ylabel(r'Z [km]')
    if(center_x is not None):
        pdb.set_trace()
        ax.set_xticks(xlabs)
        ax.set_xticklabels(xlabs)

    if(savefig is not None):
        plt.savefig(savefig, dpi=300)
    if(return_fig):
        return fig
    
    
# -------------------------------------------------------------


def horizontal_slice(x, y, var_dict, ax=None, savefig=None, projection=ccrs.Robinson(), 
                     transform=ccrs.PlateCarree(), 
                     figsize=None, xlabel=None, ylabel=None, title=None, xlim=None, ylim=None, 
                     gridlines=True, gridlinesArgs=None, coastlines=True):
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
    ax : pyplot axis object, optional
        The axis on which to render the plot
        Defaults to None, in which case a new single-axis figure will be generated and
        returned.
    projection : cartopy ccrs object, optional
        Projection to apply to the slice, as a catropy.ccrs object instance.
        Default is ccrs.Robinson.
    transform : cartopy ccrs object, optional
        Transform to apply to plots on the slice, as a catropy.ccrs object instance.
        Default is ccrs.PlateCarree.
    figsize : tuple
        tuple giving the figure size width, height in inches
        Default is None, in which case the figsize is defaulted by matplotlib.
    savefig : string, optional
        String giving the destination to write out a png of the plot.
        Defaults to None, in which case no plot if saved to file.
    ylim : list of length 2, optional
        ylimits to impose on figure.
        Defaults to None, in which case the default chosen by the plotting call is used.
    xlim : list of length 2, optional
        xlimits to impose on figure.
        Defaults to None, in which case the default chosen by the plotting call is used.
    xlabel : string, optional
        The label fot the x-axis. Fontsize will default to 12.
        Default is None, in which case the label is blank.
    ylabel : string, optional
        The label for the y-axis. Fontsize will default to 12.
        Default is None, in which case the label is blank.
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

    Returns
    -------
    matplotlib Figure object
        Returns a figure if ax=None.
    '''
    
    # -------- prepare --------
    if ax is None:
        if figsize is None: fig = plt.figure()
        else: fig = plt.figure(figsize)
        ax = fig.add_subplot(111, projection=projection)
        return_fig = True
    else:
        fig = ax.get_figure()
        return_fig = False

    X, Y = np.meshgrid(x, y)
    
    # -------- define default plot args --------
    default_args = {
                    'contour'  :  {'levels':12, 'colors':'k', 'extend':'both'},
                    'contourf' :  {'levels':12, 'cmap':'rainbow', 'extend':'both'},
                    'clabel'   :  {'inline':True, 'fmt':'%.2f', 'fontsize':9},
                    'colorbar' :  {'ax':ax, 'orientation':'vertical', 'extend':'both', 
                                   'extendrect':True, 'format':'%.0f'},
                    'gridlines':  {'draw_labels':True, 'dms':True, 'x_inline':False, 'y_inline':False, 
                                   'color':'k', 'lw':0.5, 'alpha':0.5, 'linestyle':':', 'crs':transform}
                   }
    color_formatters = {
                        'contour'  : 'clabel',
                        'contourf' : 'colorbar',
                       }
    
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
    plots = np.empty(len(var_dict), dtype=object)
    for i in range(len(var_dict)):
        d = var_dict[i]
        plotter = getattr(ax, d['plotType'])
        plots[i] = plotter(X, Y, d['var'], transform=transform, **d['plotArgs'])
        
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
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(title, fontsize=14)
    if(gridlines):
       gl = ax.gridlines(**gridlinesArgs)
       gl.xlabels_top = False
       gl.ylabels_right = False
    if(coastlines):
        ax.coastlines(resolution='110m', color='k', linestyle='-', alpha=0.5)

    if(savefig is not None):
        plt.savefig(savefig, dpi=300)
    if(return_fig):
        return fig
