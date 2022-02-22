import os
import sys
import pdb
import glob
import warnings
import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature


cmap = plt.cm.rainbow

# ==========================================================================================


def vertical_slice(var_dict, x, y, ax=None, savefig=None, inverty=True, yscale='log', xscale='linear',
                   xlabel=None, ylabel=None, title=None, xlim=None, ylim=None):
    '''
    Plot the 2D vertical slice of a variable

    Parameters
    ----------
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
    ax : pyplot axis object, optional
        The axis on which to render the plot.
        Defaults to None, in which case a new single-axis figure will be generated and
        returned.
    savefig : string, optional
        String giving the destination to write out a png of the plot.
        Defaults to None, in which case no plot if saved to file.
    inverty : bool, optional
        Whether or not to invert the y-axis (often used for a vertical pressure coordinate). 
        Defaults to true.
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
        
    Returns
    -------
    matplotlib Figure object
        Returns a figure if ax=None.
    '''
    
    # -------- prepare --------
    if ax is None:
        fig = plt.Figure()
        ax = fig.add_sbuplot(111)
        return_fig = True
    else:
        fig = ax.get_figure()
        return_fig = False

    X, Y = np.meshgrid(x, y)
    
    # -------- define default plot args --------
    default_args = {
                    'contour'  :  {'levels':12, 'cmap':plt.cm.rainbow, 'extend':'both'}
                    'contourf' :  {'levels':12, 'colors':'k'}
                    'clabel'   :  {'inline':True, 'fmt':'%.2f', 'fontsize':9}
                    'colorbar' :  {'cax':ax, 'location':'right', 'orientation':'vertical',
                                   'extend':'both', 'extendrec':True, 'format':'%.0f'}
                   }
    color_formatters = {
                        'contour'  : 'clabel'
                        'contourf' : 'colorbar'
                       }
    
    # -------- check inputs, add missing dict fields --------
    valid_keys = ['var', 'plotType', 'plotArgs', 'colorArgs', 'colorFormatter']
    if isinstance(var_dicts, dict): var_dicts = [var_dicts]
    
    for i in range(len(var_dicts)):
        d = var_dicts[i]

        # ignore unrecognized keys
        if(d.keys()[i] not in valid_keys):
            warnings.warn('var_dict key {} not recognized; ignoring'.format(d.keys()[i]))

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
                d['colorFormatter' = color_formatters[pl]]
            else:
                d['colorFormatter'] = None
        else:
            assert isinstance(d['colorFormatter'], str) or d['colorFormatter'] is None, \
                   'colorArgs must be given as a dict, not {}'.format(type(d['colorArgs']))
        
        # merge specified kwargs with defaults
        if(pl in default_args.keys()):
            d['plotArgs'] = {**d['plotArgs'], **default_args[pl]}
        if(d['colorFormatter'] in default_args.keys()):
            d['colorArgs'] = {**d['colorArgs'], **default_args[d['colorFormatter']]]}

    # -------- plot variables --------
    plots = np.empty(len(var_dicts), dtype=object):
    for i in range(len(var_dicts)):
        d = var_dicts[i]
        plotter = getattr(ax, d['plotType'])
        plots[i] = ax.plotter(X, Y, d['var'], **d['plotArgs'])
        
    # -------- format colors --------
    for i in range(len(var)):
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

    if(savefig is not None):
        plt.savefig(savefig, dpi=300)

    if(return_fig)
        return fig
    
    
# -------------------------------------------------------------


