import pdb
import numpy as np
import xarray as xr
from datetime import datetime
from matplotlib import colors
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter
from matplotlib.ticker import FuncFormatter
from matplotlib.ticker import ScalarFormatter

# local imports
import compute_ensemble_stats as ces

# --- constants
g = 9.80665   # global average of gravity at mean sea level in m/s^2
H = 7000      # scale height
a = 6.37123e6 # radius of Earth in m


# =======================================================================
# =======================================================================


def get_daily_variable(var, q=None, pmin=None, pmax=None, latmin=None, latmax=None, 
                       tmin=None, tmax=None, year=None, month=None, overwrite=False, debug=False,
                       skip_nosrctag=False):
    '''
    returns a dictionary of the ensemble data, TEM data, TEM budget data, as well as those dataset's 
    counterfactuals, impacts, p-values, and coherence for a specfied variable with a specified 
    slicing. This variable can be stored in either the ensemble mean, TEM, or TEM budget data. They 
    will all be searched for a variable with the matching name. 
    This function calls compute_ensemble_stats.get_10daily_stats() in order to get the processed data. 
    If the data has not yet been processed for this choice of slicing, it will be computed. If it has, 
    it will simply be read and returned.
    The slicing uses plotting_utils.do_slicing() (see docs therein), with the same input arguments, 
    but with average=True. Meaning, the data will be averaged over the specified slice before the 
    stats are computed.

    Parameters
    ----------
    var : str
        the variable name
    q : str, optional
        tracer name, either 'aoa', 'e90', or None (the default), in which case
        the specified variable is not a traer quantity
    debug : bool, optional
        if True, print out available variables for each dataset. Defaults to False.
    overwrite : bool, optional
        whether or not to overwrite the data produced by a previous run of this function.
        Defaults to False, in which case the processed data is simply read and returned 
        if it already exists.
    pmin : float, optional
        minimum pressure (highest vertical level)
    pmax : float, optional
        maximum pressure (lowest vertical level)
    latmin : float, optional
        minimum latitude 
    latmax : float, optional
        maximum latitude
    tmin : cftime._cftime.Datetime, optional
        minimum time
    tmax : cftime._cftime.Datetime, optional
        maximum time
    year : int, optional
        year to slice
    month : int, optional
        month to slice
    '''
    
    # ---- get tracer index
    qi = {None:0, 'aoa':1, 'e90':2}
    qi = qi[q]
    
    # ---- get processing args
    args={'overwrite':overwrite, 'pmin':pmin, 'pmax':pmax, 'latmin':latmin, 'latmax':latmax, 
          'tmin':tmin, 'tmax':tmax, 'year':year, 'month':month, 'skip_nosrctag':skip_nosrctag}
    
    # ---- begin
    #print('getting data for variable {}...'.format(var))
    
    # ---- ensemble data
    data, cf, impact, pval, coherence = ces.get_10daily_stats(dataset='ens', qi=0, **args)
    if(debug):
        print('available data vars: {}'.format(qi, list(data.data_vars)))
    try:
        data, cf, impact, pval, coherence = [data[var], cf[var], impact[var],
                                             pval[var], coherence[var]]
        return {'ensmean':data, 'cfmean':cf, 'impact':impact, 'pval':pval, 'coherence':coherence}
    except KeyError: pass
    
    # ---- tem data
    tem_data, tem_cf, tem_impact, tem_pval, tem_coherence = \
           ces.get_10daily_stats(dataset='tem', qi=qi, **args)
    if(debug):
        print('available TEM vars for tracer {}: {}'.format(qi, list(tem_data.data_vars)))
    try:
        data, cf, impact, pval, coherence = [tem_data[var], tem_cf[var], tem_impact[var],
                                             tem_pval[var], tem_coherence[var]]
        return {'ensmean':data, 'cfmean':cf, 'impact':impact, 'pval':pval, 'coherence':coherence}
    except KeyError: pass
    
    # ---- budget data
    budget_data, budget_cf, budget_impact, budget_pval, budget_coherence = \
                        ces.get_10daily_stats(dataset='budget', qi=qi, **args)
    if(debug):
        print('available TEM budget vars for tracer{}: {}'.format(qi, list(budget_data.data_vars)))
    try:
        data, cf, impact, pval, coherence = [budget_data[var], budget_cf[var], budget_impact[var],
                                             budget_pval[var], budget_coherence[var]]
        return {'ensmean':data, 'cfmean':cf, 'impact':impact, 'pval':pval, 'coherence':coherence}
    except KeyError: pass

    # ---- if variable hasn't been found by here, it doesn't exist
    assert isinstance(data, xr.core.dataarray.DataArray), 'variable {} not found!'.format(var)

    # get tropopause data
    #data_tropp, cf_tropp, impact_tropp = [data['TROP_P']/100, cf['TROP_P']/100, impact['TROP_P']/100]

# -----------------------------------------------------------------------

def do_slicing(data, pmin=None, pmax=None, latmin=None, latmax=None, 
               tmin=None, tmax=None, month=None, year=None, average=False):
    '''
    Slices the data optionally on pressure, latitude, and time. The time slicing can be done
    either by index, or by month/year. If (tmin, tmax) are not passed, then year and month must 
    be passed, and vice-versa.

    Parameters
    ----------
    data : xarray DataArray, or list
        The data to slice. If list, each item must be an xarray DataArray, and the slicing is
        applied to each
    pmin : float, optional
        minimum pressure (highest vertical level)
    pmax : float, optional
        maximum pressure (lowest vertical level)
    latmin : float, optional
        minimum latitude 
    latmax : float, optional
        maximum latitude
    tmin : cftime._cftime.Datetime, optional
        minimum time
    tmax : cftime._cftime.Datetime, optional
        maximum time
    year : int, optional
        year to slice
    month : int, optional
        month to slice
    average : bool, optional
        whether or not to average the data over the specified slice. Defaults to False
    '''

    # --- check inputs
    assert (np.array([pmin, pmax]) == None).sum() != 1,     "must pass both or neither or pmin, pmax"
    assert (np.array([latmin, latmax]) == None).sum() != 1, "must pass both or neither or pmin, pmax"
    assert (np.array([tmin, tmax]) == None).sum() != 1,     "must pass both or neither or pmin, pmax"
    assert (np.array([tmin, tmax, year]) == None).sum() > 0, \
                                                     "must pass either (tmin,tmax) or year, not both"
    if(year is None): assert month is None,                  "year must be passed if month is passed"
    assert(type(data) == type(xr.DataArray()) or type(data) == type(xr.Dataset()) or type(data) == list), \
                                                "data must be a DataArray, Dataset, or list of those"
    if(type(data) == list): assert np.sum([type(d) for d in data]) == len(data),\
                                                "if data is a list, each element must be a DataArray"
    else:
        data = [data]
    if(year is not None):
        assert(type(year)==int), "year must be an int"
    if(month is not None):
        assert(type(month)==int or type(month)==list), "month must be an int or list of int"
        if(type(month) == int): 
            month = [month]

    # --- do slicing and averaging
    for i in range(len(data)):
       
        # pressure
        if(pmin is not None):
            data[i] = data[i].sel(plev = slice(pmin, pmax))
            if(average):
                # weight by the pressure thickness of each level
                weights = np.gradient(data[i].plev)
                weights.name = "weights"
                data[i] = data[i].weighted(weights)
                data[i] = data[i].mean(("plev"))
        # latitude
        if(latmin is not None): 
            data[i] = data[i].sel(lat = slice(latmin, latmax))
            if(average):
                # weight by cos(lat)
                weights = np.cos(np.deg2rad(data[i].lat))
                weights.name = "weights"
                data[i] = data[i].weighted(weights)
                data[i] = data[i].mean(("lat"))
        # time
        if(tmin is not None):   
            data[i] = data[i].sel(time = slice(tmin, tmax))
            if(average):
                data[i] = data[i].mean('time')
        if(year is not None):
            data[i] = data[i].where(data[i]['time.year']==year, drop=True)
            if(month is not None):
                data[i] = data[i].where(data[i]['time.month']==month, drop=True)
            if(average):
                data[i] = data[i].mean('time')
    
    # if a list was passed, return a list. Otherwise return a DataArray or Dataset
    if(len(data) == 1): return data[0]
    else: return data

# -----------------------------------------------------------------------

def shift_integrated_tendency(x, dxdti, icdate):
    '''
    Shifts an integrated tendencies in time, given an initial condition date
    
    Parameters
    ----------
    x : xarray DataArray
        data to sample initial condition from
    dxdti : xarray DataArray
        integrated time tendency of the data x
    icdate : cftime._cftime.Datetime
        date to use as initial condition
    '''
    targs = {'time':icdate, 'method':'nearest'}
    x0 = x.sel(**targs)
    x_integrated = x0 + dxdti - dxdt.sel(**targs)
    # return only the data after the initial condition date
    x_integrated = x_integrated.sel(time=slice(icdate, None))
    return x_integrated

# -----------------------------------------------------------------------

def scale_EP_flux_vectors(fig, ax, epfy, epfz, dslat=None, dsp=None, p0=100000):
    '''
    scales the EP flux vectors for graphical display, according to Jucker (2021), 
    their Table 2. The scaling of the vectors depends on the aspect ratio of the
    figure, and so this must not change after these scalings are computed, else
    the displayed vectors will be incorrect. The following important assumptions
    are made:
     - it is assumed that the limits of the data coordiantes are identical to
       the limits of the coordinate axes in the plot. That is, any slicing of 
       the data must be done FIRST
     - it is assumed that pressure decreases upward on the y-axis
     - it is assumed that latitude increases to the right on the x-axis 

    Parameters
    ----------
    fig : matplotlib figure object
        the figure on which the plotting axis is located.
    ax : matplotlib axis object
        The axis on which the vectors will be plotted.
    epfy : xarray DataArray
        the northward component of the EP flux in log-pressure cordinate
        (epfy from PyTEMDiags, or Eq. A13 from Gerber+Manzini 2016).
        Must have coordinates 'plev' and 'lat'
    epfz : xarray DataArray
        the upward component of the EP flux in log-pressure cordinate
        (epfz from PyTEMDiags, or Eq. A14 from Gerber+Manzini 2016).
        Must have coordinates 'plev' and 'lat'
    dslat : float, optional
        latitude downsampling factor. If None (the default), no downsampling 
        is done. Otherwise, the data will be downsampled in latitude by dslat
        times, starting at the dslat element. That is, if dslat=4, then the
        data is subsampled via slice(4, None, 4)
    dsp : float, optional
        pressure downsampling factor. If None (the default), no downsampling 
        is done. Otherwise, do the following:
        For a linear pressure axis,the data will be downsampled in pressure 
        by dsp times, starting at the dslat element. That is, if dsp=4, then 
        the data is subsampled via slice(4, None, 4)
        For a logarithmic pressure axis, the data will be downsampled to a
        total of nlev/dsp points, with indices logarithmically spaced from 0
        to nlev (where nlev is the number of vertical levels in the data)
    p0 : float, optional
        reference pressure, in Pa. Defaults to 1000 hPa
    '''

    # get display dimensions
    # in matplotlib, we need get get the figure size, and then
    # get the portion of that figure that is covered by the axes
    fig_width, fig_height    = fig.get_size_inches()
    _, _, ax_wsize, ax_hsize = ax.get_position().bounds
    X, Y = fig_width * ax_wsize, fig_height * ax_hsize
    
    # get vertical scale type
    pscale = ax.get_yscale()
    
    # downsample if requested
    if(dslat is not None):
        epfy = epfy.isel(lat=slice(dslat, None, dslat))
        epfz = epfz.isel(lat=slice(dslat, None, dslat))
    if(dsp):
        if(pscale == 'linear'):
            epfy = epfy.isel(plev=slice(dsp, None, dsp))
            epfz = epfz.isel(plev=slice(dsp, None, dsp))
        elif(pscale == 'log'):
            idx = np.logspace(0, np.log10(len(epfy.plev)-dsp), num=int(len(epfy.plev)/dsp), base=10)
            idx = np.unique(np.round(idx).astype(int))
            epfy = epfy.isel(plev=idx)
    
    # get dimensions
    plev       = epfy.plev
    lat        = epfy.lat
    acoslat    = a*np.cos(np.deg2rad(lat))
    p0, p1     = plev.max(), plev.min()
    lat0, lat1 = lat.min(), lat.max()

    # transform from log-pressure to pressure coordinate
    Fphi = epfy * p0/plev
    Fp   = epfz * -p0/H

    # get EP flux components, i.e. Jucker Eq.1-2
    fphi = Fphi / acoslat
    fp   = Fp   / acoslat

    # apply Edmon (1980) scaling, i.e. Jucker Eq. 4
    hFphi = 2*np.pi/g * acoslat**2 * (fphi)
    hFp   = 2*np.pi/g * acoslat**2 * (a * fp)

    # set the type of pressure scaling
    if(pscale == 'linear'): pdiff = p1-p0
    elif(pscale == 'log'):  pdiff = plev * np.log(p1/p0)
    
    # do scaling
    Fx = hFphi * (X/Y) / ((lat1-lat0) * np.pi/180)
    Fy = hFp * 1/pdiff
    
    return Fx, Fy

# -----------------------------------------------------------------------

def scale_resvel_vectors(fig, ax, vtem, wtem, dslat=None, dsp=None, p0=100000, norm=False):
    '''
    scales residual velocity vectors for graphical display 

    Parameters
    ----------
    '''

    # get display dimensions
    # in matplotlib, we need get get the figure size, and then
    # get the portion of that figure that is covered by the axes
    fig_width, fig_height    = fig.get_size_inches()
    _, _, ax_wsize, ax_hsize = ax.get_position().bounds
    X, Y = fig_width * ax_wsize, fig_height * ax_hsize
    
    # get vertical scale type
    pscale = ax.get_yscale()
    
    # downsample if requested
    if(dslat is not None):
        epfy = epfy.isel(lat=slice(dslat, None, dslat))
        epfz = epfz.isel(lat=slice(dslat, None, dslat))
    if(dsp):
        if(pscale == 'linear'):
            epfy = epfy.isel(plev=slice(dsp, None, dsp))
            epfz = epfz.isel(plev=slice(dsp, None, dsp))
        elif(pscale == 'log'):
            idx = np.logspace(0, np.log10(len(epfy.plev)-dsp), num=int(len(epfy.plev)/dsp), base=10)
            idx = np.unique(np.round(idx).astype(int))
            epfy = epfy.isel(plev=idx)
    
    # get dimensions
    plev       = epfy.plev
    lat        = epfy.lat
    acoslat    = a*np.cos(np.deg2rad(lat))
    p0, p1     = plev.max(), plev.min()
    lat0, lat1 = lat.min(), lat.max()

    # transform from log-pressure to pressure coordinate
    Fphi = epfy * p0/plev
    Fp   = epfz * -p0/H

    # get EP flux components, i.e. Jucker Eq.1-2
    fphi = Fphi / acoslat
    fp   = Fp   / acoslat

    # apply Edmon (1980) scaling, i.e. Jucker Eq. 4
    hFphi = 2*np.pi/g * acoslat**2 * (fphi)
    hFp   = 2*np.pi/g * acoslat**2 * (a * fp)

    # set the type of pressure scaling
    if(pscale == 'linear'): pdiff = p1-p0
    elif(pscale == 'log'):  pdiff = plev * np.log(p1/p0)
    
    # do scaling
    Fx = hFphi * (X/Y) / ((lat1-lat0) * np.pi/180)
    Fy = hFp * 1/pdiff
    
    # do normalization
    if(norm):
        pass
    
    return Fx, Fy

# -----------------------------------------------------------------------

def to_datetime(times):
    '''
    Converts an DataArray of cftime.datetime objects to datetime objects, which
    matplotlib prefers
    
    Parameters
    ----------
    times : xarray DataArray of cftime.datetime objects
        the times to convert
    
    Returns
    -------
    datetimes : array of datetime objects
        the converted times
    '''
    datetimes = [datetime(date.year, date.month, date.day) for date in times.values]
    return datetimes

# -----------------------------------------------------------------------

def season_timeticks(taxis, times, option):
    '''
    Changes the ticks on an input time axis so that they are drawn at specified months
    
    Parameters
    ----------
    taxis : matplotlib.axis.XAxis or  matplotlib.axis.YAxis object
        the time axis to modify
    times : array-like of cftime.datetime or datetime objects
        the time coordiantes of the data being plotted
    option : str
        Which ticks to draw. Options are:
        - 'month'    : draw at DJF MAM JJA SON (every month)
        - 'season'   : draw at J A J O         (center month of every season)
        - 'solstice' : draw at J J             (center month of summer, winter)
        - 'equinox'  : draw at A O             (center month of spring, fall)
    '''
    
    months = {'month'   : [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], 
              'season'  : [1, 4, 7, 10],
              'solstice': [1, 7],
              'equinox' : [4, 10]}
    ticks = []
    for year in np.unique([t.year for t in times]):
        for month in months[option]:
            ticks.append(datetime(year, month, 1))
    taxis.set_ticks(ticks)
    
    date_format = DateFormatter("%b %-d")
    taxis.set_major_formatter(date_format)
    
# -----------------------------------------------------------------------
    
def format_paxis(ax, ticks=None):
    '''
    Formats an the pressure axis of an input axis for inverted, log-scaled, scalar formatted ticks
    
    Parameters
    ----------
    ax : matplotlib.axis object
        the axes to modify. It is assumed that the pressure axis is the y-axis
    '''
    
    if(ticks is None):
        ticks=[1000, 300, 100, 30, 10, 3, 1, 0.3, 0.1, 0.03, 0.01]
        tick_str = []
        
    lim = sorted(ax.get_ylim())[::-1]
    ax.invert_yaxis()
    ax.set_yscale('log')
    ax.yaxis.set_major_formatter(ScalarFormatter())
    ax.yaxis.set_ticks(ticks)
    ax.set_ylim(lim)
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: \
                              f'{x:.3f}'.rstrip('0').rstrip('.')))

# -----------------------------------------------------------------------

def cbarfmt(v, pos):
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
def pretty_format_colorbar():
    '''
    Returns a colorbar FuncFormatter, which can be used as the "format" input
    argument to matplotlib.figure.Figure.colorbar, which formats all values grater 
    than 0.01 as decimal values with two decimal digits, with trailing zeros removed.
    All values less than 0.01 are formatted in scientific notation, with leading 
    zeros on the exponent removed, and trailing zeros on the value removed.

    Parameters
    ----------
    fmt : str
        string format for values. If not supplied, 
    '''
    return(FuncFormatter(cbarfmt))

# -----------------------------------------------------------------------

def get_cmap_norm(levels, norm):
    '''
    Construct matplotlib color normalization objects from keywords

    Parameters
    ----------
    levels : numpy array or xarray DataArray
        the contour levels
    norm : str
        string specifying the type of norm. Options are:
        'linear    : linear colormap
        'twoslope' : separate linear color scales are used on each side 
                     of zero, to utilize the full range 
        'uneven'   : the levels are non-equidistant, and colors should 
                     be sampled evenly for each level, reagrdless of the 
                     separation between levels
        'log'      : logarithmically space the colors
        'symlog'   : logarithmically space the colors symmetrically on 
                     each side of zero
    '''
    
    if(norm == 'linear' or norm == 'twoslope'):
        low, high = min(levels), max(levels)
        if(low < 0 and high > 0):
            norm = colors.TwoSlopeNorm(vmin=min(levels), vcenter=0, vmax=max(levels))
        else:
            norm = None
    if(norm == 'log'):
        norm = colors.LogNorm(vmin=min(levels), vmax=max(levels))
    if(norm == 'symlog'):
        norm = colors.SymLogNorm(linthresh=min(levels[levels>0]), linscale=1)
    if(norm == 'uneven'):
        norm = colors.BoundaryNorm(levels, plt.get_cmap('viridis').N, extend='both')
    return norm
    
# -----------------------------------------------------------------------

def auto_levels(data, pp=0.5, get_norm=False):
    '''
    Automatically choose levels for contour plots

    Parameters
    ----------
    data : DataArray
        the data
    pp : float, optional
        the percentile to use for configuring the levels, which will span
        uniformly from (pp percentile) to zero, and zero to (100-pp percentile)
        Default is 0.5
    norm : bool, optional
        whether or not to construct and return a matplotlib color normalizer object.
        If the data does not include zero, this is None. If the data does include
        zero, this is a TwoSlopeNorm.
        Defaults to False.
    
    Returns
    ------
    levels : list
        the automatically chosen levels
    norm : numpy colors norm object, optional
        the colormap normalization associated with the levels
    '''
    low, high = np.percentile(np.ravel(data), pp), np.percentile(np.ravel(data), 100-pp)
    if(low < 0 and high > 0):
        neglevels = np.linspace(low, 0, 6)
        poslevels = np.linspace(0, high, 6)[1:]
        levels    = np.hstack([neglevels, poslevels])
        norm      = colors.TwoSlopeNorm(vmin=min(levels), vcenter=0, vmax=max(levels))
    else:
        levels = np.linspace(np.percentile(np.ravel(data), pp), 
                                           np.percentile(np.ravel(data), 100-pp), 12)
        norm = None
    if(get_norm):
        return levels, norm
    else:
        return levels
    