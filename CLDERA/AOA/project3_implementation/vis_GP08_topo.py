import matplotlib.pyplot as plt
from metpy.units import units as u
import metpy.constants as const
import numpy as np
import pdb
import cartopy.crs as ccrs

from climate_artist import horizontal_slice as plths

deg2rad = np.pi/180
rad2deg = 1/deg2rad

lat = np.linspace(-90, 90, 100) * deg2rad
lon = np.linspace(0, 360, 100) * deg2rad
LON, LAT = np.meshgrid(lon, lat)

g = const.g
h0 = 3000*u.m
phi0 = 25*deg2rad
phi1 = 60*deg2rad
m = 2

GEO = g*h0*np.sin((np.abs(LAT)-phi0)/(phi1-phi0) * np.pi)**2 * np.cos(m*LON)

mask = np.logical_and(np.abs(LAT)>phi0, np.abs(LAT)<phi1)
GEO = np.ma.filled(np.ma.masked_array(GEO, ~mask), 0)

d = { 'var':GEO, 
      'plotType':'contourf', 
      'colorArgs':{'label':'Geopotential  [m^2/s^2]', 'orientation':'horizontal'}
    }
fig = plths(lon*rad2deg, lat*rad2deg, d, gridlines=True, gridlinesArgs={'rotate_labels':False}, savefig='./GP08.png')
