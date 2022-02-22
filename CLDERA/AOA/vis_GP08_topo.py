import matplotlib.pyplot as plt
from metpy.units import units as u
import metpy.constants as const
import numpy as np
import pdb
import cartopy.crs as ccrs

deg2rad = np.pi/180
rad2deg = 1/deg2rad

lat = np.linspace(-90, 90, 100) * deg2rad
lon = np.linspace(0, 360, 100) * deg2rad
LAT, LON = np.meshgrid(lat, lon)

g = const.g
h0 = 3000*u.m
phi0 = 25*deg2rad
phi1 = 60*deg2rad
m = 2

GEO = g*h0*np.sin((LAT-phi0)/(phi1-phi0) * np.pi)**2 * np.cos(m*LON)

mask = np.logical_and(LAT>phi0, LAT<phi1)
GEO = np.ma.filled(np.ma.masked_array(GEO, ~mask), 0)

fig = plt.figure()
ax = fig.add_subplot(111, projection=ccrs.Robinson())
C = ax.contourf(LON*rad2deg, LAT*rad2deg, GEO, levels=10, cmap='rainbow', transform=ccrs.PlateCarree())

cb = fig.colorbar(C, ax=ax, orientation='horizontal') 
cb.set_label(label='Geopotential  [J/kg]', fontsize=12)

ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False, color='k', lw=0.3)

#Cl = ax.contour(LON*rad2deg, LAT*rad2deg, GEO, levels=C.levels, colors='k', transform=ccrs.PlateCarree())
#ax.clabel(Cl, colors=['black'], manual=False, inline=True, fmt=' {:.0f} '.format)

plt.savefig('GP08.png', dpi=300)

