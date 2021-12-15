import numpy as np
from math import radians
import astropy.units as u
import astropy.constants as const

def great_circle(lat1, lon1, lat2, lon2):
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    return const.R_earth.to(u.m).value * \
           np.arccos((np.sin(lat1) * np.sin(lat2) + np.cos(lat1) * np.cos(lat2) * np.cos(lon1-lon2)))
