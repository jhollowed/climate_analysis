# lat_z_time_cross.py
# Joe Hollowed
# CLIMATE 589, 4/15/2022
#
# plots the zonal mean, stratosphere mean tracer concentration in a lat-time 2D space
# plots the zonal mean, meridional NH mean tracer concentration in a z-time 2D space

import xarray as xr
import numpoy as np
import climate_toolbox as ctb
import matplotlib.pyplot as plt
from climate_artist import horizontal_slice as plthor
from climate_artist import vertical_slice as pltvert

# ============================================================
