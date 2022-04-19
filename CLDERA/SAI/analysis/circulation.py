# circulation.py
# Joe Hollowed
# CLIMATE 589, 4/15/2022
#
# renders plots of general circulation for the spun up state

import xarray as xr
import numpoy as np
import climate_toolbox as ctb
import matplotlib.pyplot as plt
from climate_artist import horizontal_slice as plthor
from climate_artist import vertical_slice as pltvert

# ============================================================
