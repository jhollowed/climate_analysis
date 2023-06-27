import xarray as xr
import glob
import pdb

file = glob.glob('./*regrid*')[0]
out = 'time_avg.nc'

df = xr.open_dataset(file)
tdf = df.mean('time')
pdb.set_trace()
tdf.to_netcdf(out)
