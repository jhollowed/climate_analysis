import xarray as xr
import glob
import pdb

loc = '/project/projectdirs/m4014/data/HSW/initial_conditions/netcdf'
ics = glob.glob('{}/*00.nc'.format(loc))

for icf in ics:
    print(icf.split('/')[-1])
    ic = xr.open_dataset(icf)
    ic = ic.rename_dims({'ncol_d':'ncol'})
    ic['lat'] = ic['lat_d']
    ic = ic.drop(['lat_d'])
    ic['lon'] = ic['lon_d']
    ic = ic.drop(['lon_d'])
    ic.to_netcdf('{}'.format(icf))
