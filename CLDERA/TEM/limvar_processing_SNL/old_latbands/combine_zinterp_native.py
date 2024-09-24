'''
Joe Hollowed 2024

Providing a function for combining vertically interpolated and native limvar data.
When we perform the vertical interpolation from hybridg model levels to pure pressure
levels, we do this only on levels for which hybm > 0 to minimize needed storage space
for the interpolated datasets. This is because for levels where hybm = 0 (everything 
above ~170 hPa in E3SM), the vertical corodinate is already an isobar. The function
provided here takes an interpolated data file, searches for the associated native
dataset by the "parent_dataset" attribute, concatenates the two datasets, and returns
it as an xarray Dataset.
'''

import pdb
import xarray as xr

def combine_interp_native_data(interp_file):

    interp = xr.open_dataset(interp_file)
    native = xr.open_dataset(interp.attrs['parent_dataset'])

    # select the native data in the region where hybm = 0
    shape = native['U'].shape
    hybm_mask = native.hybm == 0
    native = native.isel(lev=(hybm_mask))

    # rename native lev -> plev, since these are isobars
    # adopt attrs from interp
    native = native.rename({'lev':'plev'})
    native.attrs = interp.attrs

    # select variables from the native dataset that were interpolated
    native      = native[list(interp.data_vars)]

    # concat in the vertical, overwriting the Dataset interp, return
    interp = xr.concat([native, interp], dim='plev', data_vars='minimal', coords='minimal')
    return interp


# -------------------------------------------------------------------


# for testing
if __name__ == '__main__':
    data = '/gpfs/cldera/data/e3sm/e3sm-cldera/CLDERA-historical/limvar_src_tag/tem_processed/interp_to_pressure_3D/v2.LR.WCYCL20TR.pmcpu.ctools.lv.3tagso4.ens4.cf.eam.h0.1991-06_pinterp.nc'
    combine_interp_native_data(data)
