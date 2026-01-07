import pdb
import glob
import numpy as np
import pandas as pd
import xarray as xr
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

data_dir = '/work/noaa/co2/aschuh/WOMBAT_stuff'
inverse_dir = f'{data_dir}/wombat-v3-inverse'
forward_dir = f'{data_dir}/wombat-v3-forward'
gc_transport_dir = f'{forward_dir}/3a_transport_gc/intermediates/runs'

mapping = pd.read_csv(f'{gc_transport_dir}/mapping.csv')

source = 'bio_gpp'
region = 10
split = 1

fname = lambda component: f'{source}_{component}_pft05_regionRegion{region}'
comp_labels = ['intercept', 'trend', 'sin1', 'sin2', 'sin3', 'cos1', 'cos2', 'cos3']
comp_names  = dict(zip(comp_labels, ['intercept', 'trend', 
                                     'sin12_1', 'sin12_2', 'sin12_3', 
                                     'cos12_1', 'cos12_2', 'cos12_3']))

comp_maps = dict(zip(comp_labels, [mapping[mapping['basis_function'] == fname(comp_names[comp])] 
                               for comp in comp_labels]))
comp_runs = dict(zip(comp_labels, [comp_maps[comp]['run'].iat[0] for comp in comp_labels]))
comp_species = dict(zip(comp_labels, [comp_maps[comp]['species'].iat[0] for comp in comp_labels]))


comp_data = {}
for i,comp in enumerate(comp_labels):
    if(i>0):continue
    data_files = sorted(glob.glob(f'{gc_transport_dir}/{comp_runs[comp]}_split0{split}'\
                                   '/OutputDir/GEOSChem.SpeciesConcThreeHourly*'))
    print(f'reading {comp} files')
    data = [0]*len(data_files)
    for i in range(len(data_files)):
        data_3hr = xr.open_dataset(data_files[i])[f'SpeciesConcVV_{comp_species[comp]}']
        data[i]  = data_3hr.coarsen(time=8, boundary='trim').mean()
        print(f'{i+1}/{len(data_files)}', end=['\r','\n'][i==len(data_files)-1])
    print(f'concatenating')
    data = xr.concat(data, dim='time')
    comp_data[comp] = data

time = comp_data['intercept'].time
lat  = comp_data['intercept'].lat
lon  = comp_data['intercept'].lon
lev  = comp_data['intercept'].lev

fig = plt.figure()
ax = fig.add_subplot(111, projection=ccrs.Mollweide())
ax.contourf(lon, lat, comp_data['intercept'].isel(time=0, lev=0), transform=ccrs.PlateCarree(), cmap='viridis')
ax.coastlines()
ax.set_global()
plt.show()

pdb.set_trace()


