import os
import functions as funct
from netCDF4 import Dataset
import numpy as np

year = input('What Year? (XXXX): ')

# Look at the raw data files
os.chdir('/Users/jmh2g09/Documents/PhD/Data/Raw/' + year)

# Cycle through each raw file
for file in os.listdir():
    print('Gridding: ' + file)
    nc = Dataset(file, 'r')
    
    lat = nc.variables['Latitude'][:]#.tolist()
    lon = nc.variables['Longitude'][:]#.tolist()
    ssh = nc.variables['Sea Surface Height'][:]#.tolist()
    ice_conc = nc.variables['Sea Ice Concentration'][:]#.tolist()
    
    nc.close()
    
    # Grid sea surface height to a 1-degree grid

    data = funct.grid(ssh, lon, lat, 1)
    grid_ssh = data['Grid']
    grid_lon = data['Lon']
    grid_lat = data['Lat']
    print('ssh data gridded...')
    
    data = funct.grid(ice_conc, lon, lat, 1)
    grid_ice_conc = data['Grid']
    print('sea ice concentration gridded...')
    
    # Put the data in a .nc file in /Users/jmh2g09/Documents/PhD/Data/Gridded

    os.chdir('/Users/jmh2g09/Documents/PhD/Data/Gridded/' + year)
    
    month = file[-9:-7]
    nc = Dataset(year + month + '_nogeoid.nc', 'w', FORMAT='NETCDF4_CLASSIC')
    nc.description = 'Data for ' + month + ', ' + year

    nc.createDimension('lat', np.size(grid_lat))
    nc.createDimension('lon', np.size(grid_lon))

    latitudes = nc.createVariable('Latitude', float, ('lat',))
    longitudes = nc.createVariable('Longitude', float, ('lon',))
    gridded_ssh = nc.createVariable('Sea Surface Height', float, ('lon','lat'))
    gridded_ice_conc = nc.createVariable('Sea Ice Concentration', float, ('lon','lat'))

    latitudes[:] = grid_lat
    longitudes[:] = grid_lon
    gridded_ssh[:] = grid_ssh
    gridded_ice_conc[:] = grid_ice_conc

    nc.close()
    print('Complete!')
    os.chdir('/Users/jmh2g09/Documents/PhD/Data/Raw/' + year)