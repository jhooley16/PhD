import os
import functions as funct
from netCDF4 import Dataset
import numpy as np

# Look at the raw data files
os.chdir('/Users/jmh2g09/Documents/PhD/Data/Raw')

# Cycle through each raw file
for file in os.listdir():
    
    nc = Dataset(file, 'r')
    
    lat = nc.variables['Latitude'][:]
    lon = nc.variables['Longitude'][:]
    ssh = nc.variables['Sea Surface Height'][:]
    ice_conc = nc.variables['Sea Ice Concentration'][:]
    
    nc.close()
    
    # Grid sea surface height to a 1-degree grid

    data = funct.grid(ssh, lon, lat, 1)
    grid_ssh = data['Grid']
    grid_lon = data['Lon']
    grid_lat = data['Lat']
    
    data = funct.grid(ice_conc, lon, lat, 1)
    grid_ice_conc = data['Grid']

    # Put the data in a .nc file in /Users/jmh2g09/Documents/PhD/Data/Gridded

    os.chdir('/Users/jmh2g09/Documents/PhD/Data/Gridded')

    nc = Dataset(year + month + '_nogeoid.nc', 'w', format='NETCDF4')
    nc.description = 'Data for ' + month + ', ' + year

    nc.createDimension('lat', np.size(grid_lat))
    nc.createDimension('lon', np.size(grid_lon))

    latitudes = nc.createVariable('Latitude', 'f4', ('lat',))
    longitudes = nc.createVariable('Longitude', 'f4', ('lon',))
    gridded_ssh = nc.createVariable('Sea Surface Height', 'f4', ('lon','lat'))
    gridded_ice_conc = nc.createVariable('Sea Ice Concentration', 'f4', ('lon','lat'))

    latitudes[:] = grid_lat
    longitudes[:] = grid_lon
    gridded_ssh[:] = grid_ssh
    gridded_ice_conc[:] = grid_ice_conc

    nc.close()