import os
import functions as funct
from netCDF4 import Dataset
import numpy as np

year = input('What Year? (xxxx): ')

# Look at the raw data files
os.chdir('/Users/jmh2g09/Documents/PhD/Data/Raw/' + year + '/MDT_track')

# Cycle through each raw file
for file in os.listdir():
    if file[-12:] == 'MDT_track.nc':
        print('Gridding: ' + file)
        nc = Dataset(file, 'r')
    
        lat = nc.variables['lat'][:]
        lon = nc.variables['lon'][:]
        mdt = nc.variables['mean_dynamic_topography'][:]
        ice_conc = nc.variables['sea_ice_concentration'][:]
    
        nc.close()
    
        # Grid mean dynamic topography to a 1-degree grid

        data = funct.grid05(mdt, lon, lat, 1)
        grid_mdt = data['Grid']
        grid_lon = data['Lon']
        grid_lat = data['Lat']
        
        data = funct.grid05(ice_conc, lon, lat, 1)
        grid_ice = data['Grid']
    
        # Put the data in a .nc file in /Users/jmh2g09/Documents/PhD/Data/Gridded

        os.chdir('/Users/jmh2g09/Documents/PhD/Data/Gridded/' + year + '/MDT')
    
        month = file[4:6]
        nc = Dataset(year + month + '_MDT.nc', 'w', format='NETCDF4_CLASSIC')

        nc.createDimension('lat', np.size(grid_lat))
        nc.createDimension('lon', np.size(grid_lon))

        latitudes = nc.createVariable('Latitude', float, ('lat',))
        longitudes = nc.createVariable('Longitude', float, ('lon',))
        gridded_mdt = nc.createVariable('mean_dynamic_topography', float, ('lon','lat'))
        gridded_ice = nc.createVariable('sea_ice_concentration', float, ('lon','lat'))

        latitudes[:] = grid_lat
        longitudes[:] = grid_lon
        gridded_mdt[:] = grid_mdt
        gridded_ice[:] = grid_ice

        nc.close()
        print('Complete!')
        os.chdir('/Users/jmh2g09/Documents/PhD/Data/Raw/' + year + '/MDT_track')