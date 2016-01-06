import os
import functions as funct
from netCDF4 import Dataset
import numpy as np

mnth = input('What month? (xx): ')
yr = input('What year? (xxxx): ')

os.chdir('/Users/jmh2g09/Documents/PhD/Data/' + yr + mnth + '_elev')
directory = os.getcwd()

month = directory[-7:-5]
year = directory[-11:-7]
print(month, year)

data = funct.month_data(directory)

lat = data['lat']
lon = data['lon']
ssh = data['ssh']
ice_conc = data['ice_conc']

# Grid sea surface height to a 1-degree grid

data = funct.grid(ssh, lon, lat, 1)
grid_ssh = data['Grid']
grid_lon = data['Lon']
grid_lat = data['Lat']

# Put the data in a .nc file in /Users/jmh2g09/Documents/PhD/Data/Gridded

os.chdir('/Users/jmh2g09/Documents/PhD/Data/Gridded')

netfile = Dataset(year + month + '.nc', 'w', format='NETCDF4')
netfile.description = 'Data for ' + month + ', ' + year

netfile.createDimension('lat', np.size(grid_lat))
netfile.createDimension('lon', np.size(grid_lon))

latitudes = netfile.createVariable('latitude', 'f4', ('lat',))
longitudes = netfile.createVariable('longitude', 'f4', ('lon',))
gridded_data = netfile.createVariable('ssh', 'f4', ('lon','lat'))

latitudes[:] = grid_lat
longitudes[:] = grid_lon
gridded_data[:] = grid_ssh

netfile.close()