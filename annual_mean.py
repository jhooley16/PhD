import numpy as np
import os
from netCDF4 import Dataset

yr = input('What year? (xxxx): ')

os.chdir('/Users/jmh2g09/Documents/PhD/Data/Gridded/' + yr + '/Geoid')

ssh = np.zeros((361, 30))

for file in os.listdir():
    if file[-8:] == 'geoid.nc':
        
        nc = Dataset(file, 'r')
        
        ssh = ssh + nc.variables['Sea Surface Height'][:]
        lat = nc.variables['Latitude'][:]
        lon = nc.variables['Longitude'][:]
        
        nc.close()

ssh = ssh / 12

nc = Dataset('2012_geoid_mean.nc', 'w', FORMAT='NETCDF4_CLASSIC')
nc.description = yr + 'annual mean'

nc.createDimension('lat', np.size(lat))
nc.createDimension('lon', np.size(lon))

latitudes = nc.createVariable('Latitude', float, ('lat',))
longitudes = nc.createVariable('Longitude', float, ('lon',))
ssh_mean = nc.createVariable('Mean Sea Surface Height', float, ('lon','lat'))

latitudes[:] = lat
longitudes[:] = lon
ssh_mean[:] = ssh

nc.close()