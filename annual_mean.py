import numpy as np
import os
from netCDF4 import Dataset

yr = input('What year? (xxxx): ')

os.chdir('/Users/jmh2g09/Documents/PhD/Data/Gridded/' + yr + '/MDT')

mdt = np.zeros((361, 30))

for file in os.listdir():
    if file[-6:] == 'MDT.nc':
        
        nc = Dataset(file, 'r')
        
        mdt = mdt + nc.variables['mean_dynamic_topography'][:]
        lat = nc.variables['Latitude'][:]
        lon = nc.variables['Longitude'][:]
        
        nc.close()

mdt = mdt / 12

nc = Dataset('2012_mdt_mean.nc', 'w', format='NETCDF4_CLASSIC')

nc.createDimension('lat', np.size(lat))
nc.createDimension('lon', np.size(lon))

latitudes = nc.createVariable('Latitude', float, ('lat',))
longitudes = nc.createVariable('Longitude', float, ('lon',))
mdt_mean = nc.createVariable('annual_mean_dynamic_topography', float, ('lon','lat'))

latitudes[:] = lat
longitudes[:] = lon
mdt_mean[:] = mdt

nc.close()