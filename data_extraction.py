import os
import functions as funct
from netCDF4 import Dataset
import numpy as np

yr = input('What year? (xxxx): ')

for mnth in range(1,13):

    if 0 < mnth < 10:
        os.chdir('/Users/jmh2g09/Documents/PhD/Data/' + yr + '0' + str(mnth) + '_elev')
    else:
        os.chdir('/Users/jmh2g09/Documents/PhD/Data/' + yr + str(mnth) + '_elev')
    directory = os.getcwd()

    month = directory[-7:-5]
    year = directory[-11:-7]
    print(month, year)

    data = funct.month_data(directory)

    lat = data['lat']
    lon = data['lon']
    ssh = data['ssh']
    ice_conc = data['ice_conc']

    # Put the data in a .nc file in /Users/jmh2g09/Documents/PhD/Data/Raw

    os.chdir('/Users/jmh2g09/Documents/PhD/Data/Raw')

    nc = Dataset(year + month + '_raw.nc', 'w', format='NETCDF4')
    nc.description = 'Raw data for ' + month + ', ' + year

    nc.createDimension('lat', np.size(lat))
    nc.createDimension('lon', np.size(lon))
    nc.createDimension('ssh', np.size(ssh))
    nc.createDimension('ice_conc', np.size(ice_conc))

    latitudes = nc.createVariable('Latitude', 'f4', ('lat',))
    longitudes = nc.createVariable('Longitude', 'f4', ('lon',))
    ssh = nc.createVariable('Sea Surface Height', 'f4', ('ssh',))
    ice_conc = nc.createVariable('Sea Ice Concentration', 'f4', ('ice_conc',))

    latitudes[:] = lat
    longitudes[:] = lon
    ssh[:] = ssh
    ice_conc[:] = ice_conc

    nc.close()