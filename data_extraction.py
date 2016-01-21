import os
import functions as funct
from netCDF4 import Dataset
import numpy as np

yr = input('What year? (xxxx): ')

for mnth in range(1,13):

    if 0 < mnth < 10:
        os.chdir('/Users/jmh2g09/Documents/PhD/Data/elev_files' + yr + '0' + str(mnth) + '_elev')
    else:
        os.chdir('/Users/jmh2g09/Documents/PhD/Data/elev_files' + yr + str(mnth) + '_elev')
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

    os.chdir('/Users/jmh2g09/Documents/PhD/Data/Raw/' + year)

    nc = Dataset(year + month + '_raw.nc', 'w', FORMAT='NETCDF4_CLASSIC')
    nc.description = 'Raw data for ' + month + ', ' + year

    nc.createDimension('lat', np.size(lat))
    nc.createDimension('lon', np.size(lon))
    nc.createDimension('ssh', np.size(ssh))
    nc.createDimension('ice_conc', np.size(ice_conc))

    latitudes = nc.createVariable('Latitude', float, ('lat',))
    longitudes = nc.createVariable('Longitude', float, ('lon',))
    sea_surface_height = nc.createVariable('Sea Surface Height', float, ('ssh',))
    sea_ice_concentration = nc.createVariable('Sea Ice Concentration', float, ('ice_conc',))

    latitudes[:] = lat
    longitudes[:] = lon
    sea_surface_height[:] = ssh
    sea_ice_concentration[:] = ice_conc

    nc.close()