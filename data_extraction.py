import os
import functions as funct
from netCDF4 import Dataset
import numpy as np

yr = input('What year? (xxxx): ')

for mnth in range(1,13):

    if 0 < mnth < 10:
        os.chdir('/Users/jmh2g09/Documents/PhD/Data/elev_files/' + yr + '0' + str(mnth) + '_elev')
    else:
        os.chdir('/Users/jmh2g09/Documents/PhD/Data/elev_files/' + yr + str(mnth) + '_elev')
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

    nc = Dataset(year + month + '_raw.nc', 'w', format='NETCDF3_CLASSIC')

    nc.createDimension('station', np.size(lat))

    latitudes = nc.createVariable('lat', float, ('station',))
    longitudes = nc.createVariable('lon', float, ('station',))
    sea_surface_height = nc.createVariable('sea_surface_height', float, ('station',))
#    sea_ice_concentration = nc.createVariable('ice_concentration', float, ('station',))
    crs = nc.createVariable('crs', 'i', ())
    
    latitudes.long_name = 'latitude'
    latitudes.standard_name = 'latitude'
    latitudes.units = 'degrees_north'
    longitudes.long_name = 'longitude'
    longitudes.standard_name = 'longitude'
    longitudes.units = 'degrees_east'
    
    sea_surface_height.standard_name = 'sea_surface_height_above_reference_ellipsoid'
    sea_surface_height.long_name = 'sea_surface_height'
    sea_surface_height.units = 'm'
    sea_surface_height.grid_mapping = 'crs'
    sea_surface_height.tide_system = 'tide_free'
    
#    sea_ice_concentration.standard_name = 'sea_ice_concentration'
#    sea_ice_concentration.long_name = 'sea_ice_concentration'
#    sea_ice_concentration.units = '%'
#    sea_ice_concentration.grid_mapping = 'crs'
    
    crs.grid_mapping_name = 'latitude_longitude'
    crs.semi_major_axis = 6378137.
    crs.inverse_flattening = 298.257222101004
    crs.earth_gravity_constant = 398600500000000.
    crs.earth_rotation_rate = 7.292115e-05

    latitudes[:] = lat
    longitudes[:] = lon
    sea_surface_height[:] = ssh
#    sea_ice_concentration[:] = ice_conc

    nc.close()
    
    # Create .nc file for the sea ice concentration data
    nc = Dataset(year + month + '_icetrack.nc', 'w', format='NETCDF3_CLASSIC')
    nc.createDimension('station', np.size(lat))
    latitudes = nc.createVariable('lat', float, ('station',))
    longitudes = nc.createVariable('lon', float, ('station',))
    sea_ice_concentration = nc.createVariable('ice_concentration', float, ('station',))
    
    latitudes.long_name = 'latitude'
    latitudes.standard_name = 'latitude'
    latitudes.units = 'degrees_north'
    longitudes.long_name = 'longitude'
    longitudes.standard_name = 'longitude'
    longitudes.units = 'degrees_east'
    
    sea_ice_concentration.standard_name = 'sea_ice_concentration'
    sea_ice_concentration.long_name = 'sea_ice_concentration'
    sea_ice_concentration.units = '%'
    sea_ice_concentration.grid_mapping = 'crs'

    latitudes[:] = lat
    longitudes[:] = lon
    sea_ice_concentration[:] = ice_conc

    nc.close()
    
    # In order for the track file to be used with GUT, it needs to be in 
    # netCDF3_classic format.
    nc = Dataset(year + month + '_track.nc', 'w', format='NETCDF3_CLASSIC')

    nc.createDimension('station', np.size(lat))

    longitudes = nc.createVariable('lon', float, ('station',))
    latitudes = nc.createVariable('lat', float, ('station',))
    crs = nc.createVariable('crs', 'i', ())

    latitudes.long_name = 'latitude'
    latitudes.standard_name = 'latitude'
    latitudes.units = 'degrees_north'
    longitudes.long_name = 'longitude'
    longitudes.standard_name = 'longitude'
    longitudes.units = 'degrees_east'
    crs.semi_major_axis = 6378137.
    crs.inverse_flattening = 298.257222101004
    crs.earth_gravity_constant = 398600500000000.
    crs.earth_rotation_rate = 7.292115e-05

    latitudes[:] = lat
    longitudes[:] = lon

    nc.close()
    
    os.system('gut geoidheight_tf -InFile GOCO05s.gfc -T tide-free -InTrack '
    + year + month + '_track.nc -OutFile ' + year + month + '_geoid_track.nc')