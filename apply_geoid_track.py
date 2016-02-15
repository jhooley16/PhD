import os
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as pl
from mpl_toolkits.basemap import Basemap

year = input('What year? (xxxx): ')

os.chdir('/Users/jmh2g09/Documents/PhD/Data/Raw/' + year)

# Run through the months and calculate the MDT_track for each month
for imnth in range(1,13):
    if imnth < 10:
        month = '0' + str(imnth)
    if imnth >= 10:
        month = str(imnth)
    print(month)

    for file in os.listdir():
        if file[4:6] == month:
            if file[-14:] == 'geoid_track.nc':
                print(file)
                nc_geoid = Dataset(file, 'r')
                
                lat = nc_geoid.variables['lat'][:]
                lon = nc_geoid.variables['lon'][:]
                geoid_height = nc_geoid.variables['geoid_height'][:]
                
                nc_geoid.close()
            if file[-6:] == 'raw.nc':
                print(file)
                nc_ssh = Dataset(file, 'r')

                sea_surface_height = nc_ssh.variables['sea_surface_height'][:]
                
                nc_ssh.close()
            if file[-11:] == 'icetrack.nc':
                print(file)
                nc_ice = Dataset(file, 'r')
                
                sea_ice_conc = nc_ice.variables['ice_concentration'][:]
                
                nc_ice.close()

    mdt = sea_surface_height - geoid_height
    
    nc_mdt = Dataset('MDT_track/' + year + month + '_MDT_track.nc', 'w', format='NETCDF3_CLASSIC')
    
    nc_mdt.createDimension('station', np.size(lat))
    
    latitudes = nc_mdt.createVariable('lat', float, ('station',))
    longitudes = nc_mdt.createVariable('lon', float, ('station',))
    mean_dynamic_topography = nc_mdt.createVariable('mean_dynamic_topography', float, ('station',))
    sea_ice_concentration = nc_mdt.createVariable('sea_ice_concentration', float, ('station',))


    latitudes.long_name = 'latitude'
    latitudes.standard_name = 'latitude'
    latitudes.units = 'degrees_north'
    longitudes.long_name = 'longitude'
    longitudes.standard_name = 'longitude'
    longitudes.units = 'degrees_east'
    mean_dynamic_topography.long_name = 'mean_dynamic_topography'
    mean_dynamic_topography.standard_name = 'sea_surface_height_above_WGS84_geoid'
    mean_dynamic_topography.units = 'm'
    sea_ice_concentration.long_name = 'sea_ice_concentration'
    sea_ice_concentration.standard_name = 'sea_ice_concentration'
    sea_ice_concentration.units = '%'
    
    latitudes[:] = lat
    longitudes[:] = lon
    mean_dynamic_topography[:] = mdt
    sea_ice_concentration[:] = sea_ice_conc
    
    nc_mdt.close()