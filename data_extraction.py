import os
import functions as funct
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as pl
from mpl_toolkits.basemap import Basemap

for year in ['2010', '2011', '2012', '2013', '2014', '2015', '2016']:

    for month in ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']:

        if os.path.isdir('/Users/jmh2g09/Documents/PhD/Data/elev_files/' + year + month + '_elev'):
            os.chdir('/Users/jmh2g09/Documents/PhD/Data/elev_files/' + year + month + '_elev')
            directory = os.getcwd()

            print(month, year)

            data = funct.month_data(directory, lead_offset=funct.ocean_lead_offset(month))
    
            lat = data['lat']
            lon = data['lon']
            ssh = data['ssh']
            ice_conc = data['ice_conc']
            surface_type = data['type']
            print('Data Extracted')
            
            # Make a plot of the data for each month
            pl.figure()
            pl.clf()
            m = Basemap(projection='spstere', boundinglat=-50, lon_0=180, resolution='l')
            m.drawmapboundary()
            m.drawcoastlines(zorder=10)
            m.fillcontinents(zorder=10)
            m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
            m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
        
            stereo_x, stereo_y = m(lon, lat)
            m.scatter(stereo_x, stereo_y, c=ssh, marker='.', edgecolors='none')
            m.colorbar()
            pl.clim(np.mean(ssh) + 3*np.std(ssh), np.mean(ssh) - 3*np.std(ssh))
            pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Raw/'+ year + '/Figures/' + year + month + '.png', 
                format='png')
            pl.close()
            
            # Put the data in a .nc file in /Users/jmh2g09/Documents/PhD/Data/Raw
            os.chdir('/Users/jmh2g09/Documents/PhD/Data/Raw/' + year)

            nc = Dataset(year + month + '_raw.nc', 'w', format='NETCDF3_CLASSIC')

            nc.createDimension('station', np.size(lat))

            latitudes = nc.createVariable('lat', float, ('station',))
            longitudes = nc.createVariable('lon', float, ('station',))
            sea_surface_height = nc.createVariable('sea_surface_height', float, ('station',))
            surface = nc.createVariable('surface_type', 'i', ('station',))
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
    
            crs.grid_mapping_name = 'latitude_longitude'
            crs.semi_major_axis = 6378137.
            crs.inverse_flattening = 298.257222101004
            crs.earth_gravity_constant = 398600500000000.
            crs.earth_rotation_rate = 7.292115e-05

            latitudes[:] = lat
            longitudes[:] = lon
            sea_surface_height[:] = ssh
            surface[:] = surface_type
    
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
            print('Generating geoid file')
            os.system('gut geoidheight_tf -InFile GOCO05s.gfc -T tide-free -InTrack '
            + year + month + '_track.nc -OutFile ' + year + month + '_geoid_track.nc')