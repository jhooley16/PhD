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

            data = funct.month_data(directory, month)
            lat = data['lat']
            lon = data['lon']
            ssh = data['ssh']
            ssh_2 = data['ssh_2']
            ssh_3 = data['ssh_3']
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
            pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Processed/'+ year + '/Figures/' + year + month + '_ssh.png', 
                format='png', dpi=300)
            pl.close()
            
            # Put the data in a .nc file in /Users/jmh2g09/Documents/PhD/Data/Processed
            os.chdir('/Users/jmh2g09/Documents/PhD/Data/Processed/' + year)
            nc = Dataset(year + month + '_raw.nc', 'w')

            nc.createDimension('station', np.size(lat))

            latitudes = nc.createVariable('latitude', float, ('station',))
            longitudes = nc.createVariable('longitude', float, ('station',))
            sea_surface_height = nc.createVariable('sea_surface_height_seasonal_offset', float, ('station',))
            sea_surface_height_2 = nc.createVariable('sea_surface_height_no_offset', float, ('station',))
            sea_surface_height_3 = nc.createVariable('sea_surface_height_constant_offset', float, ('station',))
            surface = nc.createVariable('surface_type', 'i', ('station',))
            sea_ice_concentration = nc.createVariable('ice_concentration', float, ('station',))
    
            latitudes.long_name = 'latitude'
            latitudes.standard_name = 'latitude'
            latitudes.units = 'degrees_north'
            longitudes.long_name = 'longitude'
            longitudes.standard_name = 'longitude'
            longitudes.units = 'degrees_east'
            surface.standard_name = 'surface_type'
            surface.long_name = 'ocean_(1)_or_lead_(2)_surface'
            surface.grid_mapping = 'crs'
            sea_ice_concentration.standard_name = 'sea_ice_concentration'
            sea_ice_concentration.long_name = 'sea_ice_concentration'
            sea_ice_concentration.units = '%'
            sea_ice_concentration.grid_mapping = 'crs'
            sea_surface_height.standard_name = 'sea_surface_height_above_WGS84_ellipsoid_seasonal_offset_correction'
            sea_surface_height.long_name = 'sea_surface_height'
            sea_surface_height.units = 'm'
            sea_surface_height.grid_mapping = 'crs'
            sea_surface_height.tide_system = 'tide_free'
            sea_surface_height_2.standard_name = 'sea_surface_height_above_WGS84_ellipsoid_no_offset_correction'
            sea_surface_height_2.long_name = 'sea_surface_height'
            sea_surface_height_2.units = 'm'
            sea_surface_height_2.grid_mapping = 'crs'
            sea_surface_height_2.tide_system = 'tide_free'
            sea_surface_height_3.standard_name = 'sea_surface_height_above_WGS84_ellipsoid_constant_offset_correction'
            sea_surface_height_3.long_name = 'sea_surface_height'
            sea_surface_height_3.units = 'm'
            sea_surface_height_3.grid_mapping = 'crs'
            sea_surface_height_3.tide_system = 'tide_free'

            latitudes[:] = lat
            longitudes[:] = lon
            sea_surface_height[:] = ssh
            sea_surface_height_2[:] = ssh_2
            sea_surface_height_3[:] = ssh_3
            surface[:] = surface_type
            sea_ice_concentration[:] = ice_conc
    
            nc.close()