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
            sea_surface_height = nc.createVariable('sea_surface_height', float, ('station',))
            surface = nc.createVariable('surface_type', 'i', ('station',))
            sea_ice_concentration = nc.createVariable('ice_concentration', float, ('station',))
            crs = nc.createVariable('crs', 'i', ())
    
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
            sea_surface_height.standard_name = 'sea_surface_height_above_WGS84_ellipsoid'
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
            sea_ice_concentration[:] = ice_conc
    
            nc.close()
            
            # Generate geoid points from EGM2008 script
            os.chdir('/Users/jmh2g09/Documents/PhD/Geoid/EGM2008/')
            # Create the track file for the script (2 column, lat lon file)
            with open('INPUT.DAT', 'w') as f:
                for f1, f2 in zip(lat, lon):
                    print(f1, f2, file=f)
            # Delete the OUTPUT.DAT file, so that the script works
            os.system('rm OUTPUT.DAT')
            # Run the Script, producing the OUTPUT.DAT for this month
            os.system('./EGM2008_interp.out')
            # Read output.dat and create variables
            f = open('OUTPUT.DAT', 'r')
            lat_geoid = []
            lon_geoid = []
            geoid = []
            for line in f:
                line = line.strip()
                column = line.split()
                lat_geoid.append(float(column[0]))
                lon_geoid.append(float(column[1]))
                geoid.append(float(column[2]))
            f.close()

            # Make a plot of the data for each month
            pl.figure()
            pl.clf()
            m = Basemap(projection='spstere', boundinglat=-50, lon_0=180, resolution='l')
            m.drawmapboundary()
            m.drawcoastlines(zorder=10)
            m.fillcontinents(zorder=10)
            m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
            m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
            stereo_x, stereo_y = m(lon_geoid, lat_geoid)
            m.scatter(stereo_x, stereo_y, c=geoid, marker='.', edgecolors='none')
            m.colorbar()
            pl.clim(np.mean(ssh) + 3*np.std(ssh), np.mean(ssh) - 3*np.std(ssh))
            pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Processed/'+ year + '/Figures/' + year + month + '_geoid.png', 
                format='png', dpi=300)
            pl.close()

            # Write the data into a .nc file
            nc = Dataset('/Users/jmh2g09/Documents/PhD/Data/Processed/' + year + '/' + year + month + '_geoid_track.nc', 'w')
            nc.createDimension('station', np.size(lat_geoid))

            latitudes = nc.createVariable('latitude', float, ('station',))
            longitudes = nc.createVariable('longitude', float, ('station',))
            geoid_height = nc.createVariable('geoid_height', float, ('station',))
            crs = nc.createVariable('crs', 'i', ())

            latitudes.long_name = 'latitude'
            latitudes.standard_name = 'latitude'
            latitudes.units = 'degrees_north'
            longitudes.long_name = 'longitude'
            longitudes.standard_name = 'longitude'
            longitudes.units = 'degrees_east'
            geoid_height.standard_name = 'geoid_height'
            geoid_height.long_name = 'EGM08_geoid_height_above_WGS84_ellipsoid'
            geoid_height.units = 'm'
            geoid_height.grid_mapping = 'crs'

            crs.grid_mapping_name = 'latitude_longitude'
            crs.semi_major_axis = 6378137.
            crs.inverse_flattening = 298.257222101004
            crs.earth_gravity_constant = 398600500000000.
            crs.earth_rotation_rate = 7.292115e-05

            latitudes[:] = lat_geoid
            longitudes[:] = lon_geoid
            geoid_height[:] = geoid

            nc.close()