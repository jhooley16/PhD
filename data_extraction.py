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
            surface = data['surface']
            mode = data['mode']
            time = data['time']
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

            if year == '2015':
                pl.figure()
                pl.clf()
                m = Basemap(projection='spstere', boundinglat=-50, lon_0=180, resolution='l')
                m.drawmapboundary()
                m.drawcoastlines(zorder=10)
                m.fillcontinents(zorder=10)
                m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
                m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
                stereo_x, stereo_y = m(lon, lat)
                m.scatter(stereo_x, stereo_y, c=mode, marker='.', edgecolors='none', s=20)
                m.colorbar()
                pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Processed/SeparateModes/Figures/' + month + '_modes.png', dpi=300, transparent=True)
                pl.close()

            # Put the data in a .nc file in /Users/jmh2g09/Documents/PhD/Data/Processed
            os.chdir('/Users/jmh2g09/Documents/PhD/Data/Processed/' + year)
            nc = Dataset(year + month + '_raw.nc', 'w')

            nc.createDimension('station', np.size(lat))

            latitudes = nc.createVariable('latitude', float, ('station',))
            longitudes = nc.createVariable('longitude', float, ('station',))
            sea_surface_height = nc.createVariable('sea_surface_height_seasonal_offset', float, ('station',))
            surface_type = nc.createVariable('surface_type', 'i', ('station',))
            sea_ice_concentration = nc.createVariable('ice_concentration', float, ('station',))
            time_var = nc.createVariable('time', float, ('station',))
            mode_var = nc.createVariable('mode', 'i', ('station',))

            latitudes.long_name = 'latitude'
            latitudes.standard_name = 'latitude'
            latitudes.units = 'degrees_north'
            longitudes.long_name = 'longitude'
            longitudes.standard_name = 'longitude'
            longitudes.units = 'degrees_east'
            surface_type.standard_name = 'surface_type'
            surface_type.long_name = 'ocean_(1)_or_lead_(2)_surface'
            sea_ice_concentration.standard_name = 'sea_ice_concentration'
            sea_ice_concentration.long_name = 'sea_ice_concentration'
            sea_ice_concentration.units = '%'
            sea_surface_height.standard_name = 'sea_surface_height_above_WGS84_ellipsoid_seasonal_offset_correction'
            sea_surface_height.long_name = 'sea_surface_height'
            sea_surface_height.units = 'm'
            time_var.standard_name = 'time'
            time_var.long_name = 'sea_surface_height'
            time_var.units = 'days_since_1/1/1950'
            mode_var.standard_name = 'satellite_mode_type'
            mode_var.long_name = 'LRM_(0)_or_SAR_(1)_or_SARIn_(2)'

            latitudes[:] = lat
            longitudes[:] = lon
            sea_surface_height[:] = ssh
            surface_type[:] = surface
            sea_ice_concentration[:] = ice_conc
            time_var[:] = time
            mode_var[:] = mode
    
            nc.close()