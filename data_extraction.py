import os
import functions as funct
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as pl
from mpl_toolkits.basemap import Basemap
from datetime import date

for year in ['2012']:#, '2012', '2013', '2014', '2015', '2016']:
    for month in ['08', '09']:#, '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']:
        directory = '/Volumes/My Passport/Data/elev_files/' + year + month + '_MERGE'
        if os.path.isdir(directory):
            print(month, year)
            
            data = funct.month_data(directory, month)
            lat = data['lat']
            lon = data['lon']
            ssh = data['ssh']
            ice_conc = data['ice_conc']
            surface = data['surface']
            mode = data['mode']
            time = np.array(data['time']) + date.toordinal(date(1950, 1, 1))
            print('Data Extracted')

            day = []
            for it in range(len(time)):
                d = date.fromordinal(int(time[it]))
                t = d.timetuple()
                day.append(t[2])

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
            m.scatter(stereo_x, stereo_y, c=ssh, marker='.', edgecolors='none', s=2)
            c = m.colorbar()
            c.set_label('Sea Surface Height (m)')
            pl.clim(70, -70)
            pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Processed/'+ year + '/Figures/' + year + month + '_ssh.png', transparent=True, dpi=300, bbox_inches='tight')
            pl.close()
            
            # Make a plot of the time each point was taken
            pl.figure()
            pl.clf()
            m = Basemap(projection='spstere', boundinglat=-50, lon_0=180, resolution='l')
            m.drawmapboundary()
            m.drawcoastlines(zorder=10)
            m.fillcontinents(zorder=10)
            m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
            m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
            stereo_x, stereo_y = m(lon, lat)
            m.scatter(stereo_x, stereo_y, c=day, marker='.', edgecolors='none', s=2)
            c = m.colorbar()
            c.set_label('Day of Month')
            pl.clim(30, 0)
            pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Processed/'+ year + '/Figures/' + year + month + '_days.png', transparent=True, dpi=300, bbox_inches='tight')
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
                m.scatter(stereo_x, stereo_y, c=mode, marker='.', edgecolors='none', s=2)
                pl.savefig('/Users/jmh2g09/Documents/PhD/Data/SeparateModes/Figures/' + month + '_modes.png', transparent=True, dpi=300, bbox_inches='tight')
                pl.close()
                
                pl.figure()
                pl.clf()
                m = Basemap(projection='spstere', boundinglat=-50, lon_0=180, resolution='l')
                m.drawmapboundary()
                m.drawcoastlines(zorder=10)
                m.fillcontinents(zorder=10)
                m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
                m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
                stereo_x, stereo_y = m(lon, lat)
                m.scatter(stereo_x, stereo_y, c=surface, marker='.', edgecolors='none', s=2)
                pl.savefig('/Users/jmh2g09/Documents/PhD/Data/SeparateModes/Figures/' + month + '_surfaces.png', transparent=True, dpi=300, bbox_inches='tight')
                pl.close()
        
            # Put the data in a .nc file 
            nc = Dataset('/Users/jmh2g09/Documents/PhD/Data/Processed/' + year + month + '_track.nc', 'w')
            nc.createDimension('station', np.size(lat))

            latitudes = nc.createVariable('latitude', float, ('station',))
            longitudes = nc.createVariable('longitude', float, ('station',))
            sea_surface_height = nc.createVariable('sea_surface_height', float, ('station',))
            surface_type = nc.createVariable('surface', 'i', ('station',))
            sea_ice_concentration = nc.createVariable('ice_concentration', float, ('station',))
            time_var = nc.createVariable('time', float, ('station',))
            mode_var = nc.createVariable('mode', 'i', ('station',))

            latitudes.standard_name = 'latitude'
            latitudes.units = 'degrees_north'
            longitudes.standard_name = 'longitude'
            longitudes.units = 'degrees_east'
            surface_type.standard_name = 'surface_type'
            surface_type.units = 'ocean_(1)_or_lead_(2)_surface'
            sea_ice_concentration.standard_name = 'sea_ice_concentration'
            sea_ice_concentration.units = '%'
            sea_surface_height.standard_name = 'sea_surface_height_above_WGS84'
            sea_surface_height.units = 'm'
            time_var.standard_name = 'time'
            time_var.units = 'days_since_1/1/1950'
            mode_var.standard_name = 'satellite_mode_type'
            mode_var.units = 'LRM_(0)_or_SAR_(1)_or_SARIn_(2)'

            latitudes[:] = lat
            longitudes[:] = lon
            sea_surface_height[:] = ssh
            sea_ice_concentration[:] = ice_conc
            surface_type[:] = surface
            time_var[:] = time
            mode_var[:] = mode
    
            nc.close()