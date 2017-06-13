import numpy as np
import os
from netCDF4 import Dataset
import matplotlib.pyplot as pl
from mpl_toolkits.basemap import Basemap

# Load the total mean DOT data 
nc = Dataset('/Users/jmh2g09/Documents/PhD/Data/Gridded/DOT/Means/MDT_mean.nc', 'r')
mean_year = nc.variables['mean_dynamic_topography_seasonal_offset'][:]
mean_2_year = nc.variables['mean_dynamic_topography_constant_offset'][:]
lat = nc.variables['latitude'][:]
lon = nc.variables['longitude'][:]
nc.close()

for year in ['2010', '2011', '2012', '2013', '2014', '2015', '2016']:
    os.chdir('/Users/jmh2g09/Documents/PhD/Data/Gridded/DOT/' + year)
    for month in ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']:
        file = year + month + '_DOT_filt.nc'
        if os.path.isfile(file):
            print(file)

            # Load the DOT file
            nc = Dataset(file, 'r')
            dot = nc.variables['dynamic_ocean_topography_seasonal_offset'][:]
            dot_2 = nc.variables['dynamic_ocean_topography_constant_offset'][:]
            ice_data = nc.variables['sea_ice_concentration'][:]
            nc.close()

            # Load the mean for that month
            nc = Dataset('/Users/jmh2g09/Documents/PhD/Data/Gridded/DOT/Means/MDT_mean_' + month + '.nc', 'r')
            monthly_mean = nc.variables['mean_dynamic_topography_seasonal_offset'][:]
            monthly_mean_2 = nc.variables['mean_dynamic_topography_constant_offset'][:]
            nc.close()
            
            # Calculate the month anomaly from the total mean
            dot_anomaly = dot - mean_year
            dot_2_anomaly = dot_2 - mean_2_year
            
            # Calculate the month anomaly from the month mean
            dot_anomaly_month = dot - monthly_mean
            dot_2_anomaly_month = dot_2 - monthly_mean_2
            
            # Save both anomalies
            nc = Dataset('Anomalies/' + year + month + '_DOT_anomaly.nc', 'w')
            nc.createDimension('lat', np.size(lat))
            nc.createDimension('lon', np.size(lon))
            latitudes = nc.createVariable('latitude', float, ('lat',))
            longitudes = nc.createVariable('longitude', float, ('lon',))
            dot_anom = nc.createVariable('dynamic_ocean_topography_anomaly_seasonal_offset', float, ('lat', 'lon'))
            dot_2_anom = nc.createVariable('dynamic_ocean_topography_anomaly_constant_offset', float, ('lat', 'lon'))
            dot_anom_month = nc.createVariable('monthly_dynamic_ocean_topography_anomaly_seasonal_offset', float, ('lat', 'lon'))
            dot_2_anom_month = nc.createVariable('monthly_dynamic_ocean_topography_anomaly_constant_offset', float, ('lat', 'lon'))
            ice = nc.createVariable('sea_ice_concentration', float, ('lat', 'lon'))
            latitudes[:] = lat
            longitudes[:] = lon
            dot_anom[:] = dot_anomaly
            dot_2_anom[:] = dot_2_anomaly
            dot_anom_month[:] = dot_anomaly_month
            dot_2_anom_month[:] = dot_2_anomaly_month
            ice[:] = ice_data
            nc.close()
            
            # Plot the anomaly from the total mean
            pl.figure()
            pl.clf()
            m = Basemap(projection='spstere', boundinglat=-50, lon_0=180, resolution='l')
            m.drawmapboundary()
            m.drawcoastlines(zorder=10)
            m.fillcontinents(zorder=10)
            m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
            m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
            grid_lats, grid_lons = np.meshgrid(lat, lon)
            stereo_x, stereo_y = m(grid_lons, grid_lats)
            m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(dot_anomaly)), cmap='RdBu_r')
            c = m.colorbar()
            c.set_label('DOT anomaly (m)')
            pl.clim(-.1, .1)
            m.contour(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(ice_data)), colors='k', levels=[20])
            pl.savefig('Anomalies/Figures/' + str(year) + str(month) + '_DOT_anomaly.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
            pl.close()
            
            # Plot the anomaly from the month mean
            pl.figure()
            pl.clf()
            m = Basemap(projection='spstere', boundinglat=-50, lon_0=180, resolution='l')
            m.drawmapboundary()
            m.drawcoastlines(zorder=10)
            m.fillcontinents(zorder=10)
            m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
            m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
            grid_lats, grid_lons = np.meshgrid(lat, lon)
            stereo_x, stereo_y = m(grid_lons, grid_lats)
            m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(dot_anomaly_month)), cmap='RdBu_r')
            c = m.colorbar()
            c.set_label('DOT anomaly (m)')
            pl.clim(-.1, .1)
            m.contour(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(ice_data)), colors='k', levels=[20])
            pl.savefig('Anomalies/Figures/' + str(year) + str(month) + '_month_DOT_anomaly.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
            pl.close()