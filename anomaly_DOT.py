import numpy as np
import os
from netCDF4 import Dataset
import matplotlib.pyplot as pl
from mpl_toolkits.basemap import Basemap

# Load the total mean DOT data 
nc = Dataset('/Users/jmh2g09/Documents/PhD/Data/Gridded/Means/mean.nc', 'r')
mean_dot = nc.variables['mean_dynamic_topography_seasonal_offset'][:]
mean_2_dot = nc.variables['mean_dynamic_topography_constant_offset'][:]
mean_ssh = nc.variables['mean_sea_surface_height_seasonal_offset'][:]
mean_2_ssh = nc.variables['mean_sea_surface_height_constant_offset'][:]
monthly_mean_dot = nc.variables['monthly_mean_dynamic_topography_constant_offset'][:]
monthly_mean_2_dot = nc.variables['monthly_mean_dynamic_topography_constant_offset'][:]
monthly_mean_ssh = nc.variables['monthly_mean_sea_surface_height_constant_offset'][:]
monthly_mean_2_ssh = nc.variables['monthly_mean_sea_surface_height_constant_offset'][:]
lat = nc.variables['lat'][:]
lon = nc.variables['lon'][:]
nc.close()

for year in ['2011', '2012', '2013', '2014', '2015', '2016']:
    for month in ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']:
        file = '/Users/jmh2g09/Documents/PhD/Data/Gridded/' + year + month + '_grid.nc'
        print(file)

        # Load the DOT file
        nc = Dataset(file, 'r')
        dot = nc.variables['filtered_dynamic_ocean_topography_seasonal_offset'][:]
        dot_2 = nc.variables['filtered_dynamic_ocean_topography_constant_offset'][:]
        ssh = nc.variables['filtered_sea_surface_height_seasonal_offset'][:]
        ssh_2 = nc.variables['filtered_sea_surface_height_constant_offset'][:]
        ice_data = nc.variables['filtered_sea_ice_concentration'][:]
        nc.close()
        
        # Calculate the month anomaly from the total mean
        dot_anomaly = dot - mean_dot
        dot_2_anomaly = dot_2 - mean_2_dot
        ssh_anomaly = ssh - mean_ssh
        ssh_2_anomaly = ssh_2 - mean_2_ssh
        
        # Calculate the month anomaly from the month mean
        dot_anomaly_month = dot - monthly_mean_dot[:, :, int(month) - 1]
        dot_2_anomaly_month = dot_2 - monthly_mean_2_dot[:, :, int(month) - 1]
        ssh_anomaly_month = ssh - monthly_mean_ssh[:, :, int(month) - 1]
        ssh_2_anomaly_month = ssh_2 - monthly_mean_2_ssh[:, :, int(month) - 1]
        
        # Save both anomalies
        # remove already-existing variables to overwrite
        os.system('ncks -O -x -v dynamic_ocean_topography_anomaly_seasonal_offset,\
dynamic_ocean_topography_anomaly_constant_offset,\
monthly_dynamic_ocean_topography_anomaly_seasonal_offset,\
monthly_dynamic_ocean_topography_anomaly_constant_offset,\
sea_surface_height_anomaly_seasonal_offset,\
sea_surface_height_anomaly_constant_offset,\
monthly_sea_surface_height_anomaly_seasonal_offset,\
monthly_sea_surface_height_anomaly_constant_offset ' + file + ' ' + file)
        # Save the file
        nc = Dataset(file, 'a')

        dot_anom = nc.createVariable('dynamic_ocean_topography_anomaly_seasonal_offset', float, ('lat', 'lon'))
        dot_2_anom = nc.createVariable('dynamic_ocean_topography_anomaly_constant_offset', float, ('lat', 'lon'))
        dot_anom_month = nc.createVariable('monthly_dynamic_ocean_topography_anomaly_seasonal_offset', float, ('lat', 'lon'))
        dot_2_anom_month = nc.createVariable('monthly_dynamic_ocean_topography_anomaly_constant_offset', float, ('lat', 'lon'))
        ssh_anom = nc.createVariable('sea_surface_height_anomaly_seasonal_offset', float, ('lat', 'lon'))
        ssh_2_anom = nc.createVariable('sea_surface_height_anomaly_constant_offset', float, ('lat', 'lon'))
        ssh_anom_month = nc.createVariable('monthly_sea_surface_height_anomaly_seasonal_offset', float, ('lat', 'lon'))
        ssh_2_anom_month = nc.createVariable('monthly_sea_surface_height_anomaly_constant_offset', float, ('lat', 'lon'))
    
        dot_anom.standard_name = 'sea_surface_height_above_EIGEN6c4_anomaly_seasonal_offset'
        dot_anom.units = 'm'
        dot_2_anom.standard_name = 'sea_surface_height_above_EIGEN6c4_anomaly_constant_offset'
        dot_2_anom.units = 'm'
        dot_anom_month.standard_name = 'monthly_sea_surface_height_above_EIGEN6c4_anomaly_seasonal_offset'
        dot_anom_month.units = 'm'
        dot_2_anom_month.standard_name = 'monthly_sea_surface_height_above_EIGEN6c4_anomaly_constant_offset'
        dot_2_anom_month.units = 'm'
        ssh_anom.standard_name = 'sea_surface_height_above_WGS84_anomaly_seasonal_offset'
        ssh_anom.units = 'm'
        ssh_2_anom.standard_name = 'sea_surface_height_above_WGS84_anomaly_constant_offset'
        ssh_2_anom.units = 'm'
        ssh_anom_month.standard_name = 'monthly_sea_surface_height_above_WGS84_anomaly_seasonal_offset'
        ssh_anom_month.units = 'm'
        ssh_2_anom_month.standard_name = 'monthly_sea_surface_height_above_WGS84_anomaly_constant_offset'
        ssh_2_anom_month.units = 'm'
        
        dot_anom[:] = dot_anomaly
        dot_2_anom[:] = dot_2_anomaly
        dot_anom_month[:] = dot_anomaly_month
        dot_2_anom_month[:] = dot_2_anomaly_month
        ssh_anom[:] = ssh_anomaly
        ssh_2_anom[:] = ssh_2_anomaly
        ssh_anom_month[:] = ssh_anomaly_month
        ssh_2_anom_month[:] = ssh_2_anomaly_month

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
        pl.clim(-.1,.1)
        m.contour(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(ice_data)), colors='k', levels=[20])
        pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Gridded/Figures/' + year + '/' + year + month + '_DOT_anomaly.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
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
        pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Gridded/Figures/' + year + '/' + year + month + '_DOT_month_anomaly.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
        pl.close()