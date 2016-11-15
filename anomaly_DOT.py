import numpy as np
import os
from netCDF4 import Dataset
import matplotlib.pyplot as pl
from mpl_toolkits.basemap import Basemap

for year in ['2010', '2011', '2012', '2013', '2014', '2015', '2016']:

    os.chdir('/Users/jmh2g09/Documents/PhD/Data/Gridded/DOT/' + year)

    nc = Dataset('/Users/jmh2g09/Documents/PhD/Data/Gridded/DOT/MDT_mean.nc', 'r')
    mean_year = nc.variables['mean_dynamic_topography_seasonal_offset'][:]
    mean_2_year = nc.variables['mean_dynamic_topography_no_offset'][:]
    mean_3_year = nc.variables['mean_dynamic_topography_constant_offset'][:]
    lat = nc.variables['latitude'][:]
    lon = nc.variables['longitude'][:]
    nc.close()
    for file in os.listdir():
        if file[-11:] == 'DOT_filt.nc':
            print(file)
            month = file[4:6]
            nc = Dataset(file, 'r')
        
            dot = nc.variables['dynamic_ocean_topography_seasonal_offset'][:]
            dot_2 = nc.variables['dynamic_ocean_topography_no_offset'][:]
            dot_3 = nc.variables['dynamic_ocean_topography_constant_offset'][:]
            ice_data = nc.variables['sea_ice_concentration'][:]
        
            nc.close()
            dot_anomaly = np.subtract(dot, mean_year)
            dot_2_anomaly = np.subtract(dot_2, mean_2_year)
            dot_3_anomaly = np.subtract(dot_3, mean_3_year)
            
            nc = Dataset('Anomalies/' + year + month + '_DOT_anomaly.nc', 'w')
        
            nc.createDimension('lat', np.size(lat))
            nc.createDimension('lon', np.size(lon))
        
            latitudes = nc.createVariable('latitude', float, ('lat',))
            longitudes = nc.createVariable('longitude', float, ('lon',))
            dot_anom = nc.createVariable('dynamic_ocean_topography_anomaly_seasonal_offset', float, ('lat', 'lon'))
            dot_2_anom = nc.createVariable('dynamic_ocean_topography_anomaly_no_offset', float, ('lat', 'lon'))
            dot_3_anom = nc.createVariable('dynamic_ocean_topography_anomaly_constant_offset', float, ('lat', 'lon'))
            ice = nc.createVariable('sea_ice_concentration', float, ('lat', 'lon'))
            latitudes[:] = lat
            longitudes[:] = lon
            dot_anom[:] = dot_anomaly
            dot_2_anom[:] = dot_2_anomaly
            dot_3_anom[:] = dot_3_anomaly
            ice[:] = ice_data

            nc.close()

            grid_lats, grid_lons = np.meshgrid(lat, lon)
        
            pl.figure()
            pl.clf()
            m = Basemap(projection='spstere', boundinglat=-50, lon_0=180, resolution='l')
            m.drawmapboundary()
            m.drawcoastlines(zorder=10)
            m.fillcontinents(zorder=10)
            m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
            m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
            stereo_x, stereo_y = m(grid_lons, grid_lats)
            m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(dot_anomaly)), cmap='RdBu_r')
            m.colorbar()
            pl.clim(.3, -.3)
            m.contour(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(ice_data)), colors='k', levels=[20])
            pl.savefig('Anomalies/Figures/' + str(year) + '_' + str(month) + '_DOT_anomaly.png', format='png', transparent=True, dpi=300)
            pl.close()
            
            pl.figure()
            pl.clf()
            m = Basemap(projection='spstere', boundinglat=-50, lon_0=180, resolution='l')
            m.drawmapboundary()
            m.drawcoastlines(zorder=10)
            m.fillcontinents(zorder=10)
            m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
            m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
            stereo_x, stereo_y = m(grid_lons, grid_lats)
            m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(dot_2_anomaly)), cmap='RdBu_r')
            m.colorbar()
            pl.clim(.3, -.3)
            m.contour(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(ice_data)), colors='k', levels=[20])
            pl.savefig('Anomalies/Figures/' + str(year) + '_' + str(month) + '_DOT_anomaly_no_offset.png', format='png', transparent=True, dpi=300)
            pl.close()
            
            pl.figure()
            pl.clf()
            m = Basemap(projection='spstere', boundinglat=-50, lon_0=180, resolution='l')
            m.drawmapboundary()
            m.drawcoastlines(zorder=10)
            m.fillcontinents(zorder=10)
            m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
            m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
            stereo_x, stereo_y = m(grid_lons, grid_lats)
            m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(dot_3_anomaly)), cmap='RdBu_r')
            m.colorbar()
            pl.clim(.3, -.3)
            m.contour(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(ice_data)), colors='k', levels=[20])
            pl.savefig('Anomalies/Figures/' + str(year) + '_' + str(month) + '_DOT_anomaly_constant_offset.png', format='png', transparent=True, dpi=300)
            pl.close()