import numpy as np
import os
from netCDF4 import Dataset
import matplotlib.pyplot as pl
from mpl_toolkits.basemap import Basemap

nc = Dataset('/Users/jmh2g09/Documents/PhD/Data/Gridded/DOT/MDT_mean_ice_only.nc', 'r')
mean_year = nc.variables['mean_dynamic_topography_seasonal_offset'][:]
mean_year_2 = nc.variables['mean_dynamic_topography_no_offset'][:]
mean_year_3 = nc.variables['mean_dynamic_topography_constant_offset'][:]
lat = nc.variables['latitude'][:]
lon = nc.variables['longitude'][:]
nc.close()

anom_seasonal = np.full((12, len(lat), len(lon), 7), fill_value=np.NaN)
anom_seasonal_2 = np.full((12, len(lat), len(lon), 7), fill_value=np.NaN)
anom_seasonal_3 = np.full((12, len(lat), len(lon), 7), fill_value=np.NaN)
ice_seasonal = np.full((12, len(lat), len(lon), 7), fill_value=np.NaN)
            
for year in ['2010', '2011', '2012', '2013', '2014', '2015', '2016']:
    os.chdir('/Users/jmh2g09/Documents/PhD/Data/Gridded/DOT/' + year)
    for file in os.listdir():
        if file[-6:] == 'ice.nc':
            print(file)
            month = file[4:6]
            nc = Dataset(file, 'r')
            dot = nc.variables['dynamic_ocean_topography_seasonal_offset'][:]
            dot_2 = nc.variables['dynamic_ocean_topography_no_offset'][:]
            dot_3 = nc.variables['dynamic_ocean_topography_constant_offset'][:]
            ice_data = nc.variables['sea_ice_concentration'][:]
            nc.close()
            
            nc = Dataset('../' + str(int(month)) + '_mean_DOT_filt_ice_only.nc', 'r')
            dot_month_mean = nc.variables['dynamic_ocean_topography_seasonal_offset'][:]
            dot_month_mean_2 = nc.variables['dynamic_ocean_topography_no_offset'][:]
            dot_month_mean_3 = nc.variables['dynamic_ocean_topography_constant_offset'][:]
            nc.close()
            
            dot_anomaly = dot - mean_year
            dot_anomaly_2 = dot_2 - mean_year_2
            dot_anomaly_3 = dot_3 - mean_year_3
            
            dot_month_anomaly = dot - dot_month_mean
            dot_month_anomaly_2 = dot_2 - dot_month_mean_2
            dot_month_anomaly_3 = dot_3 - dot_month_mean_3
            
            dot_anomaly_month = dot - mean_year
            dot_anomaly_month[np.isnan(dot_anomaly_month)] = 999
            dot_anomaly_month_2 = dot_2 - mean_year_2
            dot_anomaly_month_2[np.isnan(dot_anomaly_month_2)] = 999
            dot_anomaly_month_3 = dot_3 - mean_year_3
            dot_anomaly_month_3[np.isnan(dot_anomaly_month_3)] = 999
            
            anom_seasonal[int(month)-1, :, :, int(year)-2010] = dot_anomaly_month
            anom_seasonal_2[int(month)-1, :, :, int(year)-2010] = dot_anomaly_month_2
            anom_seasonal_3[int(month)-1, :, :, int(year)-2010] = dot_anomaly_month_3
            ice_seasonal[int(month)-1, :, :, int(year)-2010] = ice_data
            
            nc = Dataset('Anomalies/' + year + month + '_DOT_anomaly_ice.nc', 'w')
        
            nc.createDimension('lat', np.size(lat))
            nc.createDimension('lon', np.size(lon))
        
            latitudes = nc.createVariable('latitude', float, ('lat',))
            longitudes = nc.createVariable('longitude', float, ('lon',))
            dot_anom = nc.createVariable('dynamic_ocean_topography_anomaly_seasonal_offset', float, ('lat', 'lon'))
            dot_anom_2 = nc.createVariable('dynamic_ocean_topography_anomaly_no_offset', float, ('lat', 'lon'))
            dot_anom_3 = nc.createVariable('dynamic_ocean_topography_anomaly_constant_offset', float, ('lat', 'lon'))
            ice = nc.createVariable('sea_ice_concentration', float, ('lat', 'lon'))
            latitudes[:] = lat
            longitudes[:] = lon
            dot_anom[:] = dot_anomaly
            dot_anom_2[:] = dot_anomaly_2
            dot_anom_3[:] = dot_anomaly_3
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
            m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(dot_anomaly_2)), cmap='RdBu_r')
            m.colorbar()
            pl.clim(-.3, .3)
            pl.savefig('Anomalies/Figures/' + str(year) + '_' + str(month) + '_DOT_anomaly_ice.png', format='png', transparent=True, dpi=300)
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
            m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(dot_month_anomaly)), cmap='RdBu_r')
            m.colorbar()
            pl.clim(-.3, .3)
            pl.savefig('Anomalies/Figures/' + str(year) + '_' + str(month) + '_DOT_month_anomaly_ice.png', format='png', transparent=True, dpi=300)
            pl.close()

seasonal = np.nanmean(anom_seasonal, axis=3)
seasonal_2 = np.nanmean(anom_seasonal_2, axis=3)
seasonal_3 = np.nanmean(anom_seasonal_3, axis=3)
seasonal[seasonal > 200] = np.NaN
seasonal_2[seasonal_2 > 200] = np.NaN
seasonal_3[seasonal_3 > 200] = np.NaN
ice_season = np.nanmean(ice_seasonal, axis=3)

for imnth in range(12):
    pl.figure()
    pl.clf()
    m = Basemap(projection='spstere', boundinglat=-50, lon_0=180, resolution='l')
    m.drawmapboundary()
    m.drawcoastlines(zorder=10)
    m.fillcontinents(zorder=10)
    m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
    stereo_x, stereo_y = m(grid_lons, grid_lats)
    m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(seasonal[imnth, :, :])), cmap='RdBu_r')
    m.colorbar()
    pl.clim(-.1, .1)
    pl.savefig('../Figures/' + str(imnth + 1) + '_DOT_anomaly_ice_only.png', format='png', transparent=True, dpi=300)
    pl.close()