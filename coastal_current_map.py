import os
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as pl
from mpl_toolkits.basemap import Basemap
from geopy.distance import vincenty

dy = 0.5 * 60 * 1.862

ddot_dy_seasonal = np.full((59, 361, 7, 12), fill_value=np.NaN)
ddot_dy_constant = np.full((59, 361, 7, 12), fill_value=np.NaN)
for year in ['2010', '2011', '2012', '2013', '2014', '2015', '2016']:
    os.chdir('/Users/jmh2g09/Documents/PhD/Data/Gridded/DOT/' + year)
    for month in ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']:
        if os.path.isfile(year + month + '_DOT_filt.nc'):
            nc = Dataset(year + month + '_DOT_filt.nc', 'r')
            dot_seasonal = nc.variables['dynamic_ocean_topography_seasonal_offset'][:]
            dot_constant = nc.variables['dynamic_ocean_topography_constant_offset'][:]
            lat = nc.variables['latitude'][:]
            lon = nc.variables['longitude'][:]
            nc.close()

            ddot_dy_seasonal[:, :, int(year)-2010, int(month)-1] = np.gradient(dot_seasonal, axis=0) / dy
            ddot_dy_constant[:, :, int(year)-2010, int(month)-1] = np.gradient(dot_constant, axis=0) / dy
            
ddot_dy_seasonal_monthly = np.nanmean(ddot_dy_seasonal * 1000, axis=2)
ddot_dy_seasonal_mean = np.nanmean(ddot_dy_seasonal_monthly, axis=2)

ddot_dy_seasonal_faster = ddot_dy_seasonal_monthly[:, :, 2:7]
ddot_dy_seasonal_slower_1 = ddot_dy_seasonal_monthly[:, :, :2]
ddot_dy_seasonal_slower_2 = ddot_dy_seasonal_monthly[:, :, 7:]

ddot_dy_seasonal_slower = np.concatenate((ddot_dy_seasonal_slower_1, ddot_dy_seasonal_slower_2), axis=2)

ddot_dy_seasonal_faster_mean = np.nanmean(ddot_dy_seasonal_faster, axis=2)
ddot_dy_seasonal_slower_mean = np.nanmean(ddot_dy_seasonal_slower, axis=2)

ddot_dy_seasonal_anomaly = ddot_dy_seasonal_faster_mean - ddot_dy_seasonal_slower_mean

pl.figure()
pl.clf()
m = Basemap(projection='spstere', boundinglat=-60, lon_0=180, resolution='l')
m.drawmapboundary()
m.drawcoastlines(zorder=10)
m.fillcontinents(zorder=10)
m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
        
grid_lats, grid_lons = np.meshgrid(lat, lon)
stereo_x, stereo_y = m(grid_lons, grid_lats)
        
m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(ddot_dy_seasonal_mean)), cmap='RdBu_r')
c = m.colorbar()
pl.title('$\partial$DOT/$\partial$y')
pl.clim(1, -1)
c.set_label('meridional DOT gradient (mm / km)')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/CoastalCurrent/mean_gradient_map.png', format='png', transparent=True, dpi=300)
pl.close()

pl.figure()
pl.clf()
m = Basemap(projection='spstere', boundinglat=-60, lon_0=180, resolution='l')
m.drawmapboundary()
m.drawcoastlines(zorder=10)
m.fillcontinents(zorder=10)
m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
        
grid_lats, grid_lons = np.meshgrid(lat, lon)
stereo_x, stereo_y = m(grid_lons, grid_lats)
        
m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(ddot_dy_seasonal_anomaly)), cmap='RdBu_r')
c = m.colorbar()
pl.title('$\partial$DOT/$\partial$y anomaly (MAMJJ - ASONDJF)')
pl.clim(0.15, -0.15)
c.set_label('meridional DOT gradient anomaly (mm / km)')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/CoastalCurrent/gradient_anomaly_map.png', format='png', transparent=True, dpi=300)
pl.close()