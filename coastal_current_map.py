import os
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as pl
from mpl_toolkits.basemap import Basemap

# Calculate the grid spacing in the y-direction (km)
dy = 0.5 * 60 * 1.862
# radius of earth (m)
a = 6378137.

# Cycle through the files and extract the DOT data
ddot_dy_seasonal = np.full((59, 361, 7, 12), fill_value=np.NaN)
ddot_dy_constant = np.full((59, 361, 7, 12), fill_value=np.NaN)
ddot_dx_seasonal = np.full((59, 361, 7, 12), fill_value=np.NaN)
ddot_dx_constant = np.full((59, 361, 7, 12), fill_value=np.NaN)
ice = np.full((59, 361, 7, 12), fill_value=np.NaN)
for year in ['2010', '2011', '2012', '2013', '2014', '2015', '2016']:
    os.chdir('/Users/jmh2g09/Documents/PhD/Data/Gridded/DOT/' + year)
    for month in ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']:
        if os.path.isfile(year + month + '_DOT_filt.nc'):
            nc = Dataset(year + month + '_DOT_filt.nc', 'r')
            dot_seasonal = nc.variables['dynamic_ocean_topography_seasonal_offset'][:]
            dot_constant = nc.variables['dynamic_ocean_topography_constant_offset'][:]
            ice[:, :, int(year)-2010, int(month)-1] = nc.variables['sea_ice_concentration'][:]
            lat = nc.variables['latitude'][:]
            lon = nc.variables['longitude'][:]
            nc.close()
            
            dx = (np.pi / 180) * a * np.cos(lat * np.pi / 180) / 1000
            dx_2d = np.transpose(np.tile(dx, (361, 1)))
            
            # Put constant and seasonal offset data in arrays arranged by year and month
            # Calculate meridional gradient (dSSH/dy)
            ddot_dy_seasonal[:, :, int(year)-2010, int(month)-1] = np.gradient(dot_seasonal, axis=0) / dy
            ddot_dy_constant[:, :, int(year)-2010, int(month)-1] = np.gradient(dot_constant, axis=0) / dy
            ddot_dx_seasonal[:, :, int(year)-2010, int(month)-1] = np.gradient(dot_seasonal, axis=1) / dx_2d
            ddot_dx_constant[:, :, int(year)-2010, int(month)-1] = np.gradient(dot_constant, axis=1) / dx_2d

# Calculate the monthly mean for seasonality
ddot_dy_seasonal_monthly = np.nanmean(ddot_dy_seasonal * 1000, axis=2)
ddot_dx_seasonal_monthly = np.nanmean(ddot_dx_seasonal * 1000, axis=2)
ice_monthly = np.nanmean(ice, axis=2)
# Calculate the overall mean, for mapping
ddot_dy_seasonal_mean = np.nanmean(ddot_dy_seasonal_monthly, axis=2)
ddot_dx_seasonal_mean = np.nanmean(ddot_dx_seasonal_monthly, axis=2)

# Based on the seasonality of the zonal DOT, split between 'faster'
# and 'slower' months
ddot_dy_seasonal_faster = ddot_dy_seasonal_monthly[:, :, 4:7]
ddot_dy_seasonal_slower_1 = ddot_dy_seasonal_monthly[:, :, :2]
ddot_dy_seasonal_slower_2 = ddot_dy_seasonal_monthly[:, :, 11:]

ddot_dx_seasonal_faster = ddot_dx_seasonal_monthly[:, :, 4:7]
ddot_dx_seasonal_slower_1 = ddot_dx_seasonal_monthly[:, :, :2]
ddot_dx_seasonal_slower_2 = ddot_dx_seasonal_monthly[:, :, 11:]

# Concatenate slower months, as they span the end/start of year
ddot_dy_seasonal_slower = np.concatenate((ddot_dy_seasonal_slower_1, ddot_dy_seasonal_slower_2), axis=2)
ddot_dx_seasonal_slower = np.concatenate((ddot_dx_seasonal_slower_1, ddot_dx_seasonal_slower_2), axis=2)

# Take the mean of each
ddot_dy_seasonal_faster_mean = np.nanmean(ddot_dy_seasonal_faster, axis=2)
ddot_dy_seasonal_slower_mean = np.nanmean(ddot_dy_seasonal_slower, axis=2)

ddot_dx_seasonal_faster_mean = np.nanmean(ddot_dx_seasonal_faster, axis=2)
ddot_dx_seasonal_slower_mean = np.nanmean(ddot_dx_seasonal_slower, axis=2)

# Calculate the DOT difference between faster and slower months
ddot_dy_seasonal_anomaly = ddot_dy_seasonal_faster_mean - ddot_dy_seasonal_slower_mean
ddot_dx_seasonal_anomaly = ddot_dx_seasonal_faster_mean - ddot_dx_seasonal_slower_mean


# Plot the mean gradient map, for location of negative gradient, where the 
# Coastal Current is likely to be
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
m.contour(stereo_x, stereo_y, np.transpose(ice_monthly[:, :, 10]), [10,], colors='k')
m.contour(stereo_x, stereo_y, np.transpose(ice_monthly[:, :, 2]), [10,], colors='k')
m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(ddot_dy_seasonal_mean)), cmap='RdBu_r')
c = m.colorbar()
pl.title('$\partial$DOT/$\partial$y')
pl.clim(1, -1)
c.set_label('mm / km')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/CoastalCurrent/mean_y_gradient_map.png', format='png', transparent=True, dpi=300)
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
m.contour(stereo_x, stereo_y, np.transpose(ice_monthly[:, :, 10]), [10,], colors='k')
m.contour(stereo_x, stereo_y, np.transpose(ice_monthly[:, :, 2]), [10,], colors='k')
m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(ddot_dx_seasonal_mean)), cmap='RdBu_r')
c = m.colorbar()
pl.title('$\partial$DOT/$\partial$x')
pl.clim(1, -1)
c.set_label('mm / km')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/CoastalCurrent/mean_x_gradient_map.png', format='png', transparent=True, dpi=300)
pl.close()

# Plot the faster - slower anomaly
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
m.contour(stereo_x, stereo_y, np.transpose(ice_monthly[:, :, 10]), [10,], colors='k')
m.contour(stereo_x, stereo_y, np.transpose(ice_monthly[:, :, 2]), [10,], colors='k')
m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(ddot_dy_seasonal_anomaly)), cmap='RdBu_r')
c = m.colorbar()
pl.title('$\partial$DOT/$\partial$y anomaly (AMJ - NDJF)')
pl.clim(0.25, -0.25)
c.set_label('mm / km')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/CoastalCurrent/gradient_y_anomaly_map.png', format='png', transparent=True, dpi=300)
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
m.contour(stereo_x, stereo_y, np.transpose(ice_monthly[:, :, 10]), [10,], colors='k')
m.contour(stereo_x, stereo_y, np.transpose(ice_monthly[:, :, 2]), [10,], colors='k')
m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(ddot_dx_seasonal_anomaly)), cmap='RdBu_r')
c = m.colorbar()
pl.title('$\partial$DOT/$\partial$x anomaly (AMJ - NDJF)')
pl.clim(0.25, -0.25)
c.set_label('mm / km')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/CoastalCurrent/gradient_x_anomaly_map.png', format='png', transparent=True, dpi=300)
pl.close()