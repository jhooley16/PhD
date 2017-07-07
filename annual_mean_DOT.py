import numpy as np
import os
from netCDF4 import Dataset
import matplotlib.pyplot as pl
from mpl_toolkits.basemap import Basemap

dot = []
dot_2 = []
dot_3 = []
ice = []
# Initiate the arrays to be filled, separating out the years and the months
dot_array = np.full((59, 361, 12, 7), fill_value=np.NaN)
dot_array_2 = np.full((59, 361, 12, 7), fill_value=np.NaN)
ice_array = np.full((59, 361, 12, 7), fill_value=np.NaN)

# Cycle through the months and arrange the grids by year and month
for year in ['2010', '2011', '2012', '2013', '2014', '2015', '2016']:
    os.chdir('/Users/jmh2g09/Documents/PhD/Data/Gridded/DOT/' + year)
    for month in ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']:
        file = year + month + '_DOT_filt.nc'
        if os.path.isfile(file):
            nc = Dataset(file, 'r')
            dot_array[:, :, int(month)-1, int(year)-2010] = nc.variables['dynamic_ocean_topography_seasonal_offset'][:]
            dot_array_2[:, :, int(month)-1, int(year)-2010] = nc.variables['dynamic_ocean_topography_constant_offset'][:]
            ice_array[:, :, int(month)-1, int(year)-2010] = nc.variables['sea_ice_concentration'][:]
            lat = nc.variables['latitude'][:]
            lon = nc.variables['longitude'][:]
            nc.close()

# Take an average of all years and all months for a total mean
dot_mean = np.nanmean(np.nanmean(dot_array, axis=3), axis=2)
dot_2_mean = np.nanmean(np.nanmean(dot_array_2, axis=3), axis=2)

# Take the mean over only the years, leaving a monthly mean for each
dot_mean_months = np.nanmean(dot_array, axis=3)
dot_mean_months_2 = np.nanmean(dot_array_2, axis=3)
# Do the same for the ice, producing a grid of mean ice concentration
ice_mean_months = np.nanmean(ice_array, axis=3)

# Save the all month mean
nc = Dataset('/Users/jmh2g09/Documents/PhD/Data/Gridded/DOT/Means/MDT_mean.nc', 'w')
nc.createDimension('lat', np.size(lat))
nc.createDimension('lon', np.size(lon))
latitudes = nc.createVariable('latitude', float, ('lat',))
longitudes = nc.createVariable('longitude', float, ('lon',))
dot_annual_mean = nc.createVariable('mean_dynamic_topography_seasonal_offset', float, ('lat','lon'))
dot_2_annual_mean = nc.createVariable('mean_dynamic_topography_constant_offset', float, ('lat','lon'))
latitudes[:] = lat
longitudes[:] = lon
dot_annual_mean[:] = dot_mean
dot_2_annual_mean[:] = dot_2_mean
nc.close()

# Plot
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
m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(dot_mean)), cmap='RdBu_r')
c = m.colorbar()
c.set_label('MDT (m)')
pl.clim(0, -2.25)
m.contour(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(ice_mean_months[:, :, 8])), [20, ], colors='k')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Gridded/DOT/Figures/mean_DOT.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
pl.close()

# Cycle through each month, save the grid and plot
mnths = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']
for it in range(12):
    mnth = mnths[it]
    nc = Dataset('/Users/jmh2g09/Documents/PhD/Data/Gridded/DOT/Means/MDT_mean_' + mnth + '.nc', 'w')
    nc.createDimension('lat', np.size(lat))
    nc.createDimension('lon', np.size(lon))
    latitudes = nc.createVariable('latitude', float, ('lat',))
    longitudes = nc.createVariable('longitude', float, ('lon',))
    dot_annual_mean = nc.createVariable('mean_dynamic_topography_seasonal_offset', float, ('lat','lon'))
    dot_2_annual_mean = nc.createVariable('mean_dynamic_topography_constant_offset', float, ('lat','lon'))
    latitudes[:] = lat
    longitudes[:] = lon
    dot_annual_mean[:] = dot_mean_months[:, :, it]
    dot_2_annual_mean[:] = dot_mean_months_2[:, :, it]
    nc.close()

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
    m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(dot_mean_months[:, :, it])), cmap='RdBu_r')
    c = m.colorbar()
    c.set_label('MDT (m)')
    pl.clim(0, -2.25)
    m.contour(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(ice_mean_months[:, :, it])), [20, ], colors='k')
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Gridded/DOT/Figures/' + mnth + '_mean_DOT.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
    pl.close()
