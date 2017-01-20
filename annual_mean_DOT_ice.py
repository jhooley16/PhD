import numpy as np
import os
from netCDF4 import Dataset
import matplotlib.pyplot as pl
from mpl_toolkits.basemap import Basemap

dot = []
dot_2 = []
dot_3 = []
ice = []
for year in ['2011', '2012', '2013', '2014', '2015']:
    os.chdir('/Users/jmh2g09/Documents/PhD/Data/Gridded/DOT/' + year)
    for file in os.listdir():
        if file[-11:] == 'DOT_filt.nc':

            nc = Dataset(file, 'r')
            dot.append(nc.variables['dynamic_ocean_topography_seasonal_offset'][:])
            dot_2.append(nc.variables['dynamic_ocean_topography_no_offset'][:])
            dot_3.append(nc.variables['dynamic_ocean_topography_constant_offset'][:])
            if file[-14:-12] == '09':
                ice.append(nc.variables['sea_ice_concentration'][:])
            lat = nc.variables['latitude'][:]
            lon = nc.variables['longitude'][:]
            nc.close()

dot_mean = np.nanmean(dot, axis=0)
dot_2_mean = np.nanmean(dot_2, axis=0)
dot_3_mean = np.nanmean(dot_3, axis=0)
ice_mean = np.nanmean(ice, axis=0)

nc = Dataset('/Users/jmh2g09/Documents/PhD/Data/Gridded/DOT/MDT_mean.nc', 'w')

nc.createDimension('lat', np.size(lat))
nc.createDimension('lon', np.size(lon))

latitudes = nc.createVariable('latitude', float, ('lat',))
longitudes = nc.createVariable('longitude', float, ('lon',))
dot_annual_mean = nc.createVariable('mean_dynamic_topography_seasonal_offset', float, ('lat','lon'))
dot_2_annual_mean = nc.createVariable('mean_dynamic_topography_no_offset', float, ('lat','lon'))
dot_3_annual_mean = nc.createVariable('mean_dynamic_topography_constant_offset', float, ('lat','lon'))
latitudes[:] = lat
longitudes[:] = lon
dot_annual_mean[:] = dot_mean

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
pl.rcParams['contour.negative_linestyle'] = 'solid'
m.contour(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(dot_mean)), 15, colors='k')
m.contourf(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(dot_mean)), 15)
c = m.colorbar()
c.set_label('DOT (m)')
pl.clim(-2.25, 0)
m.contour(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(ice_mean)), [20, ], colors='w', style='dashed')

pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Gridded/DOT/mean_DOT.png', format='png', transparent=True, dpi=300)
pl.close()