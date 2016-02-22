import numpy as np
import os
from netCDF4 import Dataset
import matplotlib.pyplot as pl
from mpl_toolkits.basemap import Basemap

yr = input('What year? (xxxx): ')

os.chdir('/Users/jmh2g09/Documents/PhD/Data/Gridded/' + yr + '/MDT')

mdt = []
for file in os.listdir():
    if file[-6:] == 'MDT.nc':
        
        nc = Dataset(file, 'r')
        
        mdt.append(nc.variables['mean_dynamic_topography'][:])
        lat = nc.variables['Latitude'][:]
        lon = nc.variables['Longitude'][:]
        
        nc.close()

mdt_mean = np.ma.mean(mdt, axis=0)
mdt_month_mean = np.ma.mean(np.ma.mean(mdt, 2), 1)

nc = Dataset('2012_mdt_mean.nc', 'w', format='NETCDF4_CLASSIC')

nc.createDimension('lat', np.size(lat))
nc.createDimension('lon', np.size(lon))
nc.createDimension('months', np.size(mdt_month_mean))

latitudes = nc.createVariable('Latitude', float, ('lat',))
longitudes = nc.createVariable('Longitude', float, ('lon',))
mdt_annual_mean = nc.createVariable('annual_mean_dynamic_topography', float, ('lon','lat'))
mdt_monthly_mean = nc.createVariable('monthly_mean_dynamic_topography', float, ('months', ))

latitudes[:] = lat
longitudes[:] = lon
mdt_annual_mean[:] = mdt_mean
mdt_monthly_mean[:] = mdt_month_mean

nc.close()

pl.figure()
pl.clf
pl.plot(np.arange(1, 13), mdt_month_mean - np.mean(mdt_month_mean))
pl.xlim([1, 12])
pl.ylabel('Monthly Average MDT (m)')
pl.xlabel('Month')
pl.savefig('Figures/monthly_mean_MDT.png', format='png')
pl.close()

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
m.pcolor(stereo_x, stereo_y, mdt_mean)
m.colorbar()
pl.clim(5, -5)
pl.savefig('Figures/annual_mean_MDT_1degree_stereo.png', format='png')
pl.close()