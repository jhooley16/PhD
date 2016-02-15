import numpy as np
import os
from netCDF4 import Dataset
import matplotlib.pyplot as pl
from mpl_toolkits.basemap import Basemap

yr = input('What year? (xxxx): ')

os.chdir('/Users/jmh2g09/Documents/PhD/Data/Gridded/' + yr + '/MDT')

mdt = []
n = 0
for file in os.listdir():
    if file[-6:] == 'MDT.nc':
        
        nc = Dataset(file, 'r')
        
        mdt.append(nc.variables['mean_dynamic_topography'][:])
        lat = nc.variables['Latitude'][:]
        lon = nc.variables['Longitude'][:]
        
        nc.close()
        n = n + 1

mdt_mean = np.ma.mean(mdt, axis=0)

nc = Dataset('2012_mdt_mean.nc', 'w', format='NETCDF4_CLASSIC')

nc.createDimension('lat', np.size(lat))
nc.createDimension('lon', np.size(lon))

latitudes = nc.createVariable('Latitude', float, ('lat',))
longitudes = nc.createVariable('Longitude', float, ('lon',))
mdt_mean = nc.createVariable('annual_mean_dynamic_topography', float, ('lon','lat'))

latitudes[:] = lat
longitudes[:] = lon
mdt2[:] = mdt_mean

nc.close()

grid_lats, grid_lons = np.meshgrid(lat, lon)

pl.figure()
pl.clf()
m = Basemap(projection='spstere', boundinglat=-50, lon_0=180, resolution='l')
m.drawmapboundary()
m.drawcoastlines(zorder=10)
#m.fillcontinents(zorder=10)
m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
        
stereo_x, stereo_y = m(grid_lons, grid_lats)
m.pcolor(stereo_x, stereo_y, mdt_mean)
m.colorbar()
#pl.clim(5, -5)
#m.contour(stereo_x, stereo_y, ice_conc, [70,])
#pl.savefig('Figures/'+ year + '_' + month + '_MDT_1degree_stereo.png', format='png')
pl.show()
pl.close()