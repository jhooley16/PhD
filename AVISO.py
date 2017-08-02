import os
from netCDF4 import Dataset
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as pl

os.chdir('/Users/jmh2g09/Documents/PhD/Data/AVISO')

for file in os.listdir():
    if file[-3:] == '.nc':
        print(file)
        year = file[-11:-7]
        month = file[-5:-3]
        
        nc = Dataset(file, 'r')
        
        latitude = nc.variables['lat'][:]
        longitude = nc.variables['lon'][:]
        
        ssha = np.array(np.squeeze(nc.variables['sla'][:]))
        nc.close()

        ssha[abs(ssha) > 100] = np.NaN

        nc_test = Dataset('test.nc', 'w')
        nc_test.createDimension('lat', np.size(latitude))
        nc_test.createDimension('lon', np.size(longitude))
        latitudes = nc_test.createVariable('lat', float, ('lat',))
        longitudes = nc_test.createVariable('lon', float, ('lon',))
        ssha_test = nc_test.createVariable('ssha', float, ('lat','lon'))
        latitudes[:] = latitude
        longitudes[:] = longitude
        ssha_test[:] = ssha
        nc_test.close()
        
        os.system('gmt grdfilter test.nc -D4 -Fg1000 -Nr -f0y -f1x -Gtest_filt.nc')

        os.system('gmt grdsample test_filt.nc -Gtest_subsampled.nc -I1.0/0.5 -R0/360/-79/-50 -fg')
        os.system('rm test.nc')
        os.system('rm test_filt.nc')
        
        nc = Dataset('test_subsampled.nc', 'r')
        latitude = nc.variables['lat'][:]
        longitude = nc.variables['lon'][:]
        ssha_sampled = np.squeeze(nc.variables['z'][:])
        nc.close()
        
        print(len(longitude))
        print(np.shape(ssha_sampled))
        
        pl.figure()
        pl.clf()
        m = Basemap(projection='spstere', boundinglat=-50, lon_0=180, resolution='l')
        m.drawmapboundary()
        m.drawcoastlines(zorder=10)
        #m.fillcontinents(zorder=10)
        m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
        m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])

        grid_lats, grid_lons = np.meshgrid(latitude, longitude)
        stereo_x, stereo_y = m(grid_lons, grid_lats)

        m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(ssha_sampled)), cmap='RdBu_r')
        m.colorbar()
        pl.clim([-0.1, 0.1])
        pl.savefig('subsampled/' + year + month + '_AVISO.png', transparent=True, dpi=300, bbox_inches='tight')
        pl.close()
        pause