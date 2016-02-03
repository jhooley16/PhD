import numpy as np
import os
from netCDF4 import Dataset
import matplotlib.pyplot as pl
from mpl_toolkits.basemap import Basemap

yr = input('What year? (xxxx): ')

os.chdir('/Users/jmh2g09/Documents/PhD/Data/Gridded/' + yr + '/MDT')

nc = Dataset('2012_mdt_mean.nc', 'r')
mean_2012 = nc.variables['annual_mean_dynamic_topography'][:]
lat = nc.variables['Latitude'][:]
lon = nc.variables['Longitude'][:]
nc.close()

for file in os.listdir():
    if file[-6:] == 'MDT.nc':
        year = file[:4]
        month = file[4:6]
        print(file)
        nc = Dataset(file, 'r')
        
        mdt = nc.variables['mean_dynamic_topography'][:]
        #ice_data = nc.variables['Sea Ice Concentration'][:]
        
        nc.close()
        
        mdt_anomaly = mdt - mean_2012
        
        nc = Dataset('Anomalies/' + year + month + '_mdt_anomaly.nc', 'w', FORMAT='NETCDF4_CLASSIC')
        
        nc.createDimension('lat', np.size(lat))
        nc.createDimension('lon', np.size(lon))
        
        latitudes = nc.createVariable('Latitude', float, ('lat',))
        longitudes = nc.createVariable('Longitude', float, ('lon',))
        mdt_anom = nc.createVariable('Mean Dynamic Height Anomaly', float, ('lon','lat'))
        latitudes[:] = lat
        longitudes[:] = lon
        mdt_anom[:] = mdt_anomaly
        
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
        m.pcolor(stereo_x, stereo_y, mdt_anomaly)
        m.colorbar()
        A = np.std(mdt_anomaly)
        B = np.mean(mdt_anomaly)
        pl.clim(B-2*A, B+2*A)
        #m.contour(stereo_x, stereo_y, ice_data, colors='k', levels=[70])
        pl.savefig('Anomalies/Figures/' + str(year) + '_' + str(month) + '_ssh_anomaly_above_GOCO05s_1degree_stereo.png', format='png', transparent=True)
        pl.close()