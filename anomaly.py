import numpy as np
import os
from netCDF4 import Dataset
import matplotlib.pyplot as pl
from mpl_toolkits.basemap import Basemap

yr = input('What year? (xxxx): ')

os.chdir('/Users/jmh2g09/Documents/PhD/Data/Gridded/' + yr + '/Geoid')

nc = Dataset('2012_geoid_mean.nc', 'r')

mean_2012 = nc.variables['Mean Sea Surface Height'][:]
lat = nc.variables['Latitude'][:]
lon = nc.variables['Longitude'][:]

nc.close()

for file in os.listdir():
    if file[-8:] == 'geoid.nc':
        year = file[:4]
        month = file[4:6]
        print(file)
        nc = Dataset(file, 'r')
        
        ssh = nc.variables['Sea Surface Height'][:]
        
        nc.close()
        
        ssh_anomaly = ssh - mean_2012
        os.chdir('/Users/jmh2g09/Documents/PhD/Data/Gridded/' + yr + '/Geoid/Anomalies')
        nc = Dataset(year + month + '_geoid_anomaly.nc', 'w', FORMAT='NETCDF4_CLASSIC')
        nc.description = 'Sea Surface Height Anomaly data for ' + month + ', ' + year
        
        nc.createDimension('lat', np.size(lat))
        nc.createDimension('lon', np.size(lon))
        
        latitudes = nc.createVariable('Latitude', float, ('lat',))
        longitudes = nc.createVariable('Longitude', float, ('lon',))
        ssh_anom = nc.createVariable('Sea Surface Height Anomaly', float, ('lon','lat'))
        
        latitudes[:] = lat
        longitudes[:] = lon
        ssh_anom[:] = ssh_anomaly
        
        nc.close()
        os.chdir('/Users/jmh2g09/Documents/PhD/Data/Gridded/' + yr + '/Geoid')
        grid_lats, grid_lons = np.meshgrid(lat, lon)
        
        pl.figure(str(year) + '_' + str(month) + '_ssh_anomaly_above_GOCO05s_1degree_stereo')
        pl.clf()
        m = Basemap(projection='spstere', boundinglat=-50, lon_0=180, resolution='l')
        m.drawmapboundary()
        m.drawcoastlines(zorder=10)
        m.fillcontinents(zorder=10)
        m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
        m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
        stereo_x, stereo_y = m(grid_lons, grid_lats)
        m.pcolor(stereo_x, stereo_y, ssh_anomaly)
        m.colorbar()
        A = np.std(ssh_anomaly)
        B = np.mean(ssh_anomaly)
        pl.clim(B-2*A, B+2*A)
        pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Gridded/' + yr + '/Geoid/Anomalies/Figures/' + str(year) + '_' + str(month) + '_ssh_anomaly_above_GOCO05s_1degree_stereo.png', format='png', transparent=True)
        pl.close()