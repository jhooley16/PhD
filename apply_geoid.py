import os
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as pl
from mpl_toolkits.basemap import Basemap

# Open the Geoid dataset
os.chdir('/Users/jmh2g09/Documents/PhD/Geoid')

nc = Dataset('geoidheight.nc', 'r')

lat_geoid = nc.variables['lat'][:]
#print('Shape of geoid lat is ' + str(np.shape(lat_geoid)))
#print('The min of the geoid lat is ' + str(np.min(lat_geoid)))
#print('The max of the geoid lat is ' + str(np.max(lat_geoid)))

lon_geoid = nc.variables['lon'][:]
# The lon_geoid is from 0:360, the lon_data is from -180:180
lon_geoid = lon_geoid - 180.
#print('Shape of geoid lon is ' + str(np.shape(lon_geoid)))
#print('The min of the geoid lon is ' + str(np.min(lon_geoid)))
#print('The max of the geoid lon is ' + str(np.max(lon_geoid)))

geoid_height = nc.variables['geoid_height'][:]
# The geoid shape needs to be transposed
geoid_height = np.transpose(geoid_height)
#print('Shape of geoid_height is ' + str(np.shape(geoid_height)))

nc.close()

yr = input('What year? (xxxx): ')

# Open the gridded ssh data
os.chdir('/Users/jmh2g09/Documents/PhD/Data/Gridded/' + yr)
for file in os.listdir():
    if file[-3:] == '.nc':
        month = file[-13:-11]
        year = file[-17:-13]
        print(file)

        nc = Dataset(file, 'r')

        lat_data = nc.variables['Latitude'][:]
        #print('The shape of the data lat is ' + str(np.shape(lat_data)))
        #print('The min of the data lat is ' + str(np.min(lat_data)))
        #print('The max of the data lat is ' + str(np.max(lat_data)))

        lon_data = nc.variables['Longitude'][:]
        #print('The shape of the data lon is ' + str(np.shape(lon_data)))
        #print('The min of the data lon is ' + str(np.min(lon_data)))
        #print('The max of the data lon is ' + str(np.max(lon_data)))

        ssh_data = nc.variables['Sea Surface Height'][:]
        #print('The shape of the ssh data is ' + str(np.shape(ssh_data)))

        nc.close()

        # Generate a meshgrid of lat,lon
        grid_lats, grid_lons = np.meshgrid(lat_data, lon_data)



        # Both the geoid height and the altimetry data are taken reference to the
        # WGS84 reference ellipsoid. So the geoid height needs to be taken away 
        # from the altimetry data to get the altimetry data above the geoid.

        ssh_above_geoid = ssh_data - geoid_height
        #print('The shape of the ssh above geoid is ' + str(np.shape(ssh_above_geoid)))

        # Plot the ssh with no geoid corrections
        pl.figure(str(year) + '_' + str(month) + '_ssh_nogeoid_1degree_stereo')
        pl.clf()
        m = Basemap(projection='spstere', boundinglat=-50, lon_0=180, resolution='l')
        m.drawmapboundary()
        m.drawcoastlines(zorder=10)
        m.fillcontinents(zorder=10)
        m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
        m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
        stereo_x, stereo_y = m(grid_lons, grid_lats)
        m.pcolor(stereo_x, stereo_y, ssh_data)
        m.colorbar()
        pl.clim(-50, 50)
        pl.savefig('Figures/'+ str(year) + '_' + str(month) + '_ssh_nogeoid_1degree_stereo.png', 
            format='png', transparent=True)
        pl.close()

        pl.figure(str(year) + '_' + str(month) + '_ssh_above_geoid_1degree_stereo')
        pl.clf()
        m = Basemap(projection='spstere', boundinglat=-50, lon_0=180, resolution='l')
        m.drawmapboundary()
        m.drawcoastlines(zorder=10)
        m.fillcontinents(zorder=10)
        m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
        m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
        stereo_x, stereo_y = m(grid_lons, grid_lats)
        m.pcolor(stereo_x, stereo_y, ssh_above_geoid)
        m.colorbar()
        A = np.std(ssh_above_geoid)
        B = np.mean(ssh_above_geoid)
        pl.clim(B-2*A, B+2*A)
        pl.savefig('Figures/' + str(year) + '_' + str(month) + '_ssh_above_geoid_1degree_stereo.png', 
            format='png', transparent=True)
        pl.close()

# Plot the Geoid
pl.figure('geoid_1degree_stereo')
pl.clf()
m = Basemap(projection='spstere', boundinglat=-50, lon_0=180, resolution='l')
m.drawmapboundary()
m.drawcoastlines(zorder=10)
m.fillcontinents(zorder=10)
m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
stereo_x, stereo_y = m(grid_lons, grid_lats)
m.pcolor(stereo_x, stereo_y, geoid_height)
m.colorbar()
#pl.clim(-150, 150)
pl.savefig('Figures/geoid_1degree_stereo.png', format='png', transparent=True)
pl.close()