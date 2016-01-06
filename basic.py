# TODO: Change this script to complete a certain task

import os
import decimal
import numpy as np
import matplotlib.pyplot as pl
from mpl_toolkits.basemap import Basemap
import functions as funct
import test

os.chdir('/Users/jmh2g09/Desktop/201209_MERGE')
directory = os.getcwd()

month = directory[-8:-6]
year = directory[-12:-8]

data = funct.month_data(directory)

lat = data['lat']
lon = data['lon']
ssh = data['ssh']
ice_conc = data['ice_conc']

# Plot a Histogram of the ssh data

b = np.arange(np.round(np.min(ssh), -1), np.round(np.max(ssh), -1), 10)
pl.figure(1)
pl.clf()
pl.hist(ssh, b)
#pl.ylim(0, 100)
pl.xlabel('Sea Surface Height (m)')
pl.ylabel('Frequency')
pl.title(str(month) + ' ' + str(year))
pl.show()

# Plot all the data as a scatter on a polar steographic map

pl.figure(str(month) + '_' + str(year) + '_alongtrack_stereo')
pl.clf()
m = Basemap(projection='spstere', boundinglat=-50, lon_0=180, resolution='l')
m.drawmapboundary()
m.drawcoastlines(zorder=10)
m.fillcontinents(zorder=10)
m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
x_map, y_map = m(lon, lat)
m.scatter(x_map, y_map, c=ssh, edgecolor='none')
m.colorbar()
pl.clim(-50, 50)
pl.title(str(month) + ' ' + str(year))
pl.show()

# Plot only the ssh points above a certain threshold

bad_ssh = []
bad_lat = []
bad_lon = []
bad_height = 50
for i in range(np.size(ssh)):
    if np.absolute(ssh[i]) >= bad_height:
        bad_ssh.append(ssh[i])
        bad_lat.append(lat[i])
        bad_lon.append(lon[i])

pl.figure(str(month) + '_' + str(year) + '_alongtrack_stereo_baddata')
pl.clf()
m = Basemap(projection='spstere', boundinglat=-50, lon_0=180, resolution='l')
m.drawmapboundary()
m.drawcoastlines()
m.fillcontinents(zorder=0)
m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
x1, y1 = m(bad_lon, bad_lat)
m.scatter(x1, y1, c=bad_ssh, edgecolor='none')
m.colorbar()
pl.title('CryoSat-2 points where |ssh| > ' + str(bad_height) + ' m')
#pl.clim(-50, 50)
pl.show()

# Plot all the data on a Equidistant Cylindrical map

pl.figure(str(month) + '_' + str(year) + '_alongtrack_cyl')
pl.clf()
m = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=-30,
            llcrnrlon=-180, urcrnrlon=180, resolution='c')
m.drawmapboundary()
m.drawcoastlines(zorder=10)
m.fillcontinents(zorder=10)
m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
x_map_cyl, y_map_cyl = m(lon, lat)
m.scatter(x_map_cyl, y_map_cyl, c=ssh, edgecolor='none')
m.colorbar()
pl.clim(-50, 50)
pl.show()

# Grid the data to a 1-degree grid

data = funct.grid(ssh, lon, lat, 1)
grid_data = data['Grid']
grid_lon = data['Lon']
grid_lat = data['Lat']

pl.figure(str(month) + '_' + str(year) + '_gridded_1degree_cyl')
pl.clf()
m = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=-30,
            llcrnrlon=-180, urcrnrlon=180, resolution='c')
m.drawmapboundary()
m.drawcoastlines(zorder=10)
m.fillcontinents(zorder=10)
m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
m.pcolor(grid_lon, grid_lat, np.transpose(grid_data))
m.colorbar()
#pl.clim(-50, 50)
pl.show()

pl.figure(str(month) + '_' + str(year) + '_gridded_1degree_stereo')
pl.clf()
m = Basemap(projection='spstere', boundinglat=-50, lon_0=180, resolution='l')
m.drawmapboundary()
m.drawcoastlines(zorder=10)
m.fillcontinents(zorder=10)
m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
X_range, Y_range = np.meshgrid(x_range, y_range)
x_map, y_map = m(X_range, Y_range)
m.pcolor(grid_lon, grid_lat, np.transpose(grid_data))
m.colorbar()
pl.clim(-50, 50)
pl.show()