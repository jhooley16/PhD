from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as pl
from mpl_toolkits.basemap import Basemap
import os

os.chdir('/Users/jmh2g09/Documents/PhD/Data/Gridded')

netfile = Dataset('201201.nc', 'r')
lon = netfile.variables['longitude'][:]
lat = netfile.variables['latitude'][:]
data = netfile.variables['ssh'][:]

pl.figure()
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