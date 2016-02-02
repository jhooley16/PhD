from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as pl
from mpl_toolkits.basemap import Basemap
import os

os.chdir('/Users/jmh2g09/Documents/PhD/Data/Gridded/2012/MDT')

for file in os.listdir():
    if file[-6:] == 'MDT.nc':
        year = file[0:5]
        month = file[4:6]
        nc = Dataset(file, 'r')
        lat = nc.variables['Latitude'][:]
        lon = nc.variables['Longitude'][:]
        mdt = nc.variables['mean_dynamic_topography'][:]
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
        m.pcolor(stereo_x, stereo_y, mdt)
        m.colorbar()
        pl.clim(-5, 5)
        pl.savefig('Figures/'+ year + '_' + month + '_MDT_1degree_stereo.png', 
            format='png', transparent=True)
        pl.close()