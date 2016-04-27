from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as pl
from mpl_toolkits.basemap import Basemap
import os

yr = input('What year? (xxxx): ')

os.chdir('/Users/jmh2g09/Documents/PhD/Data/Gridded/' + yr + '/MDT')

for file in os.listdir():
    if file[-6:] == 'MDT.nc':
        year = file[0:4]
        month = file[4:6]
        nc = Dataset(file, 'r')
        lat = nc.variables['Latitude'][:]
        lon = nc.variables['Longitude'][:]
        mdt = nc.variables['mean_dynamic_topography'][:]
        ice_conc = nc.variables['sea_ice_concentration'][:]
        nc.close()
        
        for ilat in range(1, np.shape(mdt)[1]):
            if np.all(np.isnan(mdt[:, ilat])):
                print(ilat)
        
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
        pl.clim(5, -5)
        m.contour(stereo_x, stereo_y, ice_conc, [60,])
        pl.savefig('Figures/'+ year + '_' + month + '_MDT_1degree_stereo.png', 
            format='png')
        pl.close()