import numpy as np
import os
from netCDF4 import Dataset
import matplotlib.pyplot as pl
from mpl_toolkits.basemap import Basemap

ssh = []
for year in ['2010', '2011', '2012', '2013', '2014', '2015', '2016']:

    os.chdir('/Users/jmh2g09/Documents/PhD/Data/Gridded/SSH/' + year)

    for file in os.listdir():
        if file[-11:] == 'SSH_filt.nc':

            nc = Dataset(file, 'r')
        
            ssh.append(nc.variables['sea_surface_height'][:])
            lat = nc.variables['latitude'][:]
            lon = nc.variables['longitude'][:]
        
            nc.close()

ssh_mean = np.nanmean(ssh, axis=0)

nc = Dataset('/Users/jmh2g09/Documents/PhD/Data/Gridded/SSH/SSH_mean.nc', 'w')

nc.createDimension('lat', np.size(lat))
nc.createDimension('lon', np.size(lon))

latitudes = nc.createVariable('latitude', float, ('lat',))
longitudes = nc.createVariable('longitude', float, ('lon',))
ssh_annual_mean = nc.createVariable('mean_sea_surface_height', float, ('lat','lon'))
latitudes[:] = lat
longitudes[:] = lon
ssh_annual_mean[:] = ssh_mean

nc.close()

pl.figure()
pl.clf()
m = Basemap(projection='spstere', boundinglat=-50, lon_0=180, resolution='l')
m.drawmapboundary()
m.drawcoastlines(zorder=10)
#m.fillcontinents(zorder=10)
m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
grid_lats, grid_lons = np.meshgrid(lat, lon)
stereo_x, stereo_y = m(grid_lons, grid_lats)
m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(ssh_mean)), cmap='RdBu_r')
m.colorbar()
#pl.clim(0, -1.9)
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Gridded/SSH/mean_DOT.png', format='png', transparent=True, dpi=300)
pl.close()