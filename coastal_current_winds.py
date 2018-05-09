import os
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as pl
from mpl_toolkits.basemap import Basemap

# Cycle through the files and extract the wind data
u_wind = np.full((59, 361, 7, 12), fill_value=np.NaN)
v_wind = np.full((59, 361, 7, 12), fill_value=np.NaN)
ice = np.full((59, 361, 7, 12), fill_value=np.NaN)
for year in ['2010', '2011', '2012', '2013', '2014', '2015', '2016']:
    os.chdir('/Users/jmh2g09/Documents/PhD/Data/Wind/')
    for month in ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']:
        if os.path.isfile(year + month + '_wind.nc'):
            nc = Dataset(year + month + '_wind.nc', 'r')
            u_wind[:, :, int(year)-2010, int(month)-1] = nc.variables['u_wind'][:]
            v_wind[:, :, int(year)-2010, int(month)-1] = nc.variables['v_wind'][:]
            lat = nc.variables['latitude'][:]
            lon = nc.variables['longitude'][:]
            nc.close()

u_wind_monthly = np.nanmean(u_wind, axis=2)

u_wind_faster = u_wind_monthly[:, :, 4:7]
u_wind_slower_1 = u_wind_monthly[:, :, :2]
u_wind_slower_2 = u_wind_monthly[:, :, 11:]

u_wind_slower = np.concatenate((u_wind_slower_1, u_wind_slower_2), axis=2)

u_wind_faster_mean = np.nanmean(u_wind_faster, axis=2)
u_wind_slower_mean = np.nanmean(u_wind_slower, axis=2)

u_wind_anomaly = u_wind_faster_mean - u_wind_slower_mean

pl.figure()
pl.clf()
m = Basemap(projection='spstere', boundinglat=-60, lon_0=180, resolution='l')
m.drawmapboundary()
m.drawcoastlines(zorder=10)
m.fillcontinents(zorder=10)
m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
grid_lats, grid_lons = np.meshgrid(lat, lon)
stereo_x, stereo_y = m(grid_lons, grid_lats)

m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(u_wind_anomaly)), cmap='RdBu_r')
c = m.colorbar()
#pl.title('10 m u-wind anomaly (AMJ - NDJF)')
pl.clim(3, -3)
c.set_label('m s$^{-1}$')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/CoastalCurrent/Figures/wind_anomaly_map.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
pl.close()