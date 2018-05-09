import os
import functions as funct
from netCDF4 import Dataset
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as pl
from mpl_toolkits.basemap import Basemap
import os
from datetime import date

# Load data
ssha = np.full((59, 361, 72), fill_value=np.NaN)
it = 0
for year in ['2011', '2012', '2013', '2014', '2015', '2016']:
    for month in ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']:
        file = '/Users/jmh2g09/Documents/PhD/Data/Gridded/' + year + month + '_grid.nc'
        if os.path.isfile(file):
            nc = Dataset(file, 'r')
            lat = nc.variables['lat'][:]
            lon = nc.variables['lon'][:]
            ssha[:, :, it] = nc.variables['sea_surface_height_anomaly_seasonal_offset'][:]
            nc.close()
            it += 1

A = np.array(range(1, 73))
slopes = np.full((59, 361), fill_value=np.NaN)

for ilat in range(59):
    for ilon in range(361):
        ts = np.array(ssha[ilat, ilon, :])
        ts_nan = np.where(np.isfinite(ts) == 1)[0]
        if np.any(np.isfinite(ts) == 1):
            ts_t = A[ts_nan] / 12
            ts_d = ts[ts_nan] * 1000
        
            slopes[ilat, ilon], intercept, r_value, p_value, std_err = stats.linregress(ts_t, ts_d)


pl.figure()
pl.clf()
m = Basemap(projection='spstere', boundinglat=-50, lon_0=180, resolution='l')
m.drawmapboundary()
m.drawcoastlines(zorder=10)
m.fillcontinents(zorder=10)
m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
grid_lats, grid_lons = np.meshgrid(lat, lon)
stereo_x, stereo_y = m(grid_lons, grid_lats)
m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(slopes)), cmap='RdBu_r')
c = m.colorbar()
#pl.clim(-5, 5)
#pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Gridded/Figures/' + year + '/' + year + month + '_DOT_filt.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
pl.show()
pl.close()